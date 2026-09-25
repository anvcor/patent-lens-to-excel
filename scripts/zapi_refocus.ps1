# zapi_refocus.ps1 —— 用本机 OpticStudio（ZOS-API，无界面）逐结构重新对焦：让轴上离焦 MTF 的峰值落在 0 离焦
#
#   powershell -File zapi_refocus.ps1 -File X_catalog.zmx -Out X.refocus.json [-Freq 50] [-Samp 2] [-Surf 5]
#   python refocus_merge.py spec.final.json X.refocus.json      # 写回 zmx.configs 的对焦间隔，再 make_zmx / make_seq
#
# 用户 2026-09 要求：定焦件做完 zmx 后重新对焦，看轴上点的离焦 MTF（FFT Through Focus，多色，50 lp/mm），峰值要在中心。
# lensmath 给的对焦间隔是**近轴**像面解，残余球差/色差让最佳 MTF 焦点偏几十 µm
# （JP2013-054269A 0.02x：D5 6.5052 时峰值在 −0.035mm，手调到 6.5300 才居中）。
# 「峰值在 0 离焦」⇔「像面上的轴上 MTF 对对焦间隔取极大」，所以直接对 MCE 里的对焦间隔做一维极大化：
#   目标 = 轴上视场 (MTFT+MTFS)/2，多色（wave 0），FFT，频率 -Freq；粗扫 ±0.2mm 找全局峰，再黄金分割细化到 1e-5。
# 之后用 OpticStudio 自己的 FFT Through Focus MTF 分析读回峰位（抛物线插值）当验收。
# 对焦间隔 = MCE 里状态为 Variable 的第一行 THIC（make_zmx 写的 focus.key_before）；-Surf 可强制。
# 双浮动对焦只动第一组（第二组的分工由 lensmath 凸轮给出，一维没法定两组）。
param([Parameter(Mandatory=$true)][string]$File,
      [string]$Out = "",
      [double]$Freq = 50,
      [int]$Samp = 2,                 # MTF 采样 1=32x32 2=64x64 3=128x128（与离焦 MTF 分析的默认 64x64 一致）
      [int]$Surf = 0,
      [double]$Scan = 0.2,
      [string]$ZOS = "E:\ANSYS Inc\v242\Zemax OpticStudio")
$ErrorActionPreference = "Stop"
Add-Type -Path "$ZOS\ZOSAPI_NetHelper.dll"
$null = [ZOSAPI_NetHelper.ZOSAPI_Initializer]::Initialize($ZOS)
Add-Type -Path "$ZOS\ZOSAPI_Interfaces.dll"
Add-Type -Path "$ZOS\ZOSAPI.dll"
$conn = New-Object ZOSAPI.ZOSAPI_Connection
$app = $conn.CreateNewApplication()
if (-not $app.IsValidLicenseForAPI) { "no ZOS-API license"; exit 2 }
$sys = $app.PrimarySystem
$path = (Resolve-Path $File).Path
$null = $sys.LoadFile($path, $false)
$MFE = $sys.MFE; $MCE = $sys.MCE; $FLDS = $sys.SystemData.Fields
$T = [ZOSAPI.Editors.MFE.MeritOperandType]
$ncfg = [math]::Max(1, $MCE.NumberOfConfigurations)

# 轴上视场号（Y=0 且 X=0 的那个）
$axf = 0
for ($k = 1; $k -le $FLDS.NumberOfFields; $k++) { $fd = $FLDS.GetField($k); if ($fd.X -eq 0 -and $fd.Y -eq 0) { $axf = $k } }
if ($axf -eq 0) { "no on-axis field"; exit 3 }

# 对焦间隔所在的 MCE 行：THIC 且某个结构的格子是 Variable
$row = $null
for ($r = 1; $r -le $MCE.NumberOfOperands; $r++) {
  $op = $MCE.GetOperandAt($r)
  if ($op.Type -ne [ZOSAPI.Editors.MCE.MultiConfigOperandType]::THIC) { continue }
  $sn = $op.Param1
  if ($Surf -gt 0) { if ($sn -eq $Surf) { $row = $op; break } else { continue } }
  if ($sn -eq 0) { continue }
  for ($c = 1; $c -le $ncfg; $c++) {
    if ($op.GetOperandCell($c).GetSolveData().Type -eq [ZOSAPI.Editors.SolveType]::Variable) { $row = $op; break }
  }
  if ($row) { break }
}
if (-not $row) { "no focus THIC row in MCE (pass -Surf)"; exit 4 }
$fsurf = $row.Param1
"focus spacing = THIC surface $fsurf   on-axis field = F$axf   freq = $Freq lp/mm   samp = $Samp"

function MtfAt([double]$x, [int]$c) {
  $row.GetOperandCell($c).DoubleValue = $x
  $null = $MCE.SetCurrentConfiguration($c)
  $mt = $MFE.GetOperandValue($T::MTFT, [int]$Samp, [int]0, [double]$axf, [double]$Freq, [double]0, [double]0, [double]0, [double]0)
  $ms = $MFE.GetOperandValue($T::MTFS, [int]$Samp, [int]0, [double]$axf, [double]$Freq, [double]0, [double]0, [double]0, [double]0)
  return 0.5 * ($mt + $ms)
}

# OpticStudio 自己的 FFT 离焦 MTF：返回轴上 (T+S)/2 峰值所在的离焦量（抛物线插值）
function PeakShift([int]$c) {
  $null = $MCE.SetCurrentConfiguration($c)
  try {
    $an = $sys.Analyses.New_Analysis([ZOSAPI.Analysis.AnalysisIDM]::FftThroughFocusMtf)
    $st = $an.GetSettings()
    $st.Frequency = $Freq; $st.DeltaFocus = 0.1; $st.NumberOfSteps = 81
    $null = $st.Field.SetFieldNumber($axf); $null = $st.Wavelength.UseAllWavelengths()
    $null = $an.ApplyAndWaitForCompletion()
    $res = $an.GetResults()
    $ds = $res.GetDataSeries(0)
    $xs = $ds.XData.Data; $ys = $ds.YData.Data
    $n = $xs.Length; $best = 0; $bv = -1.0
    for ($i = 0; $i -lt $n; $i++) { $v = 0.5 * ($ys[$i, 0] + $ys[$i, 1]); if ($v -gt $bv) { $bv = $v; $best = $i } }
    $pk = $xs[$best]
    if ($best -gt 0 -and $best -lt $n - 1) {
      $a = 0.5 * ($ys[($best-1), 0] + $ys[($best-1), 1]); $b = $bv; $d = 0.5 * ($ys[($best+1), 0] + $ys[($best+1), 1])
      $den = $a - 2 * $b + $d
      if ([math]::Abs($den) -gt 1e-12) { $pk = $xs[$best] + 0.5 * ($a - $d) / $den * ($xs[1] - $xs[0]) }
    }
    $null = $an.Close()
    return @($pk, $bv)
  } catch { return @([double]::NaN, [double]::NaN) }
}

$gr = (3 - [math]::Sqrt(5)) / 2
$rows = @()
"cfg  name              D_old       D_new       dD        MTF_old  MTF_new  peak_old(mm)  peak_new(mm)"
for ($c = 1; $c -le $ncfg; $c++) {
  $null = $MCE.SetCurrentConfiguration($c)
  $name = $MCE.GetOperandAt(1).GetOperandCell($c).Value
  $x0 = $row.GetOperandCell($c).DoubleValue
  $pb = PeakShift $c
  $m0 = MtfAt $x0 $c
  # 粗扫：±Scan，41 点，取全局最大
  $bx = $x0; $bm = $m0
  for ($k = -20; $k -le 20; $k++) {
    $x = $x0 + $Scan * $k / 20.0
    if ($x -lt 0.05) { continue }
    $m = MtfAt $x $c
    if ($m -gt $bm) { $bm = $m; $bx = $x }
  }
  # 黄金分割细化 [bx-h, bx+h]
  $h = $Scan / 20.0; $lo = [math]::Max(0.05, $bx - $h); $hi = $bx + $h
  $p = $lo + $gr * ($hi - $lo); $q = $hi - $gr * ($hi - $lo)
  $fp = MtfAt $p $c; $fq = MtfAt $q $c
  for ($it = 0; $it -lt 40 -and ($hi - $lo) -gt 1e-5; $it++) {
    if ($fp -lt $fq) { $lo = $p; $p = $q; $fp = $fq; $q = $hi - $gr * ($hi - $lo); $fq = MtfAt $q $c }
    else { $hi = $q; $q = $p; $fq = $fp; $p = $lo + $gr * ($hi - $lo); $fp = MtfAt $p $c }
  }
  $xn = [math]::Round(0.5 * ($lo + $hi), 4)
  $mn = MtfAt $xn $c
  $pa = PeakShift $c
  "{0,-4} {1,-16} {2,10:F4}  {3,10:F4}  {4,8:+0.0000;-0.0000}  {5,7:F4}  {6,7:F4}  {7,12:F4}  {8,12:F4}" -f $c, $name, $x0, $xn, ($xn - $x0), $m0, $mn, $pb[0], $pa[0]
  $rows += [ordered]@{ config = $c; name = "$name"; surface = $fsurf; old = $x0; new = $xn; mtf_old = $m0; mtf_new = $mn;
                      peak_old = $pb[0]; peak_new = $pa[0] }
}
if ($Out -ne "") {
  $js = [ordered]@{ file = $path; freq = $Freq; samp = $Samp; field = $axf; surface = $fsurf; configs = $rows } | ConvertTo-Json -Depth 5
  [IO.File]::WriteAllText([IO.Path]::GetFullPath((Join-Path (Get-Location) $Out)), $js, (New-Object Text.UTF8Encoding $false))
  "wrote $Out"
}
$app.CloseApplication()
