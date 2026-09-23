# zapi_vigfit.ps1 —— 用本机 OpticStudio（ZOS-API，无界面）当判官，把逐结构渐晕微调到 Py/Px=±1 真过得去
#
#   powershell -File zapi_vigfit.ps1 -File X_catalog.zmx -Out X.vigfit.json [-Step 0.01] [-MaxIter 150]
#   python vigfit_merge.py spec.final.json X.vigfit.json      # 写回 zmx.vignetting_cfg，再 make_zmx / make_seq
#
# 为什么需要：vignet.py 瞄的是**近轴入瞳**，Zemax 开 RAIM Real 瞄的是**真实光阑**，大视场下两者有光瞳像差，
# 前几片（面1~3）上的落点差零点几毫米。实测 JP2021-047297A：脚本自检 Py/Px±1 全 OK，
# OpticStudio 里 INF~0.25x 结构视场1~4 的 Py−1 / Px±1 仍被面1/2/3 切；连 OpticStudio 自己的
# Set Vignetting 也是「正好压在口径边上」，Px±1 照样判渐晕。只有真追迹才看得出来。
#
# 做法（只收，不放）：逐结构、逐视场追 Py=±1、Px=±1 四条实光线；
#   Py+1 被挡 → 上沿下移 step：VCY += step/2, VDY -= step/2（下沿不动）
#   Py−1 被挡 → 下沿上移 step：VCY += step/2, VDY += step/2（上沿不动）
#   Px±1 被挡 → VCX += step
# 直到四条都过。Zemax 的约定：py' = VDY + py·(1−VCY)，px' = VDX + px·(1−VCX)。
# 另外报告孔径类型、逐结构 PWFN/WFNO/EPD/PMAG/TOTR，当交付前的真机体检。
param([Parameter(Mandatory=$true)][string]$File,
      [string]$Out = "",
      # 默认原为 0.004 × 40（最多收 0.16）：广角大视场 Px±1 收不住（WO2023181666A1 24-70 GM II 广角端 F1 卡在面28/29/30）。
      # 改成 0.01 × 150：一步最多多收 0.01 光瞳，精度损失可忽略，能收 1.5。
      [double]$Step = 0.01,
      [int]$MaxIter = 150,
      [string]$ZOS = "E:\ANSYS Inc\v242\Zemax OpticStudio",
      [switch]$CheckOnly,
      [switch]$NoSetVig)
$ErrorActionPreference = "Stop"
Add-Type -Path "$ZOS\ZOSAPI_NetHelper.dll"
$null = [ZOSAPI_NetHelper.ZOSAPI_Initializer]::Initialize($ZOS)
Add-Type -Path "$ZOS\ZOSAPI_Interfaces.dll"
Add-Type -Path "$ZOS\ZOSAPI.dll"
$conn = New-Object ZOSAPI.ZOSAPI_Connection
$app = $conn.CreateNewApplication()
if (-not $app.IsValidLicenseForAPI) { "no ZOS-API license"; exit 2 }
$sys = $app.PrimarySystem
$null = $sys.LoadFile((Resolve-Path $File).Path, $false)
# ★ PowerShell 变量名不分大小写：$FLD 与循环变量 $fi 千万别起成 $F / $f（会互相覆盖）
$MFE = $sys.MFE; $MCE = $sys.MCE; $LDE = $sys.LDE; $FLD = $sys.SystemData.Fields
$T = [ZOSAPI.Editors.MFE.MeritOperandType]
$nfld = $FLD.NumberOfFields; $img = $LDE.NumberOfSurfaces - 1
$ncfg = [math]::Max(1, $MCE.NumberOfConfigurations)
$pw = $sys.SystemData.Wavelengths.GetWavelength(1).IsPrimary
for ($w = 1; $w -le $sys.SystemData.Wavelengths.NumberOfWavelengths; $w++) {
  if ($sys.SystemData.Wavelengths.GetWavelength($w).IsPrimary) { $pw = $w } }
"aperture=" + $sys.SystemData.Aperture.ApertureType + " value=" + $sys.SystemData.Aperture.ApertureValue + " rayAiming=" + $sys.SystemData.RayAiming.RayAiming + " primaryWave=$pw"

function Blocked4($fi) {
  $maxY = [math]::Abs($FLD.GetField(1).Y); if ($maxY -eq 0) { $maxY = 1 }
  $hy = $FLD.GetField($fi).Y / $maxY
  $rt = $sys.Tools.OpenBatchRayTrace()
  $bt = $rt.CreateNormUnpol(4, [ZOSAPI.Tools.RayTrace.RaysType]::Real, $img)
  foreach ($pp in @(@(0,1),@(0,-1),@(1,0),@(-1,0))) { $null = $bt.AddRay($pw, 0, $hy, $pp[0], $pp[1], [ZOSAPI.Tools.RayTrace.OPDMode]::None) }
  $null = $rt.RunAndWaitForCompletion(); $null = $bt.StartReadingResults()
  $res = @{}
  foreach ($lab in @("Py+1","Py-1","Px+1","Px-1")) {
    $rn=0; $ec=0; $vc=0; $X=0.0; $Y=0.0; $Z=0.0; $L=0.0; $M=0.0; $NN=0.0; $l2=0.0; $m2=0.0; $n2=0.0; $opd=0.0; $inten=0.0
    $okr = $bt.ReadNextResult([ref]$rn, [ref]$ec, [ref]$vc, [ref]$X, [ref]$Y, [ref]$Z, [ref]$L, [ref]$M, [ref]$NN, [ref]$l2, [ref]$m2, [ref]$n2, [ref]$opd, [ref]$inten)
    if (-not $okr -or $ec -ne 0 -or $vc -ne 0) { $res[$lab] = $vc }
  }
  $null = $rt.Close()     # Close() 有返回值，不吞掉就会混进函数输出，$b 变成数组
  return $res
}

$all = @()
# 先把每个结构的原始渐晕拍个快照：make_zmx 对「所有结构都是 0」的视场不写 MCE 行，
# 那些视场的 VDY/VCY/VCX 是全局量 —— 在结构 c 里改了会一路漏到后面所有结构（只收不放，还回不来）。
# 所以每进一个结构先按快照复位，再判再改。
$snap = @{}
for ($c = 1; $c -le $ncfg; $c++) {
  if ($MCE.NumberOfConfigurations -ge 1) { $null = $MCE.SetCurrentConfiguration($c) }
  $snap[$c] = @(); for ($fi = 1; $fi -le $nfld; $fi++) { $fd = $FLD.GetField($fi); $snap[$c] += ,@($fd.VDX, $fd.VDY, $fd.VCX, $fd.VCY) }
}
"cfg  APER     PWFN     WFNO     EPD      PMAG      TOTR      adjustments / residual"
for ($c = 1; $c -le $ncfg; $c++) {
  if ($MCE.NumberOfConfigurations -ge 1) { $null = $MCE.SetCurrentConfiguration($c) }
  for ($fi = 1; $fi -le $nfld; $fi++) {
    $fd = $FLD.GetField($fi); $v0 = $snap[$c][$fi-1]
    $fd.VDX = $v0[0]; $fd.VDY = $v0[1]; $fd.VCX = $v0[2]; $fd.VCY = $v0[3]
  }
  $cfgv = @(); $log = @()
  # 默认先用 OpticStudio 自己的 Set Vignetting 定边界（按真实光阑瞄准、用全部固定口径，正好压在边上），
  # 再往里收到 ±1 四条都过 —— 结果与 vignet.py 的近似无关，既不切多也不切少。
  # -NoSetVig：从文件里的渐晕（vignet.py 解的）出发只收不放（想对照脚本模型时用）。
  if (-not $CheckOnly -and -not $NoSetVig) { $FLD.SetVignetting(); $log += "SetVig" }
  for ($fi = 1; $fi -le $nfld; $fi++) {
    $fd = $FLD.GetField($fi)
    $it = 0
    while ($true) {
      if (-not $CheckOnly) {
        # 先按「写进 .zmx 以后的精度」取整再判：Set Vignetting 的解正好压在边上，
        # 写文件时舍入 1e-5 就可能又把边沿推回口径上（实测 RF100 0.06x 视场5 Py+1 被面23 切）。
        # VCY/VCX 向上取整（只收），VDY 四舍五入 —— 取整误差由下面的循环兜住。
        $fd.VCY = [math]::Ceiling([math]::Round($fd.VCY * 1e4, 6)) / 1e4
        $fd.VCX = [math]::Ceiling([math]::Round($fd.VCX * 1e4, 6)) / 1e4
        $fd.VDY = [math]::Round($fd.VDY, 4)
      }
      $b = Blocked4 $fi
      if ($b.Count -eq 0 -or $CheckOnly -or $it -ge $MaxIter) { break }
      if ($b.ContainsKey("Py+1")) { $fd.VCY = $fd.VCY + $Step/2; $fd.VDY = $fd.VDY - $Step/2 }
      if ($b.ContainsKey("Py-1")) { $fd.VCY = $fd.VCY + $Step/2; $fd.VDY = $fd.VDY + $Step/2 }
      if ($b.ContainsKey("Px+1") -or $b.ContainsKey("Px-1")) { $fd.VCX = $fd.VCX + $Step }
      $it++
    }
    if ($it -gt 0) { $log += "F$fi x$it" }
    if ($b.Count -gt 0) { $log += ("F{0} STILL {1}" -f $fi, (($b.Keys | ForEach-Object { "$_@S" + $b[$_] }) -join ",")) }
    $cfgv += ,@([math]::Round($fd.VDX,4), [math]::Round($fd.VDY,4), [math]::Round($fd.VCX,4), [math]::Round($fd.VCY,4))
  }
  $all += ,$cfgv
  $av = $sys.SystemData.Aperture.ApertureValue
  $pb = $MFE.GetOperandValue($T::PARB, $img, $pw, 0,0,0,1,0,0)
  $pc = $MFE.GetOperandValue($T::PARC, $img, $pw, 0,0,0,1,0,0)
  $pwfn = [math]::Abs($pc / (2.0 * $pb))
  $wfno = $MFE.GetOperandValue($T::WFNO, 0, 0, 0,0,0,0,0,0)
  $epd  = $MFE.GetOperandValue($T::EPDI, 0, 0, 0,0,0,0,0,0)
  $pmag = $MFE.GetOperandValue($T::PMAG, 0, $pw, 0,0,0,0,0,0)
  $totr = $MFE.GetOperandValue($T::TOTR, 0, 0, 0,0,0,0,0,0)
  "{0,-4} {1,-8:F4} {2,-8:F4} {3,-8:F4} {4,-8:F4} {5,-9:F5} {6,-9:F3} {7}" -f $c,$av,$pwfn,$wfno,$epd,$pmag,$totr,($(if ($log.Count) { $log -join "  " } else { "all +-1 rays pass" }))
}
if ($Out) {
  $obj = @{ file = (Resolve-Path $File).Path; step = $Step; fields = $nfld; configs = $ncfg; vig = $all }
  ($obj | ConvertTo-Json -Depth 6) | Out-File -Encoding utf8 $Out
  "wrote $Out"
}
$app.CloseApplication()
