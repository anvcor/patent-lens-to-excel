# patent-lens-to-excel 开发说明

光学镜头专利（JP / US / CN / WO）→ Excel / Zemax `.zmx` / CODE V `.seq` 的 Claude skill。
这份文档写给**继续开发这个 skill 的人**（包括以后的 Claude 会话）；怎么*用* skill 看 `SKILL.md`。

---

## 1. 一个目录，四个副本的关系

| 副本 | 位置 | 角色 |
|---|---|---|
| **开发目录（唯一源头）** | `E:\Download\patent-lens-to-excel\` | 改代码、改文档都在这里；同时是 git 仓库 |
| GitHub 备份 | https://github.com/anvcor/patent-lens-to-excel （`main`，直推） | 版本历史 |
| claude.ai 服务端 | Settings → Skills | **Claude 实际加载的权威副本**，只能上传 `.skill` 替换（版本化，旧版保留） |
| 本机安装缓存 | `C:\Users\T1791\AppData\Roaming\Claude\local-agent-mode-sessions\skills-plugin\…\skills\patent-lens-to-excel\` | 服务端同步下来的物化副本，**会被同步冲掉，不要在这里改** |

2026-09-17 之前本机有三份（`patent-lens-to-excel` 旧版、`_repo` git 副本、`_新脚本` 平铺编辑副本），
已合并成现在这一个。合并前的完整备份：`E:\Download\patent-lens-to-excel_合并前备份_20260917.zip`。

## 2. 目录结构

```
patent-lens-to-excel/
├── SKILL.md                 skill 正文：流程、红线、交付卡、各类专利的判据、已修的坑
├── DEVELOPMENT.md           本文件（不进 .skill 包）
├── README.md                GitHub 首页（不进 .skill 包）
├── scripts/                 18 个 .py + zapi_vigfit.ps1 + p2p_index.json + spec 示例
├── references/              extraction.md / spec-schema.md / zmx-format.md / codev-seq.md
├── assets/glass/            兜底玻璃目录（OHARA/HOYA/CDGM .AGF + HIKARI.csv），实际优先用用户 Glasscat
├── tools/
│   ├── build_skill.py       打 .skill 包（正斜杠、同名顶层文件夹、排除开发文件）
│   └── push-to-github.ps1   提交并推 GitHub
└── docs/                    示例图等（不进 .skill 包）
```

## 3. 文件放哪（硬约定，用户 2026-09-17 定的）

**开发目录里只放 skill 本身。** 交付物、中间件、临时文件一律不许写进来（`.gitignore` 也挡了 `*.zmx *.seq *.xlsx spec_*.json`）。

| 东西 | 放哪 |
|---|---|
| 交付的 `.zmx` / `.seq` / `.xlsx` | `E:\Download\` |
| 各专利 `spec_*.json` 及 `.matched / .vig / .final` 等中间件、一次性辅助脚本 | `E:\Download\patent_specs\` |
| 回归对照、调试输出 | 会话 scratchpad（`%TEMP%\claude\…\scratchpad`） |
| 打出来的 `.skill` 包 | `E:\Download\`（`build_skill.py` 默认） |

脚本本身不假设输出目录（`-o` / `--write` 必须显式给），所以约定靠调用方遵守。

## 4. 本机环境

| 项 | 实际情况（2026-09-17 核对） |
|---|---|
| Python | `E:\Python\python`（3.11.4）。**每次先 `export PYTHONIOENCODING=utf-8`**，否则 ★/⚠ 撞 GBK 报错 |
| Python 包 | `pypdf 6.14.2`、`Pillow 12.3.0`；**`openpyxl` 没装** → `build_workbook.py` 跑前 `E:\Python\python -m pip install openpyxl` |
| poppler | `pdftoppm` / `pdfinfo`（MiKTeX 自带）、`pdftotext`（mingw64）；**路径带空格/括号的 PDF 先复制成短名** |
| 玻璃目录 | `C:\Users\T1791\Documents\Zemax\Glasscat\`：挑各厂最新一本改名成 `<厂家>.AGF` 放进一个临时目录，`PATENT_GLASS_DIR` 指过去（**必须写 Windows 路径**）；OHARA 用 `OHARA_260701` |
| Zemax OpticStudio | 2024 R2，`E:\ANSYS Inc\v242\Zemax OpticStudio`；**ZOS-API 许可可用，可无界面加载 .zmx 验证**（`scripts/zapi_vigfit.ps1`） |
| CODE V | 2026，`E:\CODEV2026`（手册 `doc\*.pdf`、宏 `macro\`）；**COM 接口 `CodeV.Command.2026` 报无许可**（27000@127.0.0.1）→ .seq 只能 `seq2zmx.py` 回转后在 Zemax 里验 |
| git | 全局没配 user.name/email，提交时 `git -c user.name=T1791 -c user.email=t1791669914@gmail.com commit …`；凭证在 Windows 凭据管理器 |

PowerShell 写脚本的三个坑（`zapi_vigfit.ps1` 都踩过）：
含中文的 `.ps1` 必须 **UTF-8 带 BOM**（PS 5.1 否则按 GBK 读、`param()` 直接语法错）；
**变量名不分大小写**（`$F` 与 `$f`、`$n` 与 `$N` 会互相覆盖）；
.NET 方法的返回值不 `$null =` 吞掉会混进函数返回值。用 Python 改 `.ps1` 时注意它是 CRLF。

## 5. 流水线与脚本分工

```
which_example.py            哪个実施例是实物（PhotonsToPhotos 索引，本地缓存 p2p_index.json）
find_tables.py / crop_sheet.py / pypdf    定位、裁图或直读表格 → 手写 spec.json（E:\Download\patent_specs\）
lensmath.py --write         d/e 线判定、玻璃匹配（厂家优先硬约束、无铅、Offset 解）、近轴校验、对焦位置解 → configs
covercheck.py               盖板判定（佳能：专利没有就不补）
aptrace.py / figmeas.py     断面图逐面量口径（专利没印有効径时）
apcap.py --trim 0           先按几何干涉收口
vignet.py --write           孔径模型（逐结构工作 F 数）+ 逐结构渐晕（光阑参考坐标）+ 轴上满光瞳体检
  └ 有「★★ 切到轴上光瞳」→ apcap 抬口径 → 再跑一次 vignet
clearance.py                边厚/间隙体检
dropneg.py                  删负厚度的无光焦度面（遮光平面），并入前间隔、其后改号（mkspec 后、lensmath 前）
apcap.py                    干涉收口 + 轴上抬升 + 重建 fix_semi_surfaces（必须排在 vignet 之后）
build_workbook.py / make_zmx.py / make_seq.py     Excel / .zmx / .seq（输出到 E:\Download\）
zapi_vigfit.ps1             OpticStudio：Set Vignetting 起步 → 取整 → 收到 Py/Px=±1 全过 → JSON
vigfit_merge.py             写回 zmx.vignetting_cfg → 重出 make_zmx / make_seq
zapi_vigfit.ps1 -CheckOnly  交付前真机体检
seq2zmx.py                  反方向 .seq → .zmx（也用来回转验证 make_seq 的输出）
```

脚本之间的依赖：`make_seq` → `make_zmx.wfno_cfg`；`make_zmx` / `seq2zmx` → `vignet.aperture_cfg` / `make_zmx.cfg_paraxial`；
`aptrace` → `figmeas`。都按「同目录 import」，**用 `python scripts/xxx.py` 或 cwd=scripts 运行**。

## 6. 改代码时别破坏的约定

细节和来龙去脉都在 `SKILL.md`，这里只列最容易被改坏的：

1. **孔径类型固定 Paraxial Working F/#**：`.zmx` 头 `FNUM <结构1> 1` + 逐结构 `APER`；`.seq` 写 `FNO` + `ZOO FNO`。
   逐结构值来自 `vignet.aperture_cfg()`：专利印了近距各态 F 数（`zmx.fno_patent`）优先，否则物理光阑固定；
   光阑按 **∞ 结构**定。CODE V 的 FNO 是物方 sin 定义：`FNO_CV² = PWFN² + (β/2)²`。
   **不许回到 `FNUM v 0`（Image Space F/#）**——内对焦镜头近距会把光阑缩掉。
2. **近距结构的几何**：任何按结构重建系统的地方都要用 `vignet.cfg_dmap()`（补齐 key_after / key_last），
   不能只覆盖 key_before。
3. **渐晕的光瞳坐标是光阑参考的**（对齐 Zemax RAIM Real）；拿 OpticStudio 对照时 .zmx 口径必须全部固定。
   光阑半径**不是近轴半径**：R = 主波长轴上实光线瞄近轴入瞳边缘（EPD/2）时在光阑上的高度，
   (Px,Py) 对所有视场、所有波长都线性落成「主光线 + (Px,Py)·R」（2026-09-25 ZOS-API 实测，见 SKILL.md「Real 瞄准的光阑半径」）。
   `vignet.py` 追迹用 .zmx 的**主波长**（`waves_of` 与 `make_zmx` 同一约定），`axial_3d` 取全部系统波长的包络。
4. **顺序**：vignet → clearance → apcap → make_*；`vignet --write` 把渐晕定义面**并入** `fix_semi_surfaces`（原有全量保留，定义面另存 `vig_def_surfaces`），
   apcap 仍会按全部有口径的面重建一次；`--trim` 收紧之后绝不能再跑 vignet（自激塌缩）。
5. **CODE V 视场顺序与 Zemax 相反**（轴上在前），`YRI` / 渐晕 / `ZOO … F<i>` 一起倒；带 OAL 解的面不写 `ZOO THI`。
6. **玻璃**：厂家优先是硬约束；佳能 OHARA 第一；OHARA 只用 `S-`/`L-`；牌号必须在用户 Glasscat 里存在。
7. 生成器读口径用的是 `surfaces[].extra['有効径 φi']`；`zmx.semi_diameters` 只给 clearance / layout_check 用。
8. **交付约定**：定焦 INF/0.02x/0.06x/MFD 四结构、变焦 W∞/M∞/T∞/W0.06x/M0.06x/T0.06x 六结构；默认只出 `_catalog.zmx` + `.seq`。
9. **配套评价函数**（`make_zmx.merit_contrast`）必须与用户样板 `E:\Download\WO2024214585A1_Ex02_NikonZ24-70mmF28SII_catalog_OPT.zmx`
   的评价函数段逐行一致；对焦变量走 `make_zmx.focus_vars()`，.zmx（THIC 状态位 1）与 .seq（`ZOO THC 0`）必须同一组。
10. **`.zmx` 里 `CONF` 行是评价函数操作数**，MCE 之后不许再写 `CONF 1`。
11. **变焦**：结构里写全所有可变间隔，生成器不写 TOLE/OAL；`aperture_cfg` 按 `zoom` 分组；apcap/clearance 的干涉按全行程最小间隔。

## 7. 回归与验证

**改完脚本至少做这三件：**

```bash
export PYTHONIOENCODING=utf-8
cd /e/Download/patent-lens-to-excel
for f in scripts/*.py tools/*.py; do python -c "import ast,io,sys;ast.parse(io.open(sys.argv[1],encoding='utf-8').read())" "$f" || echo "SYNTAX $f"; done
```

1. **回归 spec**（都在 `E:\Download\patent_specs\`，先复制到 scratchpad 再跑，别覆盖）：

   | spec | 覆盖的分支 |
   |---|---|
   | `spec_JP2021-047297A_Ex01_CanonRF100Macro.final.json` | 双浮动对焦、专利近距 F 数（光阑近距收缩）、高倍率微距、8 结构 |
   | `spec_EP4215968A1_Ex04_CanonRF800.final.json` | 光阑在对焦组前、`zmx.epd` |
   | `spec_JP2020-173350A_Ex01_CanonRF600F11.final.json` | DOE 衍射面 |
   | `spec_US20240302626A1_Ex02_CanonRF35.json` | 双浮动 F1.46，入瞳球差大（光阑参考坐标的试金石） |
   | `spec_US20250251576A1_Ex01_CanonRF50F14.final.json` | 内对焦、近距轴上切光 |
   | `spec_JP2022-85382A_Ex03_CanonRF16.json` | 超广角 53°、光瞳像差极大 |
   | `spec_JP2023-140823A_Ex06.final.json` | 整组繰り出し |
   | `spec_JP2024-169698A_Ex02_CanonRF1200.final.json` | 萤石、只印 ∞ |
   | `spec_WO2024214585A1_Ex02_NikonZ2470f28SII.final.json` | **变焦**（W/M/T × ∞/0.06x 六结构）、链式两组对焦、BF 随变焦变、按变焦位置分组的光阑、评价函数 + 对焦变量 |

   看点：`make_zmx` 自校验 EFL 对得上专利 f、「孔径类型 Paraxial Working F/#」「与孔径模型逐结构一致 ✓」；
   `vignet` 没有 `★Py/Px±1被挡`、轴上体检「全部通过」。
2. **OpticStudio 真机**：`powershell -File scripts\zapi_vigfit.ps1 -File X_catalog.zmx -CheckOnly`
   → 逐结构 PWFN = APER、TOTR 恒定、PMAG = 设计倍率、`all +-1 rays pass`。
3. **.seq 回转**：`python scripts/seq2zmx.py X.seq -o <scratch>/rt.zmx --glassdir <gc> --gcat "OHARA=OHARA_260701,HOYA=HOYA20260707" --reverse-fields`
   → 对 `rt.zmx` 再跑一次 `-CheckOnly`，结果应与原 .zmx 逐位相同。

改动较大时，用多 agent 做一轮**对抗式审查**（2026-09-17 那轮抓到 10 个真 bug：浮点死循环、反方向外推、
∞ 锚点、缓存崩溃、DOE 自校验漏项……），每条发现都要求给出复现脚本。

## 8. 发布

```bash
cd /e/Download/patent-lens-to-excel
git add -A && git -c user.name=T1791 -c user.email=t1791669914@gmail.com commit -m "…"
git push origin main                       # 或 tools\push-to-github.ps1 "提交说明"
python tools/build_skill.py                # → E:\Download\patent-lens-to-excel_YYYYMMDD.skill
```

然后 claude.ai → Settings → Skills → 右上角 **Add ▾ → Upload skill** → 选包 → Save → **Upload and replace**
（技能行的 ⋮ 菜单里没有替换，别去那找；是版本化替换，旧版保留）。上传后本机安装缓存会被同步更新。

## 9. 已知问题 / TODO

- **2026-09-17 之前交付的多结构文件都要重出**：旧 `.zmx` 写的是 Image Space F/#（内对焦近距光阑被缩，
  如 RF600 F11 旧文件 MFD 实为 F17.8），旧 `.seq` 把 ∞ 的 F 数铺满 `ZOO FNO`（近距入瞳被放大），
  且旧 `vignet.py` 近距结构几何错（key_after 没更新）。回归时发现 RF50 F1.4 VCM、RF35 F1.4 VCM 在正确几何下
  MFD 结构面11/12 切轴上光束，重出时要抬口径。
- `.seq` 无法在 CODE V 里直接验（COM 无许可），目前靠 `seq2zmx` 回转 + OpticStudio。
- 没有 OpticStudio 的机器上，`vignet.py` 的 VCX 与 Set Vignetting 还差 ~0.01（VDY/VCY 已到 ~0.003）。
- 专利近距 F 数的中间结构插值（按 |β| 拉格朗日）只在 RF100 上用说明书验过（1.00x 插 5.79 vs 说明书 5.7）。
- `make_zmx.py` 的「zmx.aperture 已过期」告警在一次运行里会重复打印几遍（目录版/模型版/自校验各调一次）。
- `build_workbook.py` 依赖的 `openpyxl` 本机 Python 没装。
- 旁支：`E:\Download\Bokeh_Simulation\bokehsim\parsers.py` 读 `FNUM` 忽略第二个字段（单独任务在处理）。

## 10. 变更记录

- **2026-09-25**：`lensmath` Offset 基准重写为 `offset_base()`（按厂家顺序 / 球面件不用模压料 / 窗口外色散优先 / 停产无铅可用），
  `best()` 加 `mold_ok`（按元件排除模压料）；新增可选面字段 `pgf`。TS-E24 II 面7/18/19 → S-BSM4 / S-BAL3 / S-LAH60V。
  回归 10 个 spec + 另 9 个带 Offset 解的历史 spec（按各自真实厂家顺序）旧/新逐面比对：只有无等效面的 Offset 基准变（EFL 逐位不变），基准始终留在首选厂家内。

- **2026-09-17（晚）**：**变焦镜头**支持（`lensmath.zoom_configs`、`aperture_cfg` 按变焦位置分组、make_zmx/make_seq 全可变间隔进 MCE/ZOO、
  covercheck 变量 BF、apcap/clearance 全行程最小间隔）；**配套评价函数 + 对焦变量**（Contrast s+t 80lp/mm GQ3×6 逐结构，
  与用户样板逐行一致；.seq `ZOO THC 0`）；make_zmx 默认只出 catalog；删掉 MCE 后多余的 `CONF 1`；`zmx.mfd`。
  首个变焦回归件 WO2024214585A1 Ex2（Nikon Z 24-70/2.8 S II）。
  对抗审查（8 agent，10 个定焦回归 spec 新旧输出逐字节比对）另修：covercheck 守恒和未减（离焦 2mm）/整组前伸被当恒定 BF、
  make_seq 链式 OAL 丢失与 --no-oal 漏 focus2 key_after、apcap/clearance 最小间隔漏 key_after/key_last、include_near 重名。

- **2026-09-17**：孔径改 Paraxial Working F/#（`aperture_cfg`，专利近距 F 数优先）；`vignet` 近距几何修正
  （`cfg_dmap`）、光阑参考光瞳坐标；`make_seq` 写逐结构 `ZOO FNO`（sin 定义换算）、`seq2zmx` 对应反算；
  新增 `zapi_vigfit.ps1` / `vigfit_merge.py`（OpticStudio 真机渐晕与体检）；`lensmath` 近距物距求根括号、
  专利中间态收进 configs；对抗审查修 10 个 bug；三个本机副本合并为本仓库，新增 `tools/`、本文件。
- 更早的修补记录见 `SKILL.md` 的「已修的坑」表与 git 历史。
