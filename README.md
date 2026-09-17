# patent-lens-to-excel

光学镜头专利 → Excel / Zemax `.zmx` / CODE V `.seq` 的 Claude skill。

把日本特开（JP-A）、中国 CN-A、美国 US、WO 等专利里的透镜面数据表提取成 Excel，
反查 OHARA / HOYA / CDGM / HIKARI 玻璃牌号，并直接生成带多重结构、对焦位置解、
逐结构孔径（Paraxial Working F/#）与逐结构渐晕的 Zemax `.zmx` 或 CODE V `.seq`；
也支持反方向 `.seq → .zmx`。本机有 OpticStudio 时可经 ZOS-API 无界面加载交付文件做真机体检。

- 怎么**用**：`SKILL.md`
- 怎么**继续开发**（目录约定、环境、回归、发布、已知问题）：`DEVELOPMENT.md`

## 目录

```
SKILL.md              技能正文（流水线、红线、交付卡、已修的坑、读专利判据）
DEVELOPMENT.md        开发说明（不进 .skill 包）
scripts/              18 个 .py + zapi_vigfit.ps1（OpticStudio 真机渐晕/体检）+ p2p_index.json
references/           extraction.md / spec-schema.md / zmx-format.md / codev-seq.md
assets/glass/         兜底玻璃目录：OHARA / HOYA / CDGM 的 .AGF + HIKARI.csv + VERSION.txt
tools/                build_skill.py（打 .skill 包）、push-to-github.ps1
docs/                 示例图
```

交付物（`.zmx` / `.seq` / `.xlsx`）与各专利的 `spec_*.json` 中间件**不放在本仓库**，见 `DEVELOPMENT.md` §3。

## 交付流水线（顺序不能变）

```
which_example → spec.json → lensmath → covercheck → aptrace → apcap --trim 0
→ vignet（孔径模型 + 逐结构渐晕 + 轴上体检）→ [抬口径 → 再跑 vignet]
→ clearance → apcap → make_zmx / make_seq
→ zapi_vigfit（OpticStudio Set Vignetting + 收到 ±1 全过）→ vigfit_merge → 重出 → zapi_vigfit -CheckOnly
```

- 收紧之后绝不能再跑 `vignet`（口径是渐晕解的输入，会自激塌缩）；`apcap` 必须排在 `vignet` 之后
- 孔径类型固定 Paraxial Working F/#（`FNUM v 1` + 逐结构 `APER`），不用 Image Space F/#
- 详细规则与 8 道交付自查卡见 `SKILL.md`

## 脚本一览

| 脚本 | 作用 |
| --- | --- |
| `which_example.py` | 查 PhotonsToPhotos 索引：这篇专利哪个実施例是实物 |
| `find_tables.py` / `crop_sheet.py` | 定位并裁切专利里的面数据表页 |
| `glasslib.py` / `hikari_xlsx_to_csv.py` | 玻璃目录加载与 HIKARI 目录转换 |
| `lensmath.py` | d/e 线判定、玻璃匹配（厂家优先、无铅、Offset 解）、近轴校验、对焦求解（单组/双浮动/链式） |
| `aptrace.py` / `figmeas.py` | 从断面图反推逐面口径 |
| `vignet.py` | 逐结构孔径模型（工作 F 数）、逐结构渐晕（光阑参考坐标）、`axial_3d` 轴上满光瞳需求 |
| `clearance.py` / `covercheck.py` | 边缘间隙检查、传感器盖板判定 |
| `apcap.py` | 口径干涉收口 / 按光束收紧 / 轴上抬升 |
| `layout_check.py` | 叠加图目视校核 |
| `build_workbook.py` | 输出 Excel |
| `make_zmx.py` / `make_seq.py` | 出 Zemax `.zmx` / CODE V `.seq` |
| `seq2zmx.py` | 反方向：CODE V `.seq` → Zemax `.zmx` |
| `zapi_vigfit.ps1` / `vigfit_merge.py` | OpticStudio（ZOS-API）真机渐晕拟合与体检，结果写回 spec |

## 打包成 `.skill`

```bash
python tools/build_skill.py            # → E:\Download\patent-lens-to-excel_YYYYMMDD.skill
```

`.skill` 就是 zip，顶层必须是同名文件夹、路径必须是正斜杠（别用 PowerShell 的 Compress-Archive）。

## 版权

`assets/glass/` 为各厂家公开发布的玻璃目录，版本见 `assets/glass/VERSION.txt`，
版权归各玻璃厂商所有，此处仅作个人备份。
