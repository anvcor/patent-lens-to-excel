# patent-lens-to-excel

光学镜头专利 → Excel / Zemax `.zmx` / CODE V `.seq` 的 Claude skill（个人备份仓库）。

把日本特开（JP-A）、中国 CN-A、美国 US、WO 等专利里的透镜面数据表提取成 Excel，
反查 OHARA / HOYA / CDGM / HIKARI 玻璃牌号，并直接生成带多重结构与对焦位置解的
Zemax `.zmx` 或 CODE V `.seq`；也支持反方向 `.seq → .zmx`。

## 目录

```
SKILL.md              技能正文（流水线、交付卡、已修的坑、读专利判据）
scripts/              17 个脚本
references/           extraction.md / spec-schema.md / zmx-format.md / codev-seq.md
assets/glass/         玻璃目录：OHARA / HOYA / CDGM 的 .AGF + HIKARI.csv + VERSION.txt
```

## 交付流水线（顺序不能变）

```
aptrace → vignet（拿 axial_3d）→ 把偏小的面抬到轴上需求 → 再跑 vignet 出最终渐晕
→ clearance → apcap（干涉收口 + 按光束收紧 + 轴上抬升）→ make_zmx / make_seq
```

- 收紧之后绝不能再跑 `vignet`（口径是渐晕解的输入，会自激塌缩）
- `apcap` 必须排在 `vignet` 之后
- 详细规则与 7 道交付自查卡见 `SKILL.md`

## 脚本一览

| 脚本 | 作用 |
| --- | --- |
| `find_tables.py` / `crop_sheet.py` | 定位并裁切专利里的面数据表页 |
| `glasslib.py` / `hikari_xlsx_to_csv.py` | 玻璃目录加载与 HIKARI 目录转换 |
| `lensmath.py` | 近轴/玻璃匹配/对焦求解（含 Offset 解、双浮动对焦） |
| `aptrace.py` / `figmeas.py` | 从断面图反推逐面口径 |
| `vignet.py` | 逐结构渐晕求解、`axial_3d` 轴上满光瞳需求 |
| `clearance.py` / `covercheck.py` | 边缘间隙检查、传感器盖板判定 |
| `apcap.py` | 口径干涉收口 / 按光束收紧 / 轴上抬升 |
| `layout_check.py` | 叠加图目视校核 |
| `build_workbook.py` | 输出 Excel |
| `make_zmx.py` / `make_seq.py` | 出 Zemax `.zmx` / CODE V `.seq` |
| `seq2zmx.py` | 反方向：CODE V `.seq` → Zemax `.zmx` |

## 重新打包成 `.skill`

`.skill` 就是 zip，顶层必须是同名文件夹 `patent-lens-to-excel/`：

```bash
mkdir -p build/patent-lens-to-excel
cp -r SKILL.md scripts references assets build/patent-lens-to-excel/
cd build && zip -r ../patent-lens-to-excel.skill patent-lens-to-excel
```

## 版本说明

本仓库是 2026-09-13 时点的合并快照：`scripts/` 取本机最新版（含 `lensmath.solve_obj`
对焦求根修复），`SKILL.md` 与 `references/` 取 `.skill` 同步副本最新版（含「出口三」seq2zmx）。

`assets/glass/` 为各厂家公开发布的玻璃目录，版本见 `assets/glass/VERSION.txt`，
版权归各玻璃厂商所有，此处仅作个人备份。
