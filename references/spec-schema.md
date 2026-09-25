# spec.json 字段说明

`build_workbook.py` 的唯一输入。顶层：

| 字段 | 必填 | 说明 |
|---|---|---|
| `patent` | 是 | 专利号，用作默认输出文件名 |
| `title` | 是 | Excel 主表第 1 行大标题 |
| `subtitle` | 否 | 第 2 行灰色小字，写单位/波长/D 的定义等 |
| `note` | 否 | 第 3 行红字，写全局疑点提示 |
| `vendors` | 否 | 玻璃厂家列表，默认 `["HOYA","OHARA"]`，可选 `CDGM` / `HIKARI`，也可以是自己放进 `assets/glass/` 的任何厂家。顺序即列顺序 |
| `nd_decimals` | 否 | 专利 nd 实际印刷的小数位数。不给则自动判定；当所有 nd 恰好以 0 结尾时自动判定会偏保守，这时手动指定 |
| `tolerance` | 否 | `[nd容差, vd容差]`，完全手动覆盖，优先级最高 |
| `embodiments` | 是 | 实施例数组，见下 |

## embodiment

| 字段 | 必填 | 说明 |
|---|---|---|
| `name` | 是 | 实施例名，出现在主表首列与粘贴表名 |
| `surfaces` | 是 | 面数组，见下 |
| `states` | 否 | 对焦/变倍状态名数组，如 `["广角","中间","望远"]`、`["∞","500mm"]`。给了就会展开可变间隔列 |
| `variable` | 否 | `{"D1": {"广角": 2.78, "望远": 56.78}, ...}`，键与 surface 的 `D` 字符串对应 |
| `extra_columns` | 否 | 额外数值列名数组，如 `["有效径"]`，取值来自 surface 的 `extra` |
| `general` | 否 | `[["项目名", 值], ...]`，进「各种数据」表；各实施例按列并排 |
| `aspheric` | 否 | `[{"surface":"面6","K":…,"A4":…,"A6":…}, ...]`，键名自由，除 `surface` 外都当系数列 |

## surface

| 字段 | 说明 |
|---|---|
| `i` | 面号，字符串或整数。光阑写 `"STO"` |
| `R` | 曲率半径。平面写 `null`（或 `"inf"` / `"无穷大"`），输出为 `INFINITY` |
| `D` | 数值；或可变间隔的变量名字符串（须在 `variable` 里定义） |
| `nd` / `vd` | 该面之后介质的折射率/阿贝数。空气则省略 |
| `lens` | 透镜编号，如 `"L05"` |
| `group` | 组元/位置描述，如 `"G2(正)"`、`"G1F 部分対称前群"` |
| `type` | `"球面"` 或 `"非球面"`。**标成非球面才会取消模压玻璃扣分**，务必标准确 |
| `stop` | `true` 表示光阑行（`i` 写 `"STO"` 时可省） |
| `extra` | `{"有效径": 57.0}`，配合 `extra_columns` |
| `note` | 备注。以 `★` 开头会自动标红底 |

## 完整示例（含可变间隔与非球面）

```json
{
  "patent": "CN 122469501 A",
  "title": "CN 122469501 A 大光圈内合焦式定焦摄影镜头 — 实施例一",
  "subtitle": "单位 mm / 全部球面除注明外",
  "note": "★ 标记行为专利表格印刷疑点，表中保留印刷原值",
  "vendors": ["CDGM", "OHARA"],
  "embodiments": [{
    "name": "实施例一",
    "states": ["∞", "0.1m"],
    "variable": {
      "D1": {"∞": 16.78, "0.1m": 1.5},
      "D2": {"∞": 3.46,  "0.1m": 18.74}
    },
    "extra_columns": ["有效径"],
    "surfaces": [
      {"i": 1, "R": 86.2878, "D": 3.77, "nd": 1.92286, "vd": 20.88,
       "lens": "L01", "group": "G1(负)", "type": "球面",
       "extra": {"有效径": 54.0}},
      {"i": 6, "R": 100.5591, "D": 1.50, "nd": 1.81000, "vd": 40.99,
       "lens": "L04", "group": "G1(负)", "type": "非球面",
       "note": "非球面，系数见表3"},
      {"i": 12, "R": -73.1939, "D": "D1", "group": "G2(正)",
       "note": "★ 印刷把 D1 写在了面11，此处应为可变 D1"},
      {"i": "STO", "R": null, "D": 5.34, "group": "光阑", "stop": true}
    ],
    "general": [["焦距 f", 50.0], ["F No.", 1.2], ["半视场角", 23.8]],
    "aspheric": [
      {"surface": "面6", "K": 20.5880, "A4": -1.9370e-06, "A6": 6.5719e-09,
       "A8": -2.0656e-11, "A10": -5.0323e-14, "A12": 6.9678e-17}
    ]
  }]
}
```

## 输出细节

- Δnd / Δvd 写成 Excel 公式（`=牌号nd - 专利nd`），改了牌号列会自动重算。
- 生成后建议跑一次 `xlsx` skill 的 `scripts/recalc.py` 让公式落缓存值，
  否则用 pandas 或预览器读到的是空值。
- 可变间隔那几列会标红加粗，一眼能看出哪些是随对焦/变倍变化的。

## 2026-09 新增字段

surface 上：

| 字段 | 说明 |
|---|---|
| `glass_offset` | Zemax **Offset 玻璃解**。`{"base":"HIKARI Q-PSKH1S","base_nd":1.59255,"base_vd":67.8569,"base_dpgf":0.01421,"d_nd":0.00039,"d_vd":0.0431}`。由 `lensmath.py` 在「无等效牌号」时自动写入；`make_zmx.py` 见到就铺 `GLAS <base> 4 …` 行，优先于 `glass` / 模型玻璃 |
| `pgf` | 可选。专利印的部分分散比 θgF（Pg,F 本身，不是 ΔPgF）。有它 `lensmath` 挑 Offset 基准时才把 dPgF 纳入打分 |
| `extra["有効径 φi"]` | 该面的净口径**直径**（不是半径）。由 `aptrace.py` 逐面写入，`apcap.py` 可能再收口 |

`zmx` 下：

| 字段 | 说明 |
|---|---|
| `fix_semi_surfaces` | 要把 `DIAM` 固定（flag=1）的面号列表。`aptrace.py` 写全量 → `vignet.py --write` 把「渐晕定义面」并入（不覆盖，定义面另存 `vig_def_surfaces`）→ **`apcap.py` 结束时再重建全量**（`--no-fix-all` 可关）。所以 apcap 必须排在 vignet 之后 |
| `semi_3d` | `vignet.py` 写的逐面 3D 光束最大半径（列表，按面序）。`clearance.py` / `apcap.py` 拿它当「不许切光束」的下界 |
| `fno_patent` | 专利**各種データ**印的各对焦态 F 数，`{结构名: F}`（结构名 = lensmath 生成的 configs 名，如 `{"INF":2.92,"0.50x":4.49,"MFD(1.40x)":6.64}`）。有它时 `vignet.aperture_cfg()` 以它为准、中间结构按 |β| 插值；没有就按物理光阑固定反算。只印 ∞（佳能 US 常见）就别写近距键 |
| `wfno_override` | `{结构名: F}`，逐结构硬指定近轴工作 F 数，优先级最高 |
| `aperture` | `vignet.py --write` 写入：`{type:"paraxial_working_fno", stop_semi_paraxial, wfno_cfg[], stop_semi_cfg[], epd_cfg[], wfno_fixed_stop[], source[], configs[]}`。`make_zmx.py` 据此写 `FNUM <wfno_cfg[0]> 1` + 逐结构 `APER`，`make_seq.py` 写 `FNO` + `ZOO FNO` |
| `vignetting_source` | `vigfit_merge.py` 写回 OpticStudio 真追迹微调后的渐晕时留的标记 |
| `mfd` | **定焦**：产品标称最短撮影距離（像面起算 mm）。专利没印近距态时写它，`lensmath.py` 自动当 `--mfd` 用，凑齐 INF/0.02x/0.06x/MFD 四结构 |
| `field_type` / `max_angle` | `"angle"` = 视场按物方角度（度），FTYP 0 / CODE V `YAN`（>90° 自动 `WID Y`）。半视场 >90° 的鱼眼必须用它；`max_angle` = 专利 ω，按 6 等分出视场 |
| `fields_cfg` | `{结构名: [视场值…]}`（Zemax 顺序由大到小）：逐结构改视场值 → MCE `YFIE` / `ZOO YAN`（或 `ZOO YRI`）。鱼眼近距结构平物面到不了 >90°，收到 89° |
| `focus_vars` | 可选，`["D18","D21"]`：配套评价函数的变量（MCE THIC 状态位 1 / .seq `ZOO THC 0`）。不写时定焦取 focus/focus2 的 key_before，变焦取「同一变焦位置内会变的间隔去掉面序最后一个」 |
| `zoom` | **变焦镜头**。`{"positions":[{"name":"W","label":"Wide","inf":"W-INF","near":["W-MFD"],"f":24.70,"fno":2.91}, …], "betas":[0.06], "include_near":false, "near_d0_printed":165.0}`。`inf`/`near` 是 `embodiments[].states` 里的状态名，`variable` 里**所有可变间隔**都要有这些状态的值。写了 `zoom` 就不要写 `focus`/`focus2`（没有近距态时才用单组 `focus` 兜底）。`lensmath.py` 据此生成 `configs`（每个结构带 `zoom` 键 + 所有可变间隔 + `beta`/`efl`），`vignet.aperture_cfg` 按 `zoom` 分组、各用 `positions[].fno` 定光阑 |
