# CODE V 序列文件 (.seq) 格式 —— 从 Zemax 模型直接生成

参照物：用户把 `JP2025052870A(nikon z3512 e1).zmx` 导入 CODE V 2026 再导出的 `Z3512.seq`。
`scripts/make_seq.py` 就是照它一比一复刻的，生成结果与 CODE V 自己导出的**逐行结构一致**。

## 骨架

```
RDM;LEN       "VERSION: 2026       LENS VERSION: 92       Creation Date:  8-Sep-2026"
TITLE '"INF"'                 ← 基准结构名，注意是 单引号包双引号
FNO   1.23                    ← 所用共轭下的近轴工作 F 数（= Zemax Paraxial Working F/#）；多结构时另写 ZOO FNO 逐位置给值
DIM   M                       ← 毫米
WL    656.3 587.6 546.1 486.1 435.8     ← **降序**
REF   3                       ← 主波长在上面这行里的序号（546.1 排第 3）
WTW   3 22 30 12 3            ← 权重，与 WL 同序
INI   '   '
XRI   0.0 ×6
YRI   21.7 17.36 13.02 8.68 4.34 0.0    ← 实际像高（对应 Zemax FTYP 3）
WTF   1.0 ×6
VUX / VLX / VUY / VLY  各 6 个（基准结构的渐晕）
DOR   1.15 1.05 1.1
SO    0.0 0.1e11              ← 物面：半径 0，厚度 1e10 = 无限远
S     <R> <THI> [玻璃]        ← 每个面一行
  CIR <半口径>                ← CODE V 的 CIR 是实口径，会挡光（等价于 Zemax 固定的 DIAM）
  STO                         ← 只有光阑面写
  ASP / K / CUF / A B C D / E F G H / J
SI    0.0 0.0                 ← 像面
THI   S26 OAL S18..27 45.688  ← 位置解
UID   S0 "…"                  ← 每面一个 GUID，纯标识，可以整段不写，CODE V 会自己补
ZOO   4  …                    ← 多重结构
GO
```

## Zemax → CODE V 对照（都实测验过）

| Zemax | CODE V | 说明 |
|---|---|---|
| `WAVM i λ w` + `PWAV k` | `WL`(降序) + `REF` + `WTW` | 顺序要重排，权重跟着一起排 |
| `FTYP 3` + `YFLN` | `YRI` | 实际像高；`XRI` 全 0、`WTF` 全 1。**顺序与 Zemax 相反：CODE V 第 1 行是轴上、第 6 行才是最大视场**，四行渐晕和 `ZOO … F<i>` 必须一起倒 |
| `VDY` / `VCY` | `VUY = VCY − VDY`<br>`VLY = VCY + VDY` | 实证 INF 视场1：−0.1670 / 0.4315 → 0.5985 / 0.2645 |
| `VCX`（`VDX` 恒 0） | `VUX = VLX = VCX` | |
| `DIAM <sd> 1`（固定） | `CIR <sd>` | |
| `STOP` | `STO` | |
| `TYPE EVENASPH` + `CONI` + `PARM 2..8` | `ASP` + `K` + `CUF 0.0` + `A B C D` / `E F G H` / `J` | A…J 依次是 r⁴ r⁶ r⁸ r¹⁰ r¹² r¹⁴ r¹⁶ r¹⁸ r²⁰；**CODE V 的 K 与 Zemax 的 CONI 同义** |
| `TOLE <vb> <len>`（写在面 va） | `THI S<va> OAL S<vb>..<va+1> <len>` | 面号两边是同一套（光阑也占一号） |
| `MNUM n` + `LTTL` | `ZOO n` + `ZOO TIT` + `TIT Zk "名字"` | **写入用单层引号**：CODE V 导出时写两层（`TIT Z1 ""INF""`、`TITLE '"INF"'`），但它自己读不回来，结构名会全变成 `' '` —— 实测过 |
| `THIC <面> <cfg> <值>` | `ZOO THI S<面> <各结构值>` + `ZOO THC S<面> 100 …` | THC 100 = 该厚度逐结构独立 |
| `FVCY/FVDY/FVCX`（视场号,结构号） | `ZOO VUY F<i>` / `ZOO VLY F<i>` / `ZOO VUX F<i>` / `ZOO VLX F<i>` | CODE V 按**视场**分组、VUY/VLY 交替写 |

## 排版细节

- 续行：行尾 `&`，下一行以**一个空格**开头。CODE V 大约 72 列折行。
- 很小的数写成 `0.19013e-5` 这种「尾数 <1」的形式（`-1.90130e-06` → `-0.19013e-5`）。
- 非球面系数分三组写：`A;B;C;D` / `E;F;G;H` / `J`。**只要写了 E…H 这组，后面一定跟一行 `J`**
  （哪怕是 0）；高次项全为零时就只写 `A…D` 一行。
- 面号：`S1` 起，物面是 `SO`、像面是 `SI`；光阑占一个面号，与 Zemax 完全对齐。

## 玻璃

CODE V 的牌号 = **去掉所有非字母数字 + `_` + 厂家**：
`J-LAK01` → `JLAK01_HIKARI`，`FDS24-W` → `FDS24W_HOYA`，`S-NBH8` → `SNBH8_OHARA`。

> ⚠ **CODE V 没有 Zemax 的 Offset 玻璃解。** 从 Zemax 导入时，Offset 面会退化成**基准目录玻璃**，
> 两个偏移量直接丢掉。实测 JP2025052870A：EFL 因此从 34.4045 变成 34.4369。
> `make_seq.py --glass exact` 把这些面改写成 CODE V 的模型玻璃 `S <R> <THI> <nd> <vd>`，
> EFL 与专利一致，代价是丢掉真实玻璃的色散曲线。两者按用途选，脚本默认 `catalog` 并打印告警。

## 对焦组：文件只负责铺好起点，位置靠优化定

CODE V 侧的习惯是**用评价函数优化**出对焦组位置，不依赖软件的 solve。所以生成的 `.seq` 只要保证：

- 各可变间隔在 `ZOO THI S<n>` 里逐结构给出起始值，且 `ZOO THC S<n> 100 …`（逐结构独立），
  这样直接就能设成优化变量；
- ⚠ **带 OAL 解的那一面不能再写 `ZOO THI S<va>`**：`.seq` 是顺序执行的命令流，
  `ZOO THI` 排在解后面会把它覆盖成定值，表现为「求解不起作用」。CODE V 自己导出时两行都写
  （那只是当前值转储，**不是合法输入**），别照抄。
- `THI S<va> OAL S<vb>..<va+1> <len>` 默认保留 —— 它编码的是「前群与像面固定 ⇒ 三段间隔之和恒定」
  这条物理约束，不妨碍把其余间隔设成变量；`make_seq.py --no-oal` 可以去掉它，让三个间隔全自由。
- CODE V 允许同一面**既有 OAL 解、又有 `ZOO THI` 逐结构值**（它自己导出时也这么写），不冲突。

## 反方向：`.seq` → `.zmx`（`scripts/seq2zmx.py`）

上表整张倒着用即可，四处不是简单互换：

| 方向 | 要点 |
|---|---|
| `VUY/VLY` → `VDY/VCY` | `VCY=(VUY+VLY)/2`、**`VDY=(VLY−VUY)/2`（负号）**。上光瞳切得多 ⇒ 光瞳中心偏下 ⇒ VDY 负 |
| `FNO` / `ZOO FNO` → `FNUM v 1` + 逐结构 `APER` | CODE V 的 FNO = 所用共轭下的近轴工作 F 数（LensSystemSetupRM p.27–29），与 Zemax **Paraxial Working F/#** 同义，原样搬。**别写 `FNUM v 0`**（Image Space F/#，按 ∞ 共轭 EFL/EPD 定义，内对焦近距会把光阑缩掉） |
| `CIR` → `DIAM …1` + `CLAP` | CODE V 的 CIR 挡光，而 Zemax 的 semi-diameter 不挡，必须补 `CLAP 0 <v> 0` |
| 牌号 `XXX_厂家` → 目录名 | 去掉厂家后缀，按「去掉所有非字母数字」与 `.AGF` 的 NM 行比对（`DQK3L_CDGM`→`D-QK3L`）。查不到就退出，别退化成模型玻璃 |

`.seq` 里的 `CCY`/`THC`（优化变量）、`PIK`（拾取解）、`CMP`/`DSX`/`DSY`/`BTX`/`BTY`（公差与补偿器）
不转；`H`(r¹⁸)/`J`(r²⁰) 非零时 EVENASPH 放不下，要改 Extended Asphere。

**CODE V 自带的「导出 Zemax」不可用**：实测 A2628 玻璃全成 `___BLANK 1 0 1.5 0`、
五个结构的物距与对焦间隔全写成同一个值、VDY 整列反号。

## 往返验证

把脚本生成的 `.seq` 在 CODE V 2026 里打开再导出，规范化后比对（忽略 `UID` / `RDM;LEN` /
浮点尾巴 / 续行位置）——**语义零差异**，CODE V 只加了 `UID` 那一段。改 `make_seq.py` 之后都该这么验一次。

## 不要直接生成 `.len`

`.len` 是**版本号绑定的专有二进制**：实测 91,365 字节，头部 `]\x00\x00\x00\\\x00\x00\x00…`，
内容是定长记录 + 小端 IEEE754 双精度数组 + 80 字节定长的结构名字段 +
`THI THI THI THI THC THC THC THC` 这类操作数令牌表 + 36 个 UUID，文件头标着 `LENS VERSION: 92`。
从零写它性价比极低、且换个 CODE V 版本就可能失效。

正路是 `.seq`：纯文本命令流，CODE V 打开即执行，已往返验证语义零差异；`.len` 让 CODE V 自己另存。
`.seq` 表达不了的设置，写成 CODE V 宏命令追加在文件末尾即可（`.seq` 本来就是命令序列）。
