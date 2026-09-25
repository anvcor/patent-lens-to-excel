#!/usr/bin/env python3
"""spec.json（经 lensmath.py --write 后）→ Zemax .zmx，目录玻璃版 + 模型玻璃版，
自带多重结构与对焦组位置解，并在生成后反解析自校验。

    python3 make_zmx.py spec.matched.json -o WO2025253787A1_Ex01

spec 需要的 zmx 块:
  "zmx": {"fno":1.86, "fields_y":[...], "waves_um":[...],
          "focus":{"var_before":13,"var_after":26,"sum":10.26,
                   "key_before":"D13","key_after":"D26"},
          "configs":[{"name":"INF","d0":"INFINITY","D13":6.43}, ...]}
面上可选 "glass"（"HOYA FCD1"）、"extra"["有効径 φi"]（直径，脚本取一半写 DIAM）。
"""
import re
import json, argparse, math

ASPKEYS = ['A4','A6','A8','A10','A12','A14','A16']

# 默认光谱：可见光加权组，主波长 e 线(0.5461)。(波长 um, 权重)
WAVES_DEFAULT = [(0.4861, 12.0), (0.5461, 30.0), (0.6563, 3.0), (0.5876, 22.0), (0.4358, 3.0)]
PWAV_DEFAULT = 2
# 默认视场：最大像高的 6 等分点，由大到小（Real Image Height）
FIELD_FRACS = [1.0, 0.8, 0.6, 0.4, 0.2, 0.0]


ASPKEYS_FULL = ASPKEYS + ['A18', 'A20']
ASPKEYS_HIGH = ['A18', 'A20']          # EVENASPH 装不下的那两项

# Extended Asphere（OpticStudio 里叫 Extended Asphere，.zmx 里 TYPE 是 XASPHERE）的
# Extra Data 行。实证来源：用户机器上 5 个不同版本存出来的真文件
#   E:/Download/"135 1.4 ART DG.zmx" / "FE 16-35mm F2.8 GM II.zmx"
#   E:/Download/"GF 110mm F5.6 TS Macro.zmx" / "12mm F1.4 DC.zmx"
#   Documents/Zemax/Autosave/649/000.zmx
# 以及用户自己的 CODE V 宏 cv2zmx（github.com/anvcor/cv2zmx）的写法，两边完全一致：
#   XDAT 1 = 项数 N（真文件里恒为 10）
#   XDAT 2 = 归一化半径 Rn（恒为 1 —— 取 1 时系数就是专利印的 A2..A20，不用换算）
#   XDAT 3 = ρ² 项（= r² 项，专利没有这一项，必须留 0，否则近轴曲率被改掉）
#   XDAT 4..12 = r⁴ r⁶ r⁸ r¹⁰ r¹² r¹⁴ r¹⁶ r¹⁸ r²⁰ 的系数
# 面型式子与 Even Asphere 同源：z = cr²/(1+√(1-(1+k)c²r²)) + Σ αi·(r/Rn)^(2i)，
# 所以 CURV / CONI 原样照写，只是多项式搬到 Extra Data 里、能一直到 r²⁰。
XDAT_FMT = '  XDAT %d %.12E 0 0 1.000000000000E+00 0.000000000000E+00 0 ""'
XDAT_NTERMS = 10

# Extended Odd Asphere（.zmx 里 TYPE 是 XOSPHERE）—— 佳能这类「A3..A15 含奇数次」的
# 非球面，Even Asphere / Extended Asphere 都装不下（它们只有偶数次），Odd Asphere
# (ODDASPHE) 又只到 r^8（PARM 1..8，实证见 Documents/Zemax/.../老蛙视频_1.ZMX）。
# Extra Data 排布来自 OpticStudio 2024R2 用户手册 §2.3.1.2.28 的参数表
# （E:/ANSYS Inc/v242/Zemax OpticStudio/OpticStudio_UserManual_en.pdf，第 235 页）：
#   参数 13/14/15/16/17/18.. → XDAT 1/2/3/4/5/6..（偏移恒为 12）
#   XDAT 1 = 最大项号 N
#   XDAT 2 = 归一化半径 Rn（取 1.0，系数就是专利印的 A_n，不用换算）
#   XDAT 3 = ρ^1、XDAT 4 = ρ^2、…、XDAT k = ρ^(k-2)
# 同一张表在 Extended Asphere 那节给出的 13/14/15(ρ²)/16(ρ⁴) 与本机 5 个真 .zmx 文件
# 的 XDAT 1/2/3/4 完全对应，所以这个偏移是验证过的、不是猜的。
# 面型式子：z = cr²/(1+√(1-(1+k)c²r²)) + Σ αi·(r/Rn)^i —— CURV / CONI 原样照写。
XO_BASE = 2                 # XDAT (i + XO_BASE) = ρ^i 的系数

def _apow(A):
    """把非球面字典里的 A<n> 取成 {次数: 系数}，支持奇数次。"""
    out = {}
    for k, v in (A or {}).items():
        m = re.fullmatch(r'[Aa]\s*(\d+)', str(k))
        if m and v: out[int(m.group(1))] = float(v)
    return out


def _wls(cols, f, w):
    """加权最小二乘，MGS 正交化（列是 u^4..u^16，u∈[0,1]，正规方程会病态到没法用）。"""
    n, k = len(f), len(cols)
    sw = [math.sqrt(x) for x in w]
    Q = [[cols[j][i] * sw[i] for i in range(n)] for j in range(k)]
    y = [f[i] * sw[i] for i in range(n)]
    R = [[0.0] * k for _ in range(k)]
    for j in range(k):
        for i in range(j):
            dd = sum(Q[i][s] * Q[j][s] for s in range(n))
            R[i][j] = dd
            qi = Q[i]; qj = Q[j]
            for s in range(n): qj[s] -= dd * qi[s]
        nr = math.sqrt(sum(x * x for x in Q[j]))
        R[j][j] = nr
        if nr > 0:
            qj = Q[j]
            for s in range(n): qj[s] /= nr
    b = [sum(Q[j][s] * y[s] for s in range(n)) for j in range(k)]
    c = [0.0] * k
    for j in range(k - 1, -1, -1):
        s = b[j] - sum(R[j][m] * c[m] for m in range(j + 1, k))
        c[j] = s / R[j][j] if R[j][j] else 0.0
    return c


def refit_even(A, rmax, npts=1501, iters=50):
    """A18/A20 不为零时，把 9 项多项式在 [0, rmax] 上重新拟合进 Even Asphere 的 7 项
    (r^4..r^16)，曲率与 conic 原样不动（所以近轴一点不变）。

    为什么必须拟合而不是直接丢掉 A18：这类强非球面的各阶项在边缘是巨额相消 ——
    本篇面7 的 A18*r^18 在净口径处就有 2.1mm，丢掉它轴上 TA-RMS 会从 0.008mm 崩到 1.33mm。
    IRLS 逼近 minimax，残差比普通最小二乘小 3~4 倍。返回 (7 个系数, 最大矢高残差 mm)。"""
    co = [float(A.get(k) or 0.0) for k in ASPKEYS_FULL]
    if abs(co[7]) < 1e-30 and abs(co[8]) < 1e-30: return None
    us = [i / (npts - 1.0) for i in range(npts)]
    f = [sum(co[j] * rmax ** (4 + 2 * j) * u ** (4 + 2 * j) for j in range(9)) for u in us]
    cols = [[u ** (4 + 2 * j) for u in us] for j in range(7)]
    w = [1.0] * npts; best = None
    for _ in range(iters):
        c = _wls(cols, f, w)
        res = [sum(c[j] * cols[j][i] for j in range(7)) - f[i] for i in range(npts)]
        m = max(abs(x) for x in res)
        if best is None or m < best[0]: best = (m, c[:])
        w = [w[i] * (abs(res[i]) / m + 1e-3) ** 0.5 for i in range(npts)]
        sc = sum(w) / npts
        w = [x / sc for x in w]
    m, c = best
    return [c[j] / rmax ** (4 + 2 * j) for j in range(7)], m



def waves_of(zx):
    """zmx.waves 支持 [[um, weight], ...]；zmx.waves_um 是只给波长的老写法。"""
    if zx.get('waves'):
        return [(float(w[0]), float(w[1]) if len(w) > 1 else 1.0) for w in zx['waves']]
    if zx.get('waves_um'):
        return [(float(w), 1.0) for w in zx['waves_um']]
    return list(WAVES_DEFAULT)


def is_angle_field(zx):
    """zmx.field_type = "angle"：视场按物方角度（度）给，FTYP 0。
    半视场 >90° 的鱼眼必须这样 —— OpticStudio 的 Real Image Height 解不出 >90° 的主光线
    （US20220221688A1 RF5.2 Dual Fisheye：像高 8.55 ↔ 95.49°，FTYP 3 下视场1 整条追不出来），
    而 Angle 视场到 95°/100° 都能追（官方样例 Wide angle lens 200 degree field.zmx 就是 FTYP 0）。"""
    return str(zx.get('field_type', '')).lower() in ('angle', 'object_angle', 'yan')


def fields_cfg(zx, cfgs, base):
    """zmx.fields_cfg = {结构名: [视场值…]} → 逐结构视场表（缺的结构用 base）；没有就返回 None。"""
    fc = zx.get('fields_cfg')
    if not fc:
        return None
    miss = [k for k in fc if k not in [c['name'] for c in cfgs]]
    if miss:
        print('  ⚠ zmx.fields_cfg 的键对不上结构名：%s' % miss)
    return [[float(v) for v in fc.get(c['name'], base)] for c in cfgs]


def fields_of(zx, emb):
    """zmx.fields_y 直接给就用；否则按最大像高取 6 等分点。
    最大像高优先级：zmx.max_y > 各种数据里的 Y > 像面有効径/2。"""
    if zx.get('fields_y'):
        return [float(y) for y in zx['fields_y']]
    if is_angle_field(zx):                       # 视场按物方角度给（鱼眼 >90°）：最大半角 × 6 等分
        return [round(float(zx['max_angle']) * f, 4) for f in FIELD_FRACS]
    ymax = zx.get('max_y')
    if ymax is None:
        for k, v in emb.get('general', []):
            if str(k).startswith('Y') or '像高' in str(k):
                try: ymax = float(v); break
                except (TypeError, ValueError): pass
    if ymax is None:
        img = [s for s in emb['surfaces'] if s['i'] == 'IMG']
        phi = (img[0].get('extra') or {}).get('有効径 φi') if img else None
        if phi: ymax = phi / 2.0
    if ymax is None:
        raise SystemExit('无法确定最大像高：请在 spec 的 zmx 块里给 max_y 或 fields_y')
    return [round(ymax * f, 4) for f in FIELD_FRACS]

def _gauss_rings(n):
    """OpticStudio 的 Gaussian Quadrature 光瞳环：ρ² = (1+x)/2（x = n 点 Gauss-Legendre 根），环权重 = w/2。
    n=3 → ρ 0.33571/0.70711/0.94197、权重 5/18 / 8/18 / 5/18（与用户向导生成的文件逐位对上）。"""
    xs, ws = [], []
    for i in range(1, n + 1):
        x = math.cos(math.pi * (i - 0.25) / (n + 0.5))
        for _ in range(100):
            p0, p1 = 1.0, x
            for k in range(2, n + 1):
                p0, p1 = p1, ((2 * k - 1) * x * p1 - (k - 1) * p0) / k
            dp = n * (x * p1 - p0) / (x * x - 1.0)
            dx = p1 / dp; x -= dx
            if abs(dx) < 1e-15: break
        xs.append(x); ws.append(2.0 / ((1.0 - x * x) * dp * dp))
    pr = sorted((math.sqrt((1.0 + x) / 2.0), w / 2.0) for x, w in zip(xs, ws))
    return pr


def merit_contrast(fy, wv, ncfg, freq=80.0, rings=3, arms=6, fw=None, effl=None, pwav=PWAV_DEFAULT, effl_wt=1.0):
    """默认评价函数：优化向导「Contrast」s+t、freq lp/mm、GQ rings 环 arms 臂、无空气/玻璃约束，
    **逐结构铺一份**（CONF n + MECS/MECT）。用户 2026-09 在 OpticStudio 里手工生成的就是这套，
    本函数对 WO2024214585A1 那份 _OPT.zmx 的 480 个操作数逐个复现（权重误差 0）。

    行格式（照抄 OpticStudio 存出来的）：`MECS 0 <波长> <视场> <freq> <Px> <Py> 0 <权重> 0 0`
    权重 = 环权重 × 波长权重/最大波长权重 × 视场权重/最大视场权重 × 臂角权重；
    只有 Y 视场（旋转对称）时只追半个光瞳：离轴 arms/2 条臂（θ = 90° − (k+½)·360°/arms），
    每臂 π/(arms/2)；轴上视场只追 θ=0 一条臂、权重 π。

    effl：逐结构焦距目标（None 表示该结构不约束）。变焦镜头的焦距约束**集中放在最前面、DMFS 之前**：
      CONF 1 / BLNK 说明 / EFFL … / CONF 2 / EFFL … / CONF 3 / EFFL …（第 1 行必须是 CONF）
    `EFFL 0 <主波长> 0 0 0 0 <目标> <权重> 0 0`，防止优化对焦间隔时把焦距带跑（用户 2026-09-23 要求）。
    为什么放 DMFS 前：①打开 MFE 第一屏就能看到全部结构的焦距约束（放在各结构段首时 M/T 段在
    第 494/984 行，用户以为只有广角端）；②OpticStudio 重跑优化向导只替换 DMFS 之后的部分，前面的用户操作数保留。
    """
    fw = fw or [1.0] * len(fy)
    wmax = max(w for _, w in wv) or 1.0
    fmax = max(fw) or 1.0
    rr = _gauss_rings(rings)
    half = arms // 2
    out = []
    if effl and any(effl):
        # ★ 第 1 行必须是 CONF：多重结构下第 1 个操作数不是 CONF 时，OpticStudio 算评价函数会自己在最前面
        #   插一个 CONF 1，总数不变、最后一个操作数被挤掉（实测 T 结构少了 1 个 MECT）。
        first = True
        for c in range(1, ncfg + 1):
            tg = effl[c - 1]
            if tg:
                out.append('CONF %d 0 0 0 0 0 0 0 0 0' % c)
                if first:
                    out.append('BLNK Zoom focal length targets (one EFFL per zoom configuration).')
                    first = False
                out.append('EFFL 0 %d 0 0 0 0 %s %s 0 0' % (pwav, num(float(tg), '%.10G'), num(effl_wt, '%.10G')))
    out += ['CONF 1 0 0 0 0 0 0 0 0 0', 'DMFS 0 0 0 0 0 0 0 0 0 0',
           'BLNK contrast s+t S Wgt = 1.0000 T Wgt = 1.0000 Contrast at %s lp/MM GQ %d rings %d arms'
           % (num(freq, '%.10G'), rings, arms)]
    body = []
    for fi, y in enumerate(fy, 1):
        body.append('BLNK Operands for field %d.' % fi)
        if abs(y) < 1e-12:
            angs = [(0.0, math.pi)]
        else:
            angs = [(math.radians(90.0 - (k + 0.5) * 180.0 / half), math.pi / half) for k in range(half)]
        for wi, (_lam, ww) in enumerate(wv, 1):
            for th, aw in angs:
                for rho, rw in rr:
                    px = rho * math.cos(th); py = rho * math.sin(th)
                    if abs(px) < 1e-15: px = 0.0
                    if abs(py) < 1e-15: py = 0.0
                    wt = rw * (ww / wmax) * (fw[fi - 1] / fmax) * aw
                    for op in ('MECS', 'MECT'):
                        body.append('%s 0 %d %d %s %s %s 0 %s 0 0'
                                    % (op, wi, fi, num(freq, '%.10G'), repr(px) if px else '0',
                                       repr(py) if py else '0', repr(wt)))
    for c in range(1, ncfg + 1):
        out.append('CONF %d 0 0 0 0 0 0 0 0 0' % c)
        out.append('BLNK No air or glass constraints.')
        out.extend(body)
    return out


def focus_vars(spec, emb):
    """评价函数配套的变量：**对焦间隔**（逐结构 MCE THIC 置 Variable，∞ 结构也开 = 无穷远重新对焦）。

    - zmx.focus_vars 显式给了就用它；
    - 变焦：同一变焦位置的各结构之间会变的间隔 = 对焦间隔，按面序**去掉最后一个**
      （它是补偿段；用户 2026-09 的 _OPT.zmx：D18/D21 变量、D23 不动）；
    - 定焦：focus / focus2 的 key_before（key_after / key_last 由位置解 TOLE 跟随，守恒和自动保持）。
    """
    zx = spec['zmx']
    if zx.get('focus_vars') is not None:
        return list(zx['focus_vars'])
    cfgs = zx.get('configs') or []
    var = emb.get('variable') or {}
    if zx.get('zoom'):
        order = {s['D']: k for k, s in enumerate(emb['surfaces']) if isinstance(s['D'], str)}
        fk = set()
        groups = {}
        for c in cfgs: groups.setdefault(c.get('zoom'), []).append(c)
        for g in groups.values():
            for key in var:
                vs = [float(c[key]) for c in g if key in c]
                if len(vs) > 1 and max(vs) - min(vs) > 1e-9: fk.add(key)
        fk = sorted(fk, key=lambda k: order.get(k, 1e9))
        return fk[:-1] if len(fk) > 1 else fk
    out = []
    for f in (zx.get('focus'), zx.get('focus2')):
        if f and f.get('key_before'): out.append(f['key_before'])
    return out


def num(x, fmt='%.10G'):
    return fmt % x

def ascii_(t):
    """.zmx 必须是 latin-1；中文名转成可读的 ASCII 记号。"""
    t = str(t)
    try:
        t.encode('latin-1'); return t
    except UnicodeEncodeError:
        import re
        m = re.search(r'(\d+)', t)
        base = re.sub(r'[^\x20-\x7E]+', '', t).strip(' _-')
        if not base or base.isdigit():          # 纯中文名（如「实施例1」）只剩下数字
            return 'Example' + (('_' + base) if base else '')
        return base

def wfno_cfg(spec, emb, with_stop=False, with_beta=False):
    """逐结构近轴工作 F 数（Paraxial Working F/#）—— **每次都现算**（纯近轴，毫秒级）。

    `zmx.aperture`（vignet.py --write 写的）只拿来对照：vignet 之后又改了 fno / fno_patent /
    wfno_override / 结构名，缓存就过期了（回归审查：改 wfno_override 不生效、改结构名直接崩）。
    过期时以现算为准并告警 —— 渐晕是按旧孔径解的，要重跑 vignet。
    """
    from vignet import aperture_cfg
    zx = spec['zmx']
    cfgs = zx.get('configs') or [{'name': 'INF', 'd0': 'INFINITY'}]
    _s0, rows = aperture_cfg(spec, emb, emb['states'][0], cfgs)
    wf = [r['wfno'] for r in rows]
    ss = [r['stop_semi'] for r in rows]
    bt = [r['beta'] for r in rows]
    for t in getattr(aperture_cfg, 'notes', []):
        print('  ⚠ 孔径: ' + t)
    ap = zx.get('aperture') or {}
    if ap.get('wfno_cfg'):
        same = (ap.get('configs') == [c.get('name') for c in cfgs]
                and len(ap['wfno_cfg']) == len(wf)
                and all(abs(float(a) - b) < 1e-3 for a, b in zip(ap['wfno_cfg'], wf)))
        if not same:
            print('  ★ zmx.aperture 已过期（vignet 之后改过 fno / fno_patent / wfno_override / 结构）：'
                  '以现算 %s 为准；渐晕是按旧孔径解的，要重跑 vignet.py' % [round(v, 4) for v in wf])
    if with_beta:
        return wf, ss, bt
    return (wf, ss) if with_stop else wf


def build(spec, emb, catalog, asph_mode='auto', merit=True, freq=80.0):
    zx = spec['zmx']
    # 变焦镜头：结构里写全了所有可变间隔（lensmath.zoom_configs），不用位置解 ——
    # 每个变焦位置的守恒和不同，TOLE 的长度进不了 MCE。每个可变间隔都逐结构进 MCE 的 THIC。
    zoom = bool(zx.get('zoom'))
    fc = {} if zoom else (zx.get('focus') or {})
    kb, ka = fc.get('key_before'), fc.get('key_after')
    # 链式三段浮动（d_a + d_b + d_c = const，如本篇适马 105 微距：G2 与絞り各自移动）：
    # 第三段 key_last 的 DISZ = 守恒和 − 前两段，位置解随后会覆盖它。
    kl = fc.get('key_last')
    # 第二个对焦群（双浮动对焦，如索尼 135GM）。没有就是 None，行为与以前完全一致。
    fc2 = None if zoom else zx.get('focus2')
    kb2, ka2 = (fc2['key_before'], fc2.get('key_after')) if fc2 else (None, None)
    cfgs = zx['configs']
    surfs = [s for s in emb['surfaces'] if s['i'] != 'IMG']
    img   = [s for s in emb['surfaces'] if s['i'] == 'IMG']
    asp = {a['surface'].replace('面',''): a for a in emb.get('aspheric', [])}
    gmap = (spec.get('zmx') or {}).get('gcat') or {}
    # GCAT 必须写用户机器上真实存在的目录文件名（不带扩展名）。
    # OpticStudio 自带的 HOYA.AGF / OHARA.AGF 常常是好几年前的版本，
    # 新牌号（如 OHARA S-LAL18N）不在里面，写 GCAT HOYA OHARA 会打不开。
    _gv = [s['glass'] for s in surfs if s.get('glass')] + \
          [s['glass_offset']['base'] for s in surfs if s.get('glass_offset')]
    vendors = sorted({gmap.get(g.split()[0], g.split()[0]) for g in _gv})
    # ---- 非球面面型选择：EVENASPH 只到 r^16，有 A18/A20 就换 XASPHERE ----
    # 与用户的 CODE V 宏 cv2zmx 同一套策略（^opt_asp）：逐面判断，
    # H(r^18)/J(r^20) 一旦非零就换 Extended Asphere，其余仍用 Even Asphere。
    atype, refit, rf_notes, xa_notes, xo_notes = {}, {}, [], [], []
    xo_nterms = XDAT_NTERMS
    for s in surfs:
        key = str(s['i']); A = asp.get(key)
        if not A: continue
        pw = _apow(A)
        odd = sorted(n for n in pw if n % 2)
        if odd:
            # 奇数次项存在 → 只有 Extended Odd Asphere 能精确表达。
            # 不许丢项、不许拟合：这类面的各阶在边缘是巨额相消。
            atype[key] = 'XOSPHERE'
            xo_notes.append(key)
            xo_nterms = max(xo_nterms, max(pw))
            continue
        hi = any(abs(float(A.get(k) or 0.0)) > 0.0 for k in ASPKEYS_HIGH)
        if asph_mode == 'extended' or (asph_mode == 'auto' and hi):
            atype[key] = 'XASPHERE'
            if hi: xa_notes.append(key)
            continue
        atype[key] = 'EVENASPH'
        if not hi: continue
        # 被 --asph-type even 强按回 Even Asphere：只能在净口径上把 9 项重拟合进 7 项。
        phi = (s.get('extra') or {}).get('有効径 φi')
        if not phi:
            R_ = s['R']
            phi = abs(float(R_)) if R_ not in (None, 0, 'inf') else 0.0
        rr = refit_even(A, float(phi) / 2.0)
        if rr:
            refit[key] = rr[0]
            rf_notes.append((key, float(phi) / 2.0, rr[1],
                             float(A.get('A18') or 0.0), float(A.get('A20') or 0.0)))
    L=[]; a=L.append
    a('VERS 190513 693 1 L000001')
    a('MODE SEQ')
    a('NAME %s' % (ascii_(zx['name']) if zx.get('name')
                   else '%s %s' % (ascii_(spec['patent']), ascii_(emb['name']))))
    for n in spec.get('zmx_notes', []): a('NOTE 0 %s' % ascii_(n))
    if xa_notes:
        a('NOTE 0 %s' % ascii_(
            'Aspheres: surfaces %s carry r^18 / r^20 terms, which Even Asphere cannot hold'
            % ', '.join(xa_notes)))
        a('NOTE 0 %s' % ascii_(
            '  (it stops at r^16), so they are Extended Asphere (TYPE XASPHERE). Coefficients'))
        a('NOTE 0 %s' % ascii_(
            '  are in the Extra Data Editor: 1=number of terms (10), 2=norm radius (1.0),'))
        a('NOTE 0 %s' % ascii_(
            '  3=r^2 (kept 0), 4..12 = r^4..r^20 = the printed A4..A20. Curvature and conic'))
        a('NOTE 0 %s' % ascii_(
            '  are unchanged, so the surfaces are bit-exact against the patent table.'))
        a('NOTE 0 %s' % ascii_(
            '  Norm radius must stay 1.0 - change it and every coefficient must be rescaled'))
        a('NOTE 0 %s' % ascii_(
            '  by Rn^(2i); term 3 must stay 0 or the paraxial curvature changes.'))
    if xo_notes:
        a('NOTE 0 %s' % ascii_(
            'Aspheres: surfaces %s use ODD as well as even powers of r (the patent prints'
            % ', '.join(xo_notes)))
        a('NOTE 0 %s' % ascii_(
            '  A3..A%d). Even Asphere and Extended Asphere hold even powers only, and Odd' % xo_nterms))
        a('NOTE 0 %s' % ascii_(
            '  Asphere stops at r^8, so these are Extended Odd Asphere (TYPE XOSPHERE).'))
        a('NOTE 0 %s' % ascii_(
            '  Extra Data: 1=max term number (%d), 2=norm radius (1.0), 3..%d = r^1..r^%d.'
            % (xo_nterms, xo_nterms + 2, xo_nterms)))
        a('NOTE 0 %s' % ascii_(
            '  Terms r^1 and r^2 are kept 0 (the patent series starts at A3); curvature and'))
        a('NOTE 0 %s' % ascii_(
            '  conic are unchanged, so the surfaces are bit-exact against the patent table.'))
        a('NOTE 0 %s' % ascii_(
            '  Norm radius must stay 1.0 - change it and every coefficient needs Rn^i.'))
    if rf_notes:
        a('NOTE 0 %s' % ascii_(
            'Aspheres: --asph-type even was forced, so surfaces with r^18/r^20 were REFITTED'))
        a('NOTE 0 %s' % ascii_(
            '  into r^4..r^16 (curvature + conic unchanged, paraxial identical):'))
        for key, rr, err, a18, a20 in rf_notes:
            a('NOTE 0 %s' % ascii_(
                '  surf %s: fitted over r<=%.3f mm, max sag error %.1f nm; exact A18 = %.5E'
                % (key, rr, err * 1e6, a18)))
    a('NOTE 1 ""')
    a('PFIL 0 0 0'); a('LANG 0'); a('UNIT MM X W X CM MR CPMM')
    # 孔径类型固定写 **Paraxial Working F/#**（FNUM <值> 1；0 = Image Space F/#）。
    # ★ 不写 Image Space F/#：它按「该结构自己的 ∞ 共轭 EFL / 入瞳直径」定义，内对焦/浮动对焦镜头
    #   近距 EFL 大幅缩短时光阑会被跟着缩小 —— 用户实测 JP2021-047297A 1.00x 结构 ENPD 15.73、
    #   状态栏 WFNO 9.149（真值 ≈4.13）。也不再写 ENPD：用户要求统一走近轴工作 F 数。
    # 逐结构的值由「物理光阑固定」算出（见 vignet.aperture_cfg），∞ 结构恰好等于专利 F 数。
    wf = wfno_cfg(spec, emb)
    a('FNUM %s 1' % num(round(wf[0], 6)))
    a('ENVD 20 1 0'); a('GFAC 0 0')
    if catalog and vendors: a('GCAT ' + ' '.join(vendors) + ' ')
    # RAIM 第 2 位 = Ray Aiming：0=Off / 1=Paraxial / 2=Real（官方样例实证，见 SKILL）。
    # 大孔径 / 大视场 / 光阑前有强负透镜时必须开 Real，否则光瞳定位是错的。
    ra = {'off': 0, 'none': 0, 'paraxial': 1, 'real': 2}.get(
        str(zx.get('ray_aiming', 'real')).lower(), 2)
    a('RAIM 0 %d 1 1 0 0 0 0 0 1' % ra); a('PUSH 0 0 0 0 0 0'); a('SDMA 0 1 0')
    a('OMMA 1 1')
    fy = fields_of(zx, emb); wv = waves_of(zx)
    a('FTYP %d 0 %d %d 0 0 0' % (0 if is_angle_field(zx) else 3, len(fy), len(wv)))
    a('ROPD 2'); a('HYPR 0'); a('PICB 1')
    n12 = max(12, len(fy))
    a('XFLN ' + ' '.join(['0']*n12))
    a('YFLN ' + ' '.join([num(y) for y in fy] + ['0']*(n12-len(fy))))
    a('FWGN ' + ' '.join(['1']*len(fy) + ['0']*(n12-len(fy))))
    vig = zx.get('vignetting')
    if vig:
        cols = list(zip(*vig))                       # [(vdx..),(vdy..),(vcx..),(vcy..)]
        for t, col in zip(('VDXN','VDYN','VCXN','VCYN'), cols):
            a(t + ' ' + ' '.join([num(x) for x in col] + ['0']*(n12-len(col))))
        a('VANN ' + ' '.join(['0']*n12))
    else:
        for t in ('VDXN','VDYN','VCXN','VCYN','VANN'): a(t + ' ' + ' '.join(['0']*n12))
    for i, (lam, wt) in enumerate(wv, 1):
        a('WAVM %d %s %s' % (i, num(lam), num(wt)))
    for i in range(len(wv)+1, 25): a('WAVM %d 0.55 0' % i)
    a('PWAV %d' % zx.get('primary_wave', PWAV_DEFAULT))
    a('POLS 1 0 1 0 0 1 0')
    a('GLRS 1 0')
    a('GSTD 0 100.000 100.000 100.000 100.000 100.000 100.000 0 1 1 0 0 1 1 1 1 1 1')
    a('NSCD 100 500 0 9.9999999999999995e-07 10 9.9999999999999995e-07 0 0 0 0 0 1 1000000 0 2')
    # COFN 必须给满 4 个数据文件名：涂层 / 散射 / ABG / 面形。
    # 只写一个的话，OpticStudio 会拿空文件名去开 PROFILES 目录，
    # 弹「Can't open File ...\\Zemax\\PROFILES\\!」。
    a('COFN QF "COATING.DAT" "SCATTER_PROFILE.DAT" "ABG_DATA.DAT" "PROFILE.GRD"')
    a('COFN COATING.DAT SCATTER_PROFILE.DAT ABG_DATA.DAT PROFILE.GRD')
    a('SURF 0'); a('  TYPE STANDARD'); a('  CURV 0.0 0 0 0 0 ""')
    a('  HIDE 0 0 0 0 0 0 0 0 0 0'); a('  MIRR 2 0')
    a('  DISZ INFINITY'); a('  DIAM 0 0 0 0 1 ""')
    a('  POPS 0 0 0 0 0 0 0 0 1 1 1 1 0 0 0 0')
    # 位置解长度 = var_before 面 到 var_after 面 之间全部厚度之和（含两端）
    def _poslen(f_):
        return f_['sum'] + sum(float(s['D']) for s in surfs
                               if isinstance(s['D'], (int, float))
                               and f_['var_before'] < _idx(s) < f_['var_after'])
    poslen = _poslen(fc) if fc.get('var_after') and 'sum' in fc else None
    poslen2 = _poslen(fc2) if fc2 and fc2.get('var_after') else None   # 链式三段：focus2 无 var_after
    varsurf = {}                                   # 可变间隔名 → 所在面号（变焦时逐个进 MCE）
    for k, s in enumerate(surfs, 1):
        i = _idx(s); key = str(s['i'])
        A = asp.get(key)
        a('SURF %d' % k)
        ty = atype.get(key, 'STANDARD') if A else 'STANDARD'
        doe = s.get('doe')
        if doe:
            # 衍射面 → Binary 2（实证：Samples/Sequential/Diffractive components/Achromatic singlet.zmx）
            #   PARM 0 = 衍射级次 M；PARM 1..8 = 偶次非球面（这里 0）
            #   XDAT 1 = 最大项号 N；XDAT 2 = 归一化半径 Rn；XDAT 2+i = ρ^(2i) 的相位系数（rad）
            # 专利 ψ = 2π/λ0·Σ C_2i h^2i ⇒ A_i = 2π/λ0[mm]·C_2i·Rn^(2i)，取 Rn=1 不用换算。
            # Zemax 按每条光线自己的 λ 算偏折（λ/2π·∇Φ），色散自动正确。
            if A: raise SystemExit('面%s 同时有非球面与 DOE，Binary 2 的 PARM 需另行处理' % key)
            ty = 'BINARY_2'
        a('  TYPE %s' % ty)
        c = 0.0 if s['R'] in (None, 0, 'inf') else 1.0/float(s['R'])
        a('  CURV %s 0 0 0 0 ""' % num(c, '%.12G'))
        a('  HIDE 0 0 0 0 0 0 0 0 0 0'); a('  MIRR 2 0')
        if s.get('stop') or s['i'] == 'STO': a('  STOP')
        if ty == 'BINARY_2':
            a('  PARM 0 1')
            for j in range(1, 9): a('  PARM %d 0' % j)
            lam = float(doe.get('wl_nm', 587.56)) * 1e-6
            CC = doe['C']
            a('  XDAT 1 %.12E 0 0 0.000000000000E+00 0.000000000000E+00 0 ""' % float(len(CC)))
            a('  XDAT 2 %.12E 0 0 0.000000000000E+00 0.000000000000E+00 0 ""' % 1.0)
            for j, cval in enumerate(CC, 1):
                a('  XDAT %d %.12E 0 0 0.000000000000E+00 0.000000000000E+00 0 ""'
                  % (j + 2, 2 * math.pi / lam * float(cval)))
        elif ty == 'EVENASPH':
            a('  PARM 1 0')                      # PARM1 = r^2 项，专利没有，恒 0
            rc = refit.get(key)
            for j, kk in enumerate(ASPKEYS, 2):
                v = rc[j - 2] if rc else (A.get(kk) or 0.0)
                a('  PARM %d %s' % (j, num(v, '%.12G')))
        elif ty == 'XASPHERE':
            a(XDAT_FMT % (1, float(XDAT_NTERMS)))   # 项数
            a(XDAT_FMT % (2, 1.0))                  # 归一化半径 Rn = 1
            a(XDAT_FMT % (3, 0.0))                  # r^2 项，必须 0
            for j, kk in enumerate(ASPKEYS_FULL, 4):
                a(XDAT_FMT % (j, float(A.get(kk) or 0.0)))
        elif ty == 'XOSPHERE':
            pw = _apow(A)
            a(XDAT_FMT % (1, float(xo_nterms)))     # 最大项号 N
            a(XDAT_FMT % (2, 1.0))                  # 归一化半径 Rn = 1
            for n in range(1, xo_nterms + 1):       # XDAT n+2 = ρ^n 的系数
                a(XDAT_FMT % (n + XO_BASE, pw.get(n, 0.0)))
        D = s['D']
        if isinstance(D, str): varsurf[D] = i
        if isinstance(D, str) and zoom:
            D = cfgs[0][D] if D in cfgs[0] else float(emb['variable'][D][emb['states'][0]])
        elif isinstance(D, str):
            if D == kb:
                D = [c0 for c0 in cfgs if kb in c0][0][kb]
            elif D == ka:
                D = fc['sum'] - cfgs[0][kb]
            elif kb2 and D == kb2:
                D = [c0 for c0 in cfgs if kb2 in c0][0][kb2]
            elif ka2 and D == ka2:
                D = fc2['sum'] - cfgs[0][kb2]
            elif kl and D == kl:
                D = fc['sum'] - cfgs[0][kb] - cfgs[0][kb2]
            else:
                # 对焦组以外的第三个可变量（常见的是末面 BF=D<n>）——按基准状态取值，
                # 绝不能当成 key_after 走 sum-kb（那会把 BF 写成对焦间隔的值）
                st0 = emb.get('states', [None])[0]
                D = float(emb['variable'][D][st0])
        a('  DISZ %s' % num(float(D), '%.6G'))
        if poslen is not None and i == fc['var_after']:
            a('  TOLE %d %s' % (fc['var_before'], num(poslen, '%.6G')))
        elif poslen2 is not None and i == fc2['var_after']:
            a('  TOLE %d %s' % (fc2['var_before'], num(poslen2, '%.6G')))
        a('  CONI %s' % num((A or {}).get('k', (A or {}).get('K', 0.0))))
        if s.get('nd'):
            # Offset 玻璃解（solve code 4）：基准目录玻璃 + Nd/Vd 偏移。
            # 用来给「无等效牌号」的面保住真实色散曲线，又精确还原专利印刷的 nd/vd
            # —— 比模型玻璃好得多，用户点名要这样做。
            # 实测 OpticStudio 2024R2 存出来的格式（字段 4/5/6 = 基准玻璃的 nd/vd/dPgF，
            # 两个偏移量在**最后两个字段**，不是 4/5）：
            #   GLAS Q-PSKH1S 4 0 1.59255001 67.85687061 0.0136 0 0 0 0.00039 0.0431
            go = s.get('glass_offset')
            if catalog and go:
                a('  GLAS %s 4 0 %s %s %s 0 0 0 %s %s'
                  % (go['base'].split()[-1], num(go['base_nd'], '%.17G'),
                     num(go['base_vd'], '%.17G'), num(go.get('base_dpgf') or 0.0, '%.6G'),
                     num(go['d_nd'], '%.6G'), num(go['d_vd'], '%.6G')))
            elif catalog and s.get('glass'):
                a('  GLAS %s 0 0 %s %s 0 0 0 0 0 0'
                  % (s['glass'].split()[1], num(s['nd'],'%.6G'), num(s['vd'],'%.6G')))
            else:
                a('  GLAS ___BLANK 1 0 %s %s %s 0 0 0 0 0'
                  % (num(s['nd'],'%.6G'), num(s['vd'],'%.6G'),
                     num(s.get('dpgf') or 0.0, '%.6G')))
        phi = (s.get('extra') or {}).get('有効径 φi')
        # DIAM 的第 2 个数：0=自动（Zemax 按实际光线算），1=固定（当硬光阑用）。
        # 专利的「有効径」是刚好够用的近轴口径，零余量；固定住会把含球差的实际
        # 边缘光线切掉（实测轴上光线在专利 φ/2 的 93~96% 处，实光线会超）。
        # 所以默认自动，数值仍写专利值，翻成固定时就是专利口径。
        fixset = zx.get('fix_semi_surfaces')
        if fixset is not None:
            try: fix = 1 if int(s['i']) in fixset else 0
            except (TypeError, ValueError): fix = 0
        else:
            fix = 1 if zx.get('fix_semi') else 0
        mg = float(zx.get('semi_margin', 0.0))
        if s.get('stop') or s['i'] == 'STO' or not phi:
            a('  DIAM 0 0 0 0 1 ""'); a('  MEMA 0 0 0 0 1 ""')
        else:
            sd = num(phi / 2.0 * (1 + mg), '%.6G')
            a('  DIAM %s %d 0 0 1 ""' % (sd, fix)); a('  MEMA %s 0 0 0 1 ""' % sd)
        a('  POPS 0 0 0 0 0 0 0 0 1 1 1 1 0 0 0 0')
    n = len(surfs)+1
    a('SURF %d' % n); a('  TYPE STANDARD'); a('  CURV 0.0 0 0 0 0 ""')
    a('  HIDE 0 0 0 0 0 0 0 0 0 0'); a('  MIRR 2 0'); a('  DISZ 0'); a('  CONI 0')
    ip = (img[0].get('extra') or {}).get('有効径 φi') if img else None
    sd = num(ip/2.0,'%.6G') if ip else '0'
    a('  DIAM %s 0 0 0 1 ""' % sd); a('  MEMA %s 0 0 0 1 ""' % sd)
    a('  POPS 0 0 0 0 0 0 0 0 1 1 1 1 0 0 0 0')
    # ===== 配套评价函数（打开就能直接优化）=====
    # 优化向导 Contrast s+t 80 lp/mm GQ 3 环 6 臂，逐结构一份；变量 = 对焦间隔（下面 MCE 的 THIC 状态位 1）。
    if merit:
        # 变焦：每个 ∞ 结构加 EFFL 目标 = 专利该变焦位置的标称焦距（近距结构不约束）
        effl = None
        zpos = {p['name']: p for p in ((zx.get('zoom') or {}).get('positions') or [])}
        if zpos:
            effl = []
            for c0 in cfgs:
                inf = str(c0.get('d0', '')).upper().startswith('INF')
                zp = zpos.get(c0.get('zoom'))
                effl.append((zp.get('f') or c0.get('efl')) if (inf and zp) else None)
        L.extend(merit_contrast(fy, wv, len(cfgs), freq, effl=effl,
                                pwav=zx.get('primary_wave', PWAV_DEFAULT),
                                effl_wt=float(zx.get('effl_weight', 1.0))))
    fv = set(focus_vars(spec, emb)) if merit else set()
    a('TOL TOFF   0   0              0              0   0 0 0 0')
    a('MNUM %d 1' % len(cfgs))
    for i, c0 in enumerate(cfgs, 1):
        a('LTTL   0   %d "%s" 0 0 0 1 1 1 0 0' % (i, ascii_(c0['name'])))
    for i, c0 in enumerate(cfgs, 1):
        d0 = c0['d0']
        v = '1.00000000E+10' if str(d0).upper().startswith('INF') else num(float(d0), '%.10G')
        a('THIC   0   %d %s 0 0 0 1 1 1 0 0' % (i, v))
    if zoom:
        for key, sn in sorted(varsurf.items(), key=lambda t: t[1]):
            vals = [float(c0[key]) if key in c0 else float(emb['variable'][key][emb['states'][0]])
                    for c0 in cfgs]
            if len(cfgs) > 1 and all(abs(v - vals[0]) < 1e-9 for v in vals):
                continue                    # 所有结构都一样的就不必进 MCE
            for i, v in enumerate(vals, 1):
                # THIC 行第 4 个字段 = 状态：0 固定 / 1 变量（照 OpticStudio 存出来的 _OPT.zmx）
                a('THIC  %2d   %d %s %d 0 0 1 1 1 0 0' % (sn, i, num(v, '%.10G'), 1 if key in fv else 0))
    elif kb:
        for i, c0 in enumerate(cfgs, 1):
            a('THIC  %2d   %d %s %d 0 0 1 1 1 0 0' % (fc['var_before'], i, num(c0[kb], '%.10G'), 1 if kb in fv else 0))
    if fc2:
        for i, c0 in enumerate(cfgs, 1):
            a('THIC  %2d   %d %s %d 0 0 1 1 1 0 0' % (fc2['var_before'], i, num(c0[kb2], '%.10G'), 1 if kb2 in fv else 0))
    # 链式三段浮动（如适马 85 Art / JP2018-5099 実施例4）：三个可变间隔互不成对，
    # 位置解只能钉住其中一个，剩下的第三个必须自己进 MCE。
    # zmx.mce_extra = [{"var": 27, "key": "d27"}, ...]
    for ex in zx.get('mce_extra', []):
        for i, c0 in enumerate(cfgs, 1):
            a('THIC  %2d   %d %s %d 0 0 1 1 1 0 0'
              % (ex['var'], i, num(float(c0[ex['key']]), '%.10G'), 1 if ex['key'] in fv else 0))
    # ===== 逐结构渐晕：APER + FVCY/FVCX/FVDY/FVDX =====
    # Zemax 的 VDX/VDY/VCX/VCY 是**全局**量，只写在文件头里就等于所有结构共用一套。
    # 对焦镜头每个结构的渐晕差得很远（适马70微距 视场1 的 VDY 从 ∞ 的 +0.12 走到
    # 1:1 的 −0.34），不进多重结构的话近距结构的光瞳会整片被口径切掉 —— 用户实测打回过。
    # 操作数第 1 个参数是**视场号**（不是面号），第 2 个是结构号。
    # APER 在 Paraxial Working F/# 下就是该结构的近轴工作 F 数 —— 逐结构不同，与渐晕无关，总是要写。
    if len(cfgs) > 1:
        for i, v in enumerate(wf, 1):
            a('APER   0   %d %s 0 0 0 1 1 1 0 0' % (i, num(round(v, 6), '%.6G')))
    # ===== 逐结构视场值：YFIE（zmx.fields_cfg = {结构名: [视场值…]}，Zemax 顺序由大到小）=====
    # 鱼眼（>90° 角度视场）的近距结构：物面是平面，95° 的主光线在平物面上没有交点 ——
    # OpticStudio 会把它折成 −84.6°、像高跑到另一侧（US20220221688A1 实测）。近距结构只能取 <90°。
    fcfg = fields_cfg(zx, cfgs, fields_of(zx, emb))
    if fcfg:
        for f in range(len(fcfg[0])):
            vals = [fcfg[ci][f] for ci in range(len(cfgs))]
            if len(set(vals)) == 1:
                continue
            for ci, v in enumerate(vals, 1):
                a('YFIE  %2d   %d %s 0 0 0 1 1 1 0 0' % (f+1, ci, num(float(v), '%.6G')))
    vcfg = zx.get('vignetting_cfg')
    if vcfg and len(vcfg) == len(cfgs):
        for op, ix in (('FVCY', 3), ('FVCX', 2), ('FVDY', 1), ('FVDX', 0)):
            for f in range(len(vcfg[0])):
                vals = [vcfg[ci][f][ix] for ci in range(len(cfgs))]
                if all(abs(v) < 1e-9 for v in vals):
                    continue          # 全 0 的视场行省掉：缺行时该视场直接用文件头的全局值
                for ci, v in enumerate(vals, 1):
                    a('%s  %2d   %d %s 0 0 0 1 1 1 0 0' % (op, f+1, ci, num(float(v), '%.6G')))
    # ★ 文件末尾**不写** `CONF 1`：.zmx 里的 CONF 行是**评价函数操作数**（OpticStudio 存出来的文件
    #   只在 MFE 段里有 CONF），写在 MCE 后面会被读成评价函数末尾多出来的一个 CONF 操作数
    #   （用户 _OPT.zmx 2931 个操作数，旧写法读出 2932）。
    return ('\r\n'.join(L) + '\r\n').encode('latin-1', 'replace')

def _idx(s):
    try: return int(s['i'])
    except (TypeError, ValueError): return {'STO': None}.get(s['i'], None) or _idx.last
_idx.last = 0
def idx_of(surfs):
    n = 0
    for s in surfs:
        try: n = int(s['i'])
        except (TypeError, ValueError): n += 1
        s['_i'] = n

def _read_zmx(path):
    raw = open(path, 'rb').read()
    if raw[:2] in (b'\xff\xfe', b'\xfe\xff'): return raw.decode('utf-16')
    return raw.decode('latin-1')

def verify(path, patent_f, expect_stop=None):
    txt = _read_zmx(path).splitlines()
    S=[]; cur=None; mce={}; names=[]; aper={}; fnum=None
    for ln in txt:
        t = ln.strip()
        if t.startswith('SURF '): cur={'c':0.,'d':0.,'nd':None,'stop':False,'tole':None}; S.append(cur)
        elif cur is not None:
            if t.startswith('TYPE '): cur['type']=t.split()[1]
            elif t.startswith('XDAT '):
                cur['xdat']=cur.get('xdat',0)+1
                q=t.split()
                cur.setdefault('xv',{})[int(q[1])]=float(q[2])
                if q[1]=='1': cur['nterm']=int(float(q[2]))   # XDAT 1 = 项数/最大项号
            elif t.startswith('CURV '): cur['c']=float(t.split()[1])
            elif t.startswith('DISZ '):
                v=t.split()[1]; cur['d']=float('inf') if v.upper()=='INFINITY' else float(v)
            elif t.startswith('GLAS '):
                q=t.split()
                # Offset 玻璃解（code 4）：字段4/5 是**基准玻璃**的 nd/vd，
                # 真正生效的折射率 = 基准 nd + 最后第 2 个字段的 Nd offset。
                # 不加这一步，自校验会拿基准玻璃去算 EFL，把 Offset 的意义完全抹掉
                # （实测 JP2025052870A：34.4045 会被算成 34.4369）。
                cur['nd']=float(q[4]) + (float(q[10]) if len(q)>10 and q[2]=='4' else 0.0)
            elif t.startswith('STOP'): cur['stop']=True
            elif t.startswith('TOLE '): q=t.split(); cur['tole']=(int(q[1]),float(q[2]))
        if t.startswith('LTTL'):
            import re; names.append(re.search(r'"([^"]*)"', t).group(1))
        if t.startswith('THIC'):
            q=t.split(); mce.setdefault(int(q[1]),{})[int(q[2])]=float(q[3])
        if t.startswith('APER'):
            q=t.split(); aper[int(q[2])]=float(q[3])
        if t.startswith('FNUM '):
            q=t.split(); fnum=(float(q[1]), int(q[2]) if len(q)>2 else 0)
    body=S[1:-1]
    n=1.; y=1.; u=0.
    for k,s in enumerate(body):
        n2=s['nd'] or 1.
        c2=0.0
        if s.get('type')=='BINARY_2':
            xv=s.get('xv',{}); rn=xv.get(2,1.0)
            c2=xv.get(3,0.0)/(2*math.pi/(587.56e-6))/rn**2   # 注意：假定基准波长 587.56nm
        u=(n*u-y*s['c']*(n2-n)+2*c2*y)/n2
        if k<len(body)-1: y+=u*s['d']
        n=n2
    print('  自校验: 面数=%d 玻璃面=%d 光阑=面%s 位置解=%s' %
          (len(body), sum(1 for s in body if s['nd']),
           [i+1 for i,s in enumerate(body) if s['stop']],
           [(i+1,s['tole']) for i,s in enumerate(body) if s['tole']]))
    xo = [(i+1, s.get('xdat', 0)) for i, s in enumerate(body) if s.get('type') == 'XOSPHERE']
    if xo:
        nt = max(body[i-1].get('nterm', 0) for i, _ in xo)
        bad = [i for i, n in xo if n != nt + 2 or body[i-1].get('nterm') != nt]
        print('  自校验: Extended Odd Asphere(XOSPHERE) 面%s  XDAT 行数 %s%s'
              % ([i for i, _ in xo], sorted({n for _, n in xo}),
                 '  ★ XDAT 条数不对：面%s' % bad if bad else
                 '  (= 2 + %d 项 r^1..r^%d，正确)' % (nt, nt)))
    xa = [(i+1, s.get('xdat', 0)) for i, s in enumerate(body) if s.get('type') == 'XASPHERE']
    ev = [i+1 for i, s in enumerate(body) if s.get('type') == 'EVENASPH']
    if ev: print('  自校验: Even Asphere 面%s' % ev)
    if xa:
        bad = [i for i, n in xa if n != XDAT_NTERMS + 2]
        print('  自校验: Extended Asphere(XASPHERE) 面%s  XDAT 行数 %s%s'
              % ([i for i, _ in xa], sorted({n for _, n in xa}),
                 '  ★ 面%s 的 XDAT 条数不对' % bad if bad else '  (= 2 + 10 项，正确)'))
    print('  自校验: 反解析 EFL=%.4f (专利 %.2f)  结构=%s' % (-1/u, patent_f, names))
    # ---- 孔径：逐结构由 APER（近轴工作 F 数）反推光阑近轴半径，与孔径模型逐结构比对 ----
    if fnum:
        print('  自校验: 孔径类型 %s  值 %.4f' % ({0: 'Image Space F/#', 1: 'Paraxial Working F/#'}
                                                .get(fnum[1], '?%d' % fnum[1]), fnum[0]))
    if fnum and fnum[1] == 1:
        rows = cfg_paraxial(path)
        semis = []
        for ci, r in enumerate(rows, 1):
            F = aper.get(ci, fnum[0]) if len(rows) > 1 else fnum[0]
            semis.append(abs(r['ystop']) / (abs(r['u_img']) * 2 * F))   # 光阑半径 = |y_stop|·(1/2F)/|u'|
        spread = (max(semis) - min(semis)) / max(semis)
        if expect_stop and len(expect_stop) == len(semis):
            dev = max(abs(x - y) / y for x, y in zip(semis, expect_stop))
            if dev < 5e-4:
                verdict = '(与孔径模型逐结构一致 ✓，最大偏差 %.3f%%)' % (100*dev)
            else:
                verdict = '★ 与孔径模型不一致（最大偏差 %.2f%%）—— APER 写错了或结构间隔没对上' % (100*dev)
        elif spread < 1e-4:
            verdict = '(各结构一致，物理光阑固定)'
        else:
            verdict = ('(各结构不同，最大差 %.2f%% —— 专利近距 F 数表明镜头会收光圈时属正常；'
                       '没有孔径模型可对照)' % (100*spread))
        print('  自校验: 各结构 APER=%s → 光阑近轴半径 %s  %s'
              % ([round(aper.get(i, fnum[0]), 3) for i in range(1, len(rows)+1)],
                 [round(v, 4) for v in semis], verdict))


def cfg_paraxial(path):
    """反解析 .zmx，逐结构做近轴追迹（含 MCE THIC、位置解 TOLE、BINARY_2 衍射面）。

    返回每个结构 {'obj', 'ystop', 'u_img', 'beta', 'efl_inf'}：
      轴上物点出发的近轴边缘光线（物在 ∞ 时为 y=1 的平行光），ystop = 该光线在光阑面的高度，
      u_img = 像方斜率，beta = |u_物/u_像|（∞ 时 0）。
    ★ 衍射面的 2·C2·y 项不能漏 —— 漏了 DOE 镜头（RF600/800 F11）的光阑半径会反推错 20%。
    """
    txt = _read_zmx(path).splitlines()
    S = []; cur = None; mce = {}; names = []
    for ln in txt:
        t = ln.strip()
        if t.startswith('SURF '): cur = {'c': 0., 'd': 0., 'nd': None, 'stop': False, 'tole': None}; S.append(cur)
        elif cur is not None:
            if t.startswith('TYPE '): cur['type'] = t.split()[1]
            elif t.startswith('XDAT '):
                q = t.split(); cur.setdefault('xv', {})[int(q[1])] = float(q[2])
            elif t.startswith('CURV '): cur['c'] = float(t.split()[1])
            elif t.startswith('DISZ '):
                v = t.split()[1]; cur['d'] = float('inf') if v.upper() == 'INFINITY' else float(v)
            elif t.startswith('GLAS '):
                q = t.split()
                cur['nd'] = float(q[4]) + (float(q[10]) if len(q) > 10 and q[2] == '4' else 0.0)
            elif t.startswith('STOP'): cur['stop'] = True
            elif t.startswith('TOLE '): q = t.split(); cur['tole'] = (int(q[1]), float(q[2]))
        if t.startswith('LTTL'): names.append(t)
        if t.startswith('THIC'):
            q = t.split(); mce.setdefault(int(q[1]), {})[int(q[2])] = float(q[3])
    body = S[1:-1]
    ks = [i for i, s in enumerate(body) if s['stop']][0]
    def c2_of(s):
        if s.get('type') != 'BINARY_2': return 0.0
        xv = s.get('xv', {}); rn = xv.get(2, 1.0)
        return xv.get(3, 0.0) / (2*math.pi/(587.56e-6)) / rn**2      # 假定基准波长 587.56nm（与 EFL 自校验同）
    def trace(ds, y, u):
        n = 1.; ys = []
        for k, s in enumerate(body):
            n2 = s['nd'] or 1.
            u = (n*u - y*s['c']*(n2-n) + 2*c2_of(s)*y) / n2
            ys.append(y)
            if k < len(body)-1: y += u*ds[k]
            n = n2
        return ys, u
    ncfg = max([len(names)] + [len(v) for v in mce.values()] + [1])
    out = []
    for ci in range(1, ncfg+1):
        ds = [s['d'] for s in body]
        for row, v in mce.items():
            if row >= 1 and ci in v: ds[row-1] = v[ci]
        for j, s in enumerate(body):          # 位置解：d_j = L − Σ d(a..j-1)
            if s['tole']:
                a0, L0 = s['tole']; ds[j] = L0 - sum(ds[a0-1:j])
        o = mce.get(0, {}).get(ci, S[0]['d'])
        fin = o is not None and o < 1e9
        _ya, ua = trace(ds, 1.0, 0.0)
        if fin:
            ys, u = trace(ds, o, 1.0)
            beta = abs(1.0/u)
        else:
            ys, u = _ya, ua
            beta = 0.0
        out.append({'obj': o if fin else None, 'ystop': ys[ks], 'u_img': u, 'beta': beta,
                    'efl_inf': -1.0/ua})
    return out


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('spec'); ap.add_argument('-o',required=True)
    ap.add_argument('--emb', type=int, default=0)
    ap.add_argument('--asph-type', choices=('auto', 'even', 'extended'), default='auto',
                    help='auto=有 A18/A20 的面用 Extended Asphere(XASPHERE)、其余 Even Asphere；'
                         'extended=所有非球面都用 XASPHERE；'
                         'even=全部强按 Even Asphere（高次项会被重拟合进 r^16，有残差）')
    ap.add_argument('--no-merit', action='store_true',
                    help='不写配套评价函数、不把对焦间隔设成变量（默认写：Contrast s+t 80lp/mm GQ3×6 逐结构 + 对焦间隔变量）')
    ap.add_argument('--mf-freq', type=float, default=80.0, help='评价函数的空间频率 lp/mm（默认 80）')
    ap.add_argument('--modelglass', action='store_true',
                    help='另外再出一份模型玻璃版 <o>_modelglass.zmx（默认只出目录版）')
    a=ap.parse_args()
    spec=json.load(open(a.spec,encoding='utf-8'))
    emb=spec['embodiments'][a.emb]
    idx_of(emb['surfaces'])
    global _idx
    _idx = lambda s: s['_i']
    pf=dict(emb.get('general',[])).get('f (mm)', 0)
    # 交付默认只出目录玻璃版（用户 2026-09：「只需要输出一个 catalog zmx 和 seq」）；
    # 模型玻璃版要的时候加 --modelglass。
    for cat,tag in ((True,'catalog'),) + (((False,'modelglass'),) if a.modelglass else ()):
        p='%s_%s.zmx'%(a.o,tag)
        open(p,'wb').write(build(spec,emb,cat,a.asph_type,not a.no_merit,a.mf_freq))
        print(p); verify(p, pf, wfno_cfg(spec, emb, with_stop=True)[1])

if __name__=='__main__':
    main()
