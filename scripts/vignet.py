#!/usr/bin/env python3
"""子午实光线追迹 → 解渐晕系数 VDY/VCY，并找出真正卡光的「渐晕定义面」。

    python3 vignet.py spec.matched.json --write

依据：专利的「有効径」就是厂家设定的渐晕口径。用实光线（含非球面）在这些口径下
找出每个视场实际能通过的光瞳范围，换算成 Zemax 的 VDY/VCY；同时统计每条边缘
光线被哪个面挡住 —— 被挡次数最多的那几个面就是渐晕定义面，只需固定它们。
**逐结构求解**：对焦镜头每个结构的渐晕都不一样（适马70微距 ∞→1:1，
视场1 的 VDY 从 +0.12 一路走到 −0.34），Zemax 的 VDX/VDY/VCX/VCY 是全局量，
不进多重结构就等于所有结构共用一套 —— 必须写成 MCE 的 FVDY/FVCY/FVCX 操作数。
**解出来的光瞳还要保证真过得去**：Layout 画的是 Py=−1/0/+1 三条子午光线，
正是椭圆的 ±Y 端点；边缘卡在口径上（哪怕只差 0.05mm）那两条就会被挡掉。
纯标准库。**含斜光线**：VDY/VCY 由子午扫描解，VCX 由 3D 斜光线在光瞳中心高度上
横向扫描解。只解子午、把 VCX 留 0 是错的 —— Zemax 会按全宽 X 光瞳追迹，
光阑后各面的自动口径必然被撑大，离轴视场的光束也会比实际大一圈。
"""
import re
import json, math, argparse

# --- 非球面系数：支持奇数次项（佳能 A3..A15 等）-------------------------------
# 返回 [(指数, 系数), ...]，指数为任意正整数。偶数次专利仍然原样工作。
def acoef(a):
    out = []
    for k, v in (a or {}).items():
        m = re.fullmatch(r'[Aa]\s*(\d+)', str(k))
        if m and v:
            out.append((int(m.group(1)), float(v)))
    return sorted(out)
# -----------------------------------------------------------------------------


# --- 衍射面（DOE，佳能 DO / 尼康 PF）------------------------------------------
# 专利式：ψ(h) = 2πm/λ0 · Σ C_{2j} h^{2j}  ⇒  W(h)=Σ C h^{2j} 是「基准波长下的 OPD（mm）」。
# 追迹用 d 线折射率 ≈ 基准波长，m=1：切向动量增量 = dW/dr · r̂（投到面切平面）。
# 近轴光焦度 φ_DOE = −2·C2。spec 写法：面上 "doe": {"wl_nm": 587.56, "C": [C2, C4, C6, …]}
def doe_dw(C, r):
    return sum(2 * (j + 1) * c * r ** (2 * j + 1) for j, c in enumerate(C or []))

def _doe_kick(p, N, g):
    """p = n'·d'（折射后），N = 单位法线，g = 衍射切向增量（3 维）。返回新的单位方向。"""
    n2 = sum(v * v for v in p) ** 0.5
    pn = sum(a * b for a, b in zip(p, N))
    gn = sum(a * b for a, b in zip(g, N))
    pt = [a - pn * b + (c - gn * b) for a, b, c in zip(p, N, g)]
    q = n2 * n2 - sum(v * v for v in pt)
    if q < 0: return None
    s = math.sqrt(q) * (1 if pn >= 0 else -1)
    return [(a + s * b) / n2 for a, b in zip(pt, N)]
# -----------------------------------------------------------------------------


class Surf:
    def __init__(s, c, k, A, z, n_after, semi, is_stop, doe=None):
        s.c, s.k, s.A, s.z, s.n, s.semi, s.stop = c, k, A, z, n_after, semi, is_stop
        s.doe = doe
    def sag(s, y):
        y2 = y * y
        # r^18 / r^20 项在追迹发散时会把 y 顶到 1e15，y**20 直接 OverflowError。
        # 任何真实镜片都不会有 |y| > 1000mm，越界一律当追失。
        if not (y2 < 1e6): return None
        r = 1 - (1 + s.k) * s.c * s.c * y2
        if r < 0: return None
        z = s.c * y2 / (1 + math.sqrt(r))
        ay = abs(y)
        # 系数表已是 [(指数, 系数)]，奇数次项必须用 |y| —— 面型只依赖 r=|y|。
        for e, cc in s.A: z += cc * ay ** e
        return z
    def dsag(s, y):
        h = 1e-7
        a, b = s.sag(y + h), s.sag(y - h)
        return None if a is None or b is None else (a - b) / (2 * h)

def build(spec, emb, state, dmap=None):
    asph = {str(a['surface']).replace('面', ''): a for a in emb.get('aspheric', [])}
    S, z = [], 0.0
    rows = [s for s in emb['surfaces'] if s['i'] != 'IMG']
    img = [s for s in emb['surfaces'] if s['i'] == 'IMG']
    for s in rows:
        D = s['D']
        if isinstance(D, str): D = (dmap or {}).get(D, emb['variable'][D][state])
        a = asph.get(str(s['i']))
        A = acoef(a) if a else []
        phi = (s.get('extra') or {}).get('有効径 φi')
        S.append(Surf(0.0 if s['R'] in (None, 0) else 1.0/float(s['R']),
                      (a or {}).get('k', (a or {}).get('K', 0.0)) or 0.0, A, z,
                      s.get('nd') or 1.0, (phi/2.0) if phi else None,
                      bool(s.get('stop') or s['i'] == 'STO'),
                      (s.get('doe') or {}).get('C')))
        z += float(D)
    return S, z, ((img[0].get('extra') or {}).get('有効径 φi') if img else None)

def cfg_dmap(spec, emb, c):
    """一个结构的**全部**可变间隔。

    ★ 不能只取 key_before：configs 里只存了对焦群前面那个间隔，key_after（= 守恒和 − key_before）
    与链式第三段 key_last 不补上的话，build() 会回落到基准态的值 —— 近距结构的镜头平白变长，
    像面跟着跑掉（实测 JP2021-047297A MFD 结构总长 212.01，真值 162.37），渐晕全建立在错的几何上。
    """
    zx = spec.get('zmx') or {}
    var = emb.get('variable') or {}
    d = {}
    for k, v in c.items():
        if k not in var or isinstance(v, bool): continue
        try: d[k] = float(v)                    # 老 spec 里偶有字符串数值（'2.603'）
        except (TypeError, ValueError): pass
    fc, fc2 = zx.get('focus') or {}, zx.get('focus2') or {}
    kb1, kb2, kl = fc.get('key_before'), fc2.get('key_before'), fc.get('key_last')
    if kl and kb1 in d and kb2 in d and 'sum' in fc:
        d.setdefault(kl, float(fc['sum']) - d[kb1] - d[kb2])
    else:
        for f, kb in ((fc, kb1), (fc2, kb2)):
            ka = f.get('key_after')
            if ka and kb in d and 'sum' in f:
                d.setdefault(ka, float(f['sum']) - d[kb])
    return d


def _par(S, y0, u0):
    """近轴追迹（与 solve_state 里的 par 同一套式子）→ (各面高度, 像方斜率)。"""
    n, y, u, ys = 1.0, y0, u0, []
    for k, s in enumerate(S):
        u = (n*u - y*s.c*(s.n - n) + 2.0*((s.doe or [0.0])[0])*y) / s.n
        ys.append(y)
        if k < len(S)-1: y += u * (S[k+1].z - S[k].z)
        n = s.n
    return ys, u


def _lagr(P, t):
    r = 0.0
    for i, (xi, yi) in enumerate(P):
        term = yi
        for j, (xj, _) in enumerate(P):
            if i != j: term *= (t - xj) / (xi - xj)
        r += term
    return r


def aperture_cfg(spec, emb, state, cfgs):
    """逐结构近轴工作 F 数（Zemax: Paraxial Working F/#，逐结构 APER）+ 对应的入瞳半径。

    两层模型，按优先级：

    ① **专利印了各对焦态的 F 数（zmx.fno_patent = {结构名: F}）就以它为准**，中间结构按 |β| 插值。
       实证 JP2021-047297A（RF100mm F2.8L Macro）：专利 0.5x/1.4x 印 4.49/6.64，
       佳能官方说明书「撮影倍率と実効FNo.」表 0.5x/1.0x/1.4x = 4.5/5.7/6.6，
       the-digital-picture 实测光损 1x 2 档、1.4x 2⅓ 档 —— 三方一致；
       而「光阑全开不变」只有 3.46/4.13/4.80。说明书还写着「对焦时光圈叶片会动」：
       这支镜头**近距会主动收光圈**，固定光阑模型在近距是错的。
       插值：落在专利记载 |β| 区间内的结构，对 F 数做 n 点 n−1 次拉格朗日（不单调就退成分段线性）；
       区间外按最近那一态的「光阑收缩比」外推（光阑不再继续变）。
       本篇验证：插值给出 1.00x = 5.79，说明书 5.7（+1.6%）；0.3x = 3.89，实测 f/4。
       守门：印刷值比「光阑全开不变」的工作 F 数还亮（光阑只能比 ∞ 更大才做得到）→ 不可信，丢弃该点；
       各态印的都等于 ∞ 值（名义 F 数，不是工作 F 数）→ 整组丢弃。
    ② 否则**物理光阑固定**：在 ∞ 结构按专利 F 数（或 zmx.epd）定出光阑近轴半径，
       逐结构追近轴边缘光线求像方斜率 u'：WFNO = 1/(2|n'u'|)。∞ 结构恰好等于专利 F 数；
       整组对焦、内对焦、浮动对焦、光阑前后有无移动组，一套式子全覆盖。

    ★ 为什么绝不能用「像方空间 F/#」（Image Space F/#）：它按**该结构自己**的 ∞ 共轭 EFL / 入瞳直径定义，
      内对焦/浮动对焦近距 EFL 大幅缩短（本篇 100.81 → 46.2 @1x → 36.0 @1.4x），光阑被跟着缩小 ——
      用户实测 1.00x 结构 ENPD 15.73、状态栏 WFNO 9.149（说明书 5.7）。

    zmx.wfno_override = {结构名: F} 最高优先级，逐结构强行指定。

    ★ 变焦镜头（zmx.zoom + 结构带 'zoom' 键）：**按变焦位置分组各算一遍**。
      每个变焦位置的光阑在它自己的 ∞ 结构上、按该位置的 F 数（zoom.positions[].fno）定 ——
      恒定光圈变焦的物理光阑直径随焦距变，全变焦共用一个光阑是错的；
      同一位置的对焦结构再按「光阑固定」（或 fno_patent）算工作 F 数。
    """
    zx = spec['zmx']
    zpos = {p['name']: p for p in ((zx.get('zoom') or {}).get('positions') or [])}
    if zpos and any(c.get('zoom') for c in cfgs):
        order = []
        for c in cfgs:
            if c.get('zoom') not in order: order.append(c.get('zoom'))
        rows_all, notes_all, s_first = [None] * len(cfgs), [], None
        names_all = [c.get('name') for c in cfgs]
        ov_all, fp_all = zx.get('wfno_override') or {}, zx.get('fno_patent') or {}
        for tag, tab in (('fno_patent', fp_all), ('wfno_override', ov_all)):
            for k in tab:
                if k not in names_all:
                    notes_all.append('%s 的键 %r 对不上任何结构名 —— 已忽略' % (tag, k))
        for zn in order:
            idx = [i for i, c in enumerate(cfgs) if c.get('zoom') == zn]
            p = zpos.get(zn)
            if p is None:
                raise SystemExit('结构 %s 的 zoom=%r 在 zmx.zoom.positions 里找不到'
                                 % ([cfgs[i].get('name') for i in idx], zn))
            gnames = {cfgs[i].get('name') for i in idx}
            zx2 = {k: v for k, v in zx.items() if k not in ('zoom', 'epd')}
            # 每组只拿自己结构名的键，否则别的组会报「对不上任何结构名」的假告警
            zx2['wfno_override'] = {k: v for k, v in ov_all.items() if k in gnames}
            zx2['fno_patent'] = {k: v for k, v in fp_all.items() if k in gnames}
            if p.get('fno'): zx2['fno'] = float(p['fno'])
            if p.get('epd'): zx2['epd'] = float(p['epd'])
            spec2 = dict(spec); spec2['zmx'] = zx2
            ss, rows = aperture_cfg(spec2, emb, state, [cfgs[i] for i in idx])
            notes_all += ['[%s] %s' % (zn, t) for t in getattr(aperture_cfg, 'notes', [])]
            if s_first is None: s_first = ss
            for i, r in zip(idx, rows):
                r['zoom'] = zn; r['stop_semi_inf'] = ss; rows_all[i] = r
        aperture_cfg.notes = notes_all
        return s_first, rows_all
    ov = zx.get('wfno_override') or {}
    fpat = {k: float(v) for k, v in (zx.get('fno_patent') or {}).items()}
    notes = []
    names = [c.get('name') for c in cfgs]
    for k in list(fpat):
        if k not in names:
            notes.append('fno_patent 的键 %r 对不上任何结构名 %s —— 已忽略' % (k, names))
    for k in list(ov):
        if k not in names:
            notes.append('wfno_override 的键 %r 对不上任何结构名 —— 已忽略' % k)
    def _is_inf(c):
        d0 = c.get('d0')
        if d0 is None or str(d0).upper().startswith('INF'): return True
        try: return float(d0) >= 1e9
        except (TypeError, ValueError): return False
    # ★ 光阑必须在 **∞ 结构**上按专利 F 数定（专利 F 数就是 ∞ 的值）；结构顺序不一定是 ∞ 在前
    c_inf = next((c for c in cfgs if _is_inf(c)), None)
    if c_inf is None:
        c_inf = cfgs[0]
        notes.append('结构里没有 ∞ 物距 —— 光阑按结构 %r 的 F 数定，专利 F 数是 ∞ 值时这会偏' % c_inf.get('name'))
    S0, _z0, _ = build(spec, emb, state, cfg_dmap(spec, emb, c_inf))
    ks0 = [k for k, s in enumerate(S0) if s.stop][0]
    yA0, uA0 = _par(S0, 1.0, 0.0)
    rEP0 = 0.5*float(zx['epd']) if zx.get('epd') else (1.0/(2*zx['fno'])) / abs(uA0)
    stop_semi = rEP0 * abs(yA0[ks0])
    rows = []
    for c in cfgs:
        S, _zimg, _ = build(spec, emb, state, cfg_dmap(spec, emb, c))
        ks = [k for k, s in enumerate(S) if s.stop][0]
        yA, uA = _par(S, 1.0, 0.0)
        obj = None if _is_inf(c) else float(c.get('d0'))
        if obj is None:
            ys, u, beta = yA, uA, 0.0
        else:
            ys, u = _par(S, obj, 1.0)          # 轴上物点出发、斜率 1 的近轴边缘光线
            beta = abs(1.0 / u)                # 近轴横向放大率 |β| = |u_物 / u_像|（物方 u=1、空气）
        sc = stop_semi / abs(ys[ks])
        wf_fixed = 1.0 / (2.0*abs(u*sc))
        rows.append({'name': c.get('name'), 'efl': -1.0/uA, 'beta': beta, 'wfno_fixed': wf_fixed,
                     'yA_stop': abs(yA[ks]), 'obj': obj, 'src': '固定光阑'})
    # ---- ① 专利印刷的各态 F 数 ----
    raw = [(r['beta'], fpat[r['name']], r['wfno_fixed'], r['name']) for r in rows if r['name'] in fpat]
    close_raw = [p for p in raw if p[0] > 1e-6]
    f_inf = 0.5/(rEP0*abs(uA0)) if not zx.get('epd') else None     # = zx.fno
    f_ref = float(zx['fno']) if zx.get('fno') else (raw[0][1] if raw else None)
    pts = []
    # 名义 F 数判定要**先**做：近距各态印的都等于 ∞ 的 F 数、而光阑固定时应当明显更暗 → 整组是名义值
    if close_raw and f_ref and all(abs(p[1] - f_ref) < 0.005*f_ref for p in close_raw) \
            and any(p[2] > 1.02*p[1] for p in close_raw):
        notes.append('近距各态印的 F 数都等于 ∞ 值 %.2f → 是名义 F 数（不是工作 F 数），整组丢弃' % f_ref)
    else:
        for b_, f0, wfx, nm in raw:
            if f0 < wfx * 0.995:
                notes.append('%s 印 F%.2f 比光阑全开不变的 %.3f 还亮 → 不可信，丢弃' % (nm, f0, wfx))
                continue
            pts.append((b_, f0, wfx))
    # ∞ 锚点：专利 F 数本来就是 ∞ 的值；fno_patent 里没写 INF 键时补上，
    # 否则 |β| 比最小记载态还小的结构会按「收缩比外推」把 ∞ 的 F 数也改掉（回归实测 INF 2.92 → 3.79）
    if pts and not any(p[0] < 1e-9 for p in pts):
        r_inf = next((r for r in rows if r['beta'] < 1e-9), None)
        if r_inf is not None:
            pts.append((0.0, r_inf['wfno_fixed'], r_inf['wfno_fixed']))
    pts.sort()
    if len(pts) >= 2:
        P = [(b, f) for b, f, _ in pts]
        mono = all(_lagr(P, pts[0][0] + (pts[-1][0]-pts[0][0])*i/200) <=
                   _lagr(P, pts[0][0] + (pts[-1][0]-pts[0][0])*(i+1)/200) + 1e-9 for i in range(200))
        for r in rows:
            b = r['beta']
            if pts[0][0] - 1e-9 <= b <= pts[-1][0] + 1e-9:
                if mono:
                    ft = _lagr(P, b)
                else:
                    ft = next(f0 + (f1-f0)*(b-b0)/(b1-b0) for (b0, f0), (b1, f1) in zip(P, P[1:])
                              if b0 - 1e-9 <= b <= b1 + 1e-9)
                src = '专利印刷' if any(abs(b - q[0]) < 1e-6 for q in pts) else \
                      ('专利插值' if mono else '专利分段线性')
            else:
                q = pts[-1] if b > pts[-1][0] else pts[0]
                ft = r['wfno_fixed'] * q[1] / q[2]                    # 光阑收缩比保持最近那一态
                src = '专利外推(光阑比)'
            if ft >= r['wfno_fixed'] - 1e-9:
                r['wfno_patent'] = ft; r['src'] = src
    for r in rows:
        wf = float(ov.get(r['name'], r.get('wfno_patent', r['wfno_fixed'])))
        if r['name'] in ov: r['src'] = 'wfno_override'
        r['wfno'] = wf
        r['stop_semi'] = stop_semi * r['wfno_fixed'] / wf             # 近轴线性：光阑半径 ∝ 1/WFNO
        r['rEP'] = r['stop_semi'] / r['yA_stop']
    aperture_cfg.notes = notes
    return stop_semi, rows


def _window(y, dy, ylim, t0):
    """|y + t·dy| ≤ ylim 对应的 t 区间（面型多项式只在净口径附近才有意义）。"""
    if ylim >= 9e3 or abs(dy) < 1e-12: return t0 - 400.0, t0 + 400.0
    ta, tb = (-ylim - y) / dy, (ylim - y) / dy
    return (ta, tb) if ta < tb else (tb, ta)


def _window3(x, y, dx, dy, ylim, t0):
    """3D：r(t)² = (x+t·dx)² + (y+t·dy)² ≤ ylim² 对应的 t 区间。"""
    if ylim >= 9e3: return t0 - 400.0, t0 + 400.0
    a = dx*dx + dy*dy
    if a < 1e-18: return t0 - 400.0, t0 + 400.0
    b = 2.0*(x*dx + y*dy); c = x*x + y*y - ylim*ylim
    disc = b*b - 4*a*c
    if disc <= 0: return t0, t0          # 整条光线都在有效域外
    sq = math.sqrt(disc)
    return (-b - sq)/(2*a), (-b + sq)/(2*a)


def _solve_t(f, fp, t0, ta, tb):
    """求 t 使 f(t)=0（光线与面的交点），只在 [ta,tb] 这个「面型有意义」的窗口里找。

    2026-09 修：原来用**定点迭代** t ← (顶点z + sag(y(t)) − z)/dz，收敛条件是
    |tan(入射角)·sag′| < 1。超广角（半角 56°）打在强曲率/强非球面的前组上时这个乘子
    >2，迭代**发散或跳到伪根** —— 表现为同一视场的光瞳里「过/挡」交替成碎片
    （实测 US20260235852A1 Ex1：面1 高度在 Py=−1 处给出 −20.5，Py=−0.5 处 −27.9，
    非单调，纯属伪根），渐晕随后被解成「最长连续可行段」里的某个碎片，VCY 0.88（假）。
    改成 Newton（f′ = dz − dy·sag′，没有那个收敛限制）+ 窗口内扫描二分兜底；
    窗口 |y| ≤ 1.05·净口径 把非球面多项式在有效域外的发散根整个排除掉。
    """
    if tb <= ta: return None
    t = min(max(t0, ta), tb)
    last = t
    for _ in range(60):
        v = f(t)
        if v is None:
            t = 0.5 * (t + last)
            if abs(t - last) < 1e-12: break
            continue
        if abs(v) < 1e-11: return t
        d = fp(t)
        if d is None or abs(d) < 1e-9: break
        last = t
        t = min(max(t - v / d, ta), tb)
        if abs(t - last) < 1e-13: return t if abs(f(t) or 1) < 1e-8 else None
    pt, pv = None, None
    for i in range(121):
        tt = ta + (tb - ta) * i / 120.0
        vv = f(tt)
        if vv is not None and pv is not None and vv * pv <= 0:
            a, b, fa = pt, tt, pv
            for _ in range(70):
                m = 0.5 * (a + b); fm = f(m)
                if fm is None: break
                if fa * fm <= 0: b = m
                else: a, fa = m, fm
            return 0.5 * (a + b)
        pt, pv = tt, vv
    return None


def trace(S, zimg, y0, z0, ang, apert=True):
    """从 (y0,z0) 以角 ang(rad) 出发追子午实光线。返回 (各面高度, 被挡的面序号或None, 像高)。"""
    dy, dz = math.sin(ang), math.cos(ang)
    y, z, n = y0, z0, 1.0
    hs = []
    for k, s in enumerate(S):
        t0 = (s.z - z) / dz                      # 先交到顶点平面
        def _f(tt, _s=s, _y=y, _z=z, _dy=dy, _dz=dz):
            sg = _s.sag(_y + tt * _dy)
            return None if sg is None else (_z + tt * _dz) - (_s.z + sg)
        def _fp(tt, _s=s, _y=y, _dy=dy, _dz=dz):
            sl = _s.dsag(_y + tt * _dy)
            return None if sl is None else _dz - _dy * sl
        _yl = (1.05 * s.semi + 0.2) if s.semi is not None else 1e4
        _ta, _tb = _window(y, dy, _yl, t0)
        t = _solve_t(_f, _fp, t0, _ta, _tb)
        if t is None: return hs, k, None
        y, z = y + t * dy, z + t * dz
        hs.append(y)
        if apert and s.semi is not None and abs(y) > s.semi + 1e-9:
            return hs, k, None
        sl = s.dsag(y)
        if sl is None: return hs, k, None
        nx, nz = -sl, 1.0
        L = math.hypot(nx, nz); nx, nz = nx/L, nz/L
        if dy*nx + dz*nz > 0: nx, nz = -nx, -nz     # 法线取与入射方向相反
        mu = n / s.n
        ci = -(dy*nx + dz*nz)
        disc = 1 - mu*mu*(1 - ci*ci)
        if disc < 0: return hs, k, None
        ct = math.sqrt(disc)
        dy, dz = mu*dy + (mu*ci - ct)*nx, mu*dz + (mu*ci - ct)*nz
        if s.doe:
            gw = doe_dw(s.doe, abs(y)) * (1 if y >= 0 else -1)
            dd = _doe_kick([0.0, s.n*dy, s.n*dz], [0.0, nx, nz], [0.0, gw, 0.0])
            if dd is None: return hs, k, None
            dy, dz = dd[1], dd[2]
        n = s.n
    t = (zimg - z) / dz
    return hs, None, y + t * dy

def trace3(S, zimg, P0, d, apert=True):
    """3D 斜光线。旋转对称系统：面型只依赖 r=hypot(x,y)，法线 = (-sl*x/r, -sl*y/r, 1)。
    返回 (各面处的 r, 被挡面序号或 None, 像面 y, 像面 x)。"""
    x, y, z = P0; dx, dy, dz = d
    n = 1.0; rs = []
    for k, s in enumerate(S):
        t0 = (s.z - z) / dz
        def _f3(tt, _s=s, _x=x, _y=y, _z=z, _dx=dx, _dy=dy, _dz=dz):
            sg = _s.sag(math.hypot(_x + tt*_dx, _y + tt*_dy))
            return None if sg is None else (_z + tt*_dz) - (_s.z + sg)
        def _fp3(tt, _s=s, _x=x, _y=y, _dx=dx, _dy=dy, _dz=dz):
            xx, yy = _x + tt*_dx, _y + tt*_dy
            r = math.hypot(xx, yy)
            sl = _s.dsag(r) if r > 1e-9 else 0.0
            if sl is None: return None
            drdt = ((xx*_dx + yy*_dy) / r) if r > 1e-9 else 0.0
            return _dz - sl * drdt
        _yl = (1.05 * s.semi + 0.2) if s.semi is not None else 1e4
        _ta, _tb = _window3(x, y, dx, dy, _yl, t0)
        t = _solve_t(_f3, _fp3, t0, _ta, _tb)
        if t is None: return rs, k, None, None
        x += t*dx; y += t*dy; z += t*dz
        r = math.hypot(x, y); rs.append(r)
        if apert and s.semi is not None and r > s.semi + 1e-9:
            return rs, k, None, None
        sl = s.dsag(r) if r > 1e-9 else 0.0
        if sl is None: return rs, k, None, None
        if r > 1e-9: nx, ny, nz = -sl*x/r, -sl*y/r, 1.0
        else:        nx, ny, nz = 0.0, 0.0, 1.0
        L = math.sqrt(nx*nx + ny*ny + nz*nz); nx, ny, nz = nx/L, ny/L, nz/L
        if dx*nx + dy*ny + dz*nz > 0: nx, ny, nz = -nx, -ny, -nz
        mu = n / s.n
        ci = -(dx*nx + dy*ny + dz*nz)
        disc = 1 - mu*mu*(1 - ci*ci)
        if disc < 0: return rs, k, None, None
        ct = math.sqrt(disc)
        dx, dy, dz = (mu*dx + (mu*ci - ct)*nx,
                      mu*dy + (mu*ci - ct)*ny,
                      mu*dz + (mu*ci - ct)*nz)
        if s.doe and r > 1e-12:
            gw = doe_dw(s.doe, r)
            dd = _doe_kick([s.n*dx, s.n*dy, s.n*dz], [nx, ny, nz], [gw*x/r, gw*y/r, 0.0])
            if dd is None: return rs, k, None, None
            dx, dy, dz = dd
        n = s.n
    t = (zimg - z) / dz
    return rs, None, y + t*dy, x + t*dx




FRACS = (1.0, 0.8, 0.6, 0.4, 0.2, 0.0)
Z0 = -80.0                      # 无限远物：光线起始平面


def solve_state(S, zimg, zx, obj=None, margin=0.010, verbose=True, tag='', fit_ellipse=False, rEP_fixed=None):
    """解一个结构（一组间隔 + 一个物距）的渐晕。

    obj=None 表示物在无限远，视场参数 p = 入射半角；否则 obj 是物距（mm，面1顶点起算），
    视场参数 p = 物高。两种情形都用「入瞳高度 ye」参数化光瞳（Py = (ye-yep)/rEP）。

    margin：光瞳两端各内缩 margin×rEP 再折成 VDY/VCY。Layout 的 Py=±1 光线正落在端点上，
    不留余量时舍入误差 + Zemax 的实际光线瞄准（RAIM Real 打到真实光阑，不是近轴入瞳）
    就足以让那两条被口径切掉 —— 用户实测就是这么被挡的。
    """
    def par(y0, u0):
        n, y, u, ys = 1.0, y0, u0, []
        for k, s in enumerate(S):
            u = (n*u - y*s.c*(s.n - n) + 2.0*((s.doe or [0.0])[0])*y) / s.n   # DOE: φ=−2·C2
            ys.append(y)
            if k < len(S)-1: y += u * (S[k+1].z - S[k].z)
            n = s.n
        return ys, u
    yA, uA = par(1.0, 0.0)
    yB, _uB = par(0.0, 1.0)
    ks = [k for k, s in enumerate(S) if s.stop][0]
    rEP = (1.0/(2*zx['fno'])) / abs(uA)                 # 入瞳半径（随结构变：光阑前的 d12 在动）
    if zx.get('epd'):
        # 光阑在对焦组之前（EP4215968A1）：物理光阑不变 ⇒ 入瞳直径各结构恒定。
        # 按 F 数定的话，内对焦把近距 EFL 缩短，光阑会被跟着缩小（MFD 入瞳 66→24）。
        rEP = 0.5*float(zx['epd'])
    if rEP_fixed:
        # aperture_cfg() 按「物理光阑固定」算好的该结构入瞳半径 —— 与 .zmx 的
        # Paraxial Working F/#（逐结构 APER）是同一个模型，渐晕才对得上。
        rEP = float(rEP_fixed)
        # 光阑面本身不参与挡光：它的大小已由孔径定义给出（Zemax 开 Real 瞄准时 Py=±1 正打在
        # 光阑近轴半径上）。否则近距结构瞄近轴入瞳的实光线在光阑上略超专利有効径（本篇 15.51>15.405），
        # 会凭空在轴上写出 VCY≈0.017 的假渐晕。
        for _s in S:
            if _s.stop: _s.semi = None
    ca, cb = -yB[ks], yA[ks]
    y1 = ca*yA[0] + cb*yB[0]; u1 = cb
    zEP = -y1/u1 if u1 else 0.0
    z_obj = None if obj is None else -abs(obj)

    def start(p, ye):
        if obj is None:
            return ye + (Z0 - zEP)*math.tan(p), Z0, p
        return p, z_obj, math.atan2(ye - p, zEP - z_obj)

    def start3(p, px_mm, ye):
        if obj is None:
            return ((px_mm, ye + (Z0 - zEP)*math.tan(p), Z0),
                    (0.0, math.sin(p), math.cos(p)))
        vx, vy, vz = px_mm, ye - p, zEP - z_obj
        L = math.sqrt(vx*vx + vy*vy + vz*vz)
        return (0.0, p, z_obj), (vx/L, vy/L, vz/L)

    def shoot(p, ye, apert=True):
        y0, z0, ang = start(p, ye)
        return trace(S, zimg, y0, z0, ang, apert=apert)

    pupil_scale = 1.0
    r_sp = rEP * abs(yA[ks])
    if rEP_fixed:
        # ★ 与 Zemax「Ray Aiming = Real」对齐：Py=±1 打在光阑的**近轴半径**上，不是近轴入瞳边缘。
        # 大孔径时入瞳有球差，瞄近轴入瞳边缘的实光线在光阑上会多出 1~2%
        # （回归实测 US20240302626A1 RF35 F1.46：13.65 vs 13.348），各面轴上需求跟着虚胖，
        # 凭空写出轴上渐晕（VCY 0.013/0.024/0.057）和一串「★★ 切到轴上光瞳」假报警。
        # 做法：轴上视场二分出「实光线正好落在光阑近轴半径上」的入瞳高度，拿它当本结构的 rEP。
        # （离轴的光瞳彗差这里不管 —— 交给 zapi_vigfit.ps1 用 OpticStudio 真追迹兜底。）
        r_sp = rEP * abs(yA[ks])
        def _hstop(ye):
            hs, _b, _yi = shoot(0.0, ye, apert=False)
            return (abs(hs[ks]) - r_sp) if (hs and len(hs) > ks and _yi is not None) else None
        lo, hi = 0.7*rEP, 1.3*rEP
        flo, fhi = _hstop(lo), _hstop(hi)
        for _ in range(40):                     # 有界：hi 逼近 rEP 时浮点会卡死（回归实测 >2 万次追迹）
            if fhi is not None or hi <= rEP*(1 + 1e-9): break
            hi = 0.5*(hi + rEP); fhi = _hstop(hi)
        if flo is not None and fhi is not None and flo < 0 < fhi:
            for _ in range(50):
                mid = 0.5*(lo + hi); fm = _hstop(mid)
                if fm is None or fm > 0: hi = mid
                else: lo = mid
            pupil_scale = 0.5*(lo + hi) / rEP
            rEP = 0.5*(lo + hi)

    def shoot3(p, px_mm, ye):
        P0, d = start3(p, px_mm, ye)
        return trace3(S, zimg, P0, d)

    def shoot3n(p, px_mm, ye):
        """不受口径阻挡的追迹 —— 用来算「这条光线本来该走到多高」。"""
        P0, d = start3(p, px_mm, ye)
        return trace3(S, zimg, P0, d, apert=False)

    def chief(p):
        """解过光阑中心的实主光线，返回 (入瞳高度 yep, 像高)。"""
        def val(ye):
            hs, _b, yi = shoot(p, ye, apert=False)
            return (hs[ks] if len(hs) > ks else None), yi
        # 取「离入瞳中心最近」的那个根，不是从最负端数过来的第一个。
        # 大孔径强像差系统里 h_stop(ye) 并不单调，远端会出现**伪根**；旧写法一撞上
        # 伪根就返回一个荒唐的像高，外层按像高二分视场角的 bracket 随之塌掉 ——
        # 实测本篇 INF 的 12.98/8.65 两个视场都被钉在 8.69°（真值 14.9°/10.1°）。
        # 真主光线按定义过近轴入瞳中心，实光线只偏一点点，所以 |ye| 最小的根才是它。
        for span in (1.5, 4.0, 12.0):
            prev_m = prev_v = None; cands = []
            N = 96
            for i in range(N+1):
                ye = -span*rEP + 2*span*rEP*i/N
                v, _ = val(ye)
                if v is None: prev_m = prev_v = None; continue
                if prev_v is not None and (v == 0 or (v > 0) != (prev_v > 0)):
                    cands.append((abs(0.5*(prev_m+ye)), prev_m, ye))
                prev_m, prev_v = ye, v
            lo = hi = None
            if cands:
                cands.sort(); lo, hi = cands[0][1], cands[0][2]
                break
        if lo is None: return None, None
        for _ in range(44):
            m = 0.5*(lo+hi); v, _ = val(m)
            if v is None: break
            if v > 0: hi = m
            else: lo = m
        m = 0.5*(lo+hi); _v, yi = val(m)
        return m, yi

    def chief2(p):
        """(参考入瞳高度, 像高, 是否真主光线)。

        近距超广角的最大视场会出现**主光线物理上追不出来**的情形：过光阑中心那条
        要在前组上走到 |y|>|R|（面型压根不存在那一圈），于是整片退成 vig 全 0（满光瞳），
        Zemax 会据此把上游口径撑到造不出来 —— 这正是红线里最忌讳的兜底。
        实测 US20250389929A1 実施例2 的 MFD(200mm) 结构视场1 就是这样。
        退路：参考点仍取**近轴入瞳中心 ye=0**（Zemax 的 Py=0 就定义在这里），
        像高改由「可行光瞳区间中点」那条光线给出。VDY/VCY 随后仍在 yep±rEP 内解，
        这样 |VDY|+(1−VCY) ≤ 1 自动成立，光瞳不会被推出单位圆。
        """
        ye0, yi0 = chief(p)
        if ye0 is not None: return ye0, yi0, True
        NS, span = 240, 3.0
        okv = [shoot(p, -span*rEP + 2*span*rEP*i/NS)[1] is None for i in range(NS+1)]
        best, i = (0, -1, -1), 0
        while i <= NS:
            if okv[i]:
                j = i
                while j + 1 <= NS and okv[j+1]: j += 1
                if j - i > best[0]: best = (j - i, i, j)
                i = j + 1
            else: i += 1
        if best[2] < 0: return None, None, False
        ymid = -span*rEP + 2*span*rEP*(0.5*(best[1]+best[2]))/NS
        return 0.0, shoot(p, ymid, apert=False)[2], False

    def edge(p, yep, sign):
        """从光瞳中心往 sign 方向二分出「最后一条过得去的」入瞳高度（精确，不吃网格步长）。"""
        ok, bad = yep, yep + sign*rEP
        if shoot(p, bad)[1] is None: return bad
        for _ in range(44):
            m = 0.5*(ok+bad)
            if shoot(p, m)[1] is None: ok = m
            else: bad = m
        return ok

    def ellipse_blocked(p, yc, ax, ay, n=24):
        bad = 0
        for i in range(n):
            th = 2*math.pi*i/n
            if shoot3(p, ax*math.cos(th), yc + ay*math.sin(th))[1] is not None: bad += 1
        return bad

    ymax = zx['max_y']
    vig, beams, blk_all = [], [], {}
    maxr = [0.0]*len(S)
    # ---- 轴上满光瞳的逐面需求半径（不加任何渐晕、不受口径阻挡）----
    # 用户 2026-09 打回：口径要是切到轴上大光瞳，**F 数就变了**（光圈变小），
    # 那不是"少一点渐晕"的问题，是整只镜头的相对孔径被改掉。所以轴上这一圈是硬底线。
    axr = [0.0]*len(S)
    for _i in range(48):
        _th = 2*math.pi*_i/48
        for _g in (1.0, 0.85, 0.6):
            _rr = shoot3n(0.0, _g*rEP*math.cos(_th), _g*rEP*math.sin(_th))[0]
            for _k, _v in enumerate(_rr): axr[_k] = max(axr[_k], _v)
    if verbose:
        print('%s入瞳半径 %.4f  入瞳位置 面1前 %.3f%s%s'
              % (tag, rEP, -zEP, '' if obj is None else '   物距 %.2f' % obj,
                 '' if abs(pupil_scale - 1) < 1e-6 else
                 '   实光线瞄准缩放 ×%.4f（轴上边缘光线落在光阑近轴半径上）' % pupil_scale))
        print('  视场 Y′  实际半角      VDY      VCY      VCX   通过光瞳  Py/Px±1  主要卡光面')
    for fr in FRACS:
        Y = ymax*fr
        # --- 解视场参数 p：无限远解半角，有限共轭解物高 ---
        if Y == 0:
            p = 0.0
        elif obj is None:
            lo, hi = 1e-4, math.radians(70)
            for _ in range(34):
                m = 0.5*(lo+hi); _, yi, _r = chief2(m)
                if yi is None: hi = m
                elif abs(yi) < Y: lo = m
                else: hi = m
            p = 0.5*(lo+hi)
        else:
            # 有限共轭解物高：先用小物高探出线性放大率给出量级 p0，再在 [0, 1.6*p0] 二分。
            # 盲目从 Y 起 ×1.6 外推的老写法，一旦物高大到主光线追失就会把 hi 推到天上，
            # 最大视场整片解不出、退化成 vig 全 0（= 满光瞳）——实测 JP2025052870A
            # 的 0.06x 结构视场1 就是这样，害得 Zemax 把面9~12 的自动口径撑到 31.7。
            p_probe, yi0 = max(1.0, Y)*0.05, None
            for _ in range(24):
                _, yi0, _r = chief2(p_probe)
                if yi0: break
                p_probe *= 0.5
            if not yi0:
                vig.append([0.0, 0.0, 0.0, 0.0]); continue
            lo, hi = 0.0, 1.6*p_probe*Y/abs(yi0)
            for _ in range(34):
                m = 0.5*(lo+hi); _, yi, _r = chief2(m)
                if yi is None: hi = m
                elif abs(yi) < Y: lo = m
                else: hi = m
            p = 0.5*(lo+hi)
        yep, _yi, _real_chief = chief2(p)
        if yep is None:
            vig.append([0.0, 0.0, 0.0, 0.0]); continue
        # 主光线自己被挡时**绝不能**退回全 0 —— 全 0 = 满光瞳，Zemax 会据此把上游各面的
        # 自动口径撑到造不出来（实测 JP2025052870A：0.06x 结构视场1 退成全 0 后，
        # 面9~12 的自动口径被撑到 31.3~31.7，而真实光束只有 27.7）。
        # 正确做法：扫全入瞳取最长连续可行区间，用它的中点当搜索种子；VDY 仍相对主光线。
        yseed = yep
        if shoot(p, yep)[1] is not None:
            NS = 400
            okv = [shoot(p, yep - rEP + 2*rEP*i/NS)[1] is None for i in range(NS+1)]
            best = (0, -1, -1); i = 0
            while i <= NS:
                if okv[i]:
                    j = i
                    while j + 1 <= NS and okv[j+1]: j += 1
                    if j - i > best[0]: best = (j - i, i, j)
                    i = j + 1
                else: i += 1
            if best[2] < 0:
                vig.append([0.0, 0.0, 0.0, 0.0]); continue
            yseed = yep - rEP + 2*rEP*(0.5*(best[1]+best[2]))/NS
        # ---- 光阑参考的光瞳坐标（孔径模型启用时）----
        # Zemax RAIM Real：Py=±1 是「实光线打在光阑 主光线高度 ± 光阑近轴半径」，
        # 不是入瞳上的 yep±rEP。两者差的是光瞳像差（轴上球差 + 离轴彗差/畸变），
        # 单一缩放补不齐：回归实测 RF100 离轴要 OpticStudio 补 4~7 步，RF35 轴上凭空多出渐晕。
        # 所以逐视场二分出 Py=±1 对应的入瞳高度当搜索边界，VDY/VCY/VCX 再折回光阑坐标。
        stopref = bool(rEP_fixed)
        if stopref:
            def hstop(ye):
                hs_, _b, _yi = shoot(p, ye, apert=False)
                return hs_[ks] if (hs_ and len(hs_) > ks) else None
            s0 = hstop(yep)
            h1 = hstop(yep + 0.1*rEP)
            if s0 is None or h1 is None or h1 == s0:
                stopref = False
            else:
                gs = 1.0 if h1 > s0 else -1.0
                def ye_at(target, lo, hi):
                    # 从 lo 往 hi 分步扫（远端实光线常追不出来 → 直接两端夹会失败），找到第一个跨过 target 的小区间再二分
                    N = 64
                    xa, fa = lo, hstop(lo)
                    if fa is None: return None
                    xp = fp = None
                    for i in range(1, N+1):
                        xb = lo + (hi - lo)*i/N; fb = hstop(xb)
                        if fb is None:
                            # 追迹器对超出口径 ~5% 的光线直接判追失（窗口按 semi 开），
                            # 被前片挡掉的那一侧常常追不到光阑。那一侧反正已被挡，边界只作搜索上限用 ——
                            # 拿最后两点线性外推。★ 只许**朝着 target 的方向**外推（回归审查：
                            # 往上扫却外推出下方的解，RF16 有 24 处，轴上视场都被写出假渐晕）。
                            if xp is None or fa == fp: return None
                            if (target - fa)*(fa - fp) <= 0: return None
                            return xa + (target - fa)*(xa - xp)/(fa - fp)
                        if (fa - target)*(fb - target) <= 0:
                            for _ in range(50):
                                mid = 0.5*(xa + xb); fm = hstop(mid)
                                if fm is None: return None
                                if (fa - target)*(fm - target) <= 0: xb, fb = mid, fm
                                else: xa, fa = mid, fm
                            return 0.5*(xa + xb)
                        xp, fp = xa, fa
                        xa, fa = xb, fb
                    return None
                def ye_find(target):
                    # 按 target 在主光线哪一侧决定扫描方向（光阑高度随入瞳高度单调，朝向 gs）
                    up = gs*(target - s0) >= 0
                    r_ = ye_at(target, yep, yep + (1.6 if up else -1.6)*rEP)
                    return r_ if r_ is not None else ye_at(target, yep, yep + (-1.6 if up else 1.6)*rEP)
                def py_of(ye):
                    h = hstop(ye)
                    return None if h is None else gs*(h - s0)/r_sp
                ye_top = ye_at(s0 + gs*r_sp, yep, yep + 1.6*rEP)
                ye_bot = ye_at(s0 - gs*r_sp, yep, yep - 1.6*rEP)
                if ye_top is None or ye_bot is None: stopref = False
        top_lim = ye_top if stopref else yep + rEP
        bot_lim = ye_bot if stopref else yep - rEP
        def edge_seed(sign):
            ok_, bad_ = yseed, (top_lim if sign > 0 else bot_lim)
            if shoot(p, bad_)[1] is None: return bad_
            for _ in range(44):
                m = 0.5*(ok_+bad_)
                if shoot(p, m)[1] is None: ok_ = m
                else: bad_ = m
            return ok_
        hi_ok, lo_ok = edge_seed(+1.0), edge_seed(-1.0)
        # 统计卡光面（粗扫，只为找渐晕定义面）
        blk_cnt = {}
        for i in range(121):
            ye = yep - rEP + 2*rEP*i/120
            b = shoot(p, ye)[1]
            if b is not None:
                blk_cnt[b] = blk_cnt.get(b, 0) + 1
                blk_all[b] = blk_all.get(b, 0) + 1
        # --- 两端各内缩 margin×rEP，保证 Py=±1 真的过得去 ---
        m_ep = margin*rEP
        hi_ok = min(hi_ok - m_ep, top_lim) if hi_ok < top_lim - 1e-9 else top_lim
        lo_ok = max(lo_ok + m_ep, bot_lim) if lo_ok > bot_lim + 1e-9 else bot_lim
        if stopref:
            py_hi = 1.0 if hi_ok >= top_lim - 1e-9 else (py_of(hi_ok) if py_of(hi_ok) is not None else 1.0)
            py_lo = -1.0 if lo_ok <= bot_lim + 1e-9 else (py_of(lo_ok) if py_of(lo_ok) is not None else -1.0)
            vdy = 0.5*(py_hi + py_lo)
            vcy = 1 - 0.5*(py_hi - py_lo)
        else:
            vdy = ((hi_ok+lo_ok)/2 - yep)/rEP
            vcy = 1 - (hi_ok-lo_ok)/(2*rEP)
        vdy = 0.0 if abs(vdy) < 0.002 else vdy
        vcy = 0.0 if vcy < 0.002 else vcy
        if stopref:
            yc = ye_find(s0 + gs*vdy*r_sp)
            if yc is None: yc = 0.5*(hi_ok + lo_ok)
        else:
            yc = yep + vdy*rEP
        # VCX：光瞳中心高度上横向扫描（Zemax 的椭圆光瞳，X 半轴取在椭圆中心处）
        vcx = 0.0
        px_lim = rEP
        if stopref:
            yc_s = hstop(yc) or 0.0
            def xstop(px):
                rr_ = shoot3n(p, px, yc)[0]
                if not rr_ or len(rr_) <= ks: return None
                return math.sqrt(max(rr_[ks]**2 - yc_s**2, 0.0))
            def px_at(target):
                # 同 ye_at：从 0 往外分步扫，追不出来就拿最后两点线性外推
                xa_, fa_, xp_, fp_ = 0.0, 0.0, None, None
                for i in range(1, 65):
                    xb_ = 1.6*rEP*i/64; fb_ = xstop(xb_)
                    if fb_ is None:
                        if xp_ is None or fa_ <= fp_ or target <= fa_: return None
                        return xa_ + (target - fa_)*(xa_ - xp_)/(fa_ - fp_)
                    if fb_ >= target:
                        for _ in range(50):
                            mid = 0.5*(xa_ + xb_); fm = xstop(mid)
                            if fm is None or fm >= target: xb_ = mid
                            else: xa_ = mid
                        return 0.5*(xa_ + xb_)
                    xp_, fp_ = xa_, fa_
                    xa_, fa_ = xb_, fb_
                return None
            _pl = px_at(r_sp)
            if _pl: px_lim = _pl
        if shoot3(p, 0.0, yc)[1] is None:
            last, hit = 0.0, False
            for i in range(1, 201):
                px = i/200.0*px_lim
                if shoot3(p, px, yc)[1] is None: last = px
                else: hit = True; break
            if hit:            # X 向没挡住任何东西时 VCX 就是 0，绝不能再减余量
                lo2, hi2 = last, min(px_lim, last + px_lim/200.0)
                for _ in range(40):
                    m = 0.5*(lo2+hi2)
                    if shoot3(p, m, yc)[1] is None: lo2 = m
                    else: hi2 = m
                if stopref:
                    _pe = max(0.5*(lo2+hi2) - m_ep, 0.0)
                    xe = xstop(_pe)
                    if xe is None: xe = r_sp*_pe/px_lim          # 追不出来就按线性映射
                    vcx = max(0.0, 1 - xe/r_sp)
                else:
                    vcx = max(0.0, 1 - (0.5*(lo2+hi2) - m_ep)/rEP)
        vcx = 0.0 if vcx < 0.002 else vcx
        # --- 椭圆边界体检 + 收缩：真实通光区是几个圆的交集，椭圆会在斜方位上探出去 ---
        ax, ay = (1-vcx)*rEP, (1-vcy)*rEP
        # 真实通光区是几个圆的交集，椭圆必然在斜方位上探出去一点 —— 那是 Zemax 椭圆
        # 渐晕模型的固有近似，OpticStudio 自己也会把那几条切掉，**不要拿它去收缩光瞳**
        # （实测缩到边界零遮挡要多丢 12% 的光，VCX 从 0.02 被推到 0.16，MTF/相对照度全废）。
        # 必须保证的只有两组：Layout 画的 Py=±1（Px=0）和光瞳中心高度上的 Px=±1，
        # 这两组由上面的精确二分 + margin 保证。--fit-ellipse 才收缩。
        if fit_ellipse:
            for _ in range(60):
                if ellipse_blocked(p, yc, ax, ay) == 0: break
                if ax > 0.35*rEP: ax *= 0.985
                else: ax *= 0.985; ay *= 0.985
        vcx = round(min(0.9, max(0.0, 1 - ax/rEP)) + 5e-5, 4)
        vcy = round(min(0.9, max(0.0, 1 - ay/rEP)) + 5e-5, 4)
        vdy = round(vdy, 4)
        if vcx < 0.002: vcx = 0.0
        if vcy < 0.002: vcy = 0.0
        # 有限共轭时物高取正 → 像高为负，解的是 −Y 那个视场点；Zemax 的视场是
        # 「实际像高 +Y」，两者关于光轴镜像，所以 VDY 要翻号（VCY/VCX 对称不变）。
        if obj is not None and vdy: vdy = round(-vdy, 4)
        vig.append([0.0, vdy, vcx, vcy])
        # 最终验证：Py=±1（Layout 画的就是这两条）与椭圆边界 72 点
        ax, ay = (1-vcx)*rEP, (1-vcy)*rEP
        yc = yep + (-vdy if obj is not None else vdy)*rEP
        if stopref:
            _vd = -vdy if obj is not None else vdy          # 回到求解时的朝向
            _yt = ye_find(s0 + gs*(_vd + (1-vcy))*r_sp)
            _yb = ye_find(s0 + gs*(_vd - (1-vcy))*r_sp)
            if _yt is not None and _yb is not None:
                yc, ay = 0.5*(_yt + _yb), 0.5*abs(_yt - _yb)
                _ycs = hstop(yc) or 0.0
                _lo, _hi = 0.0, 1.6*rEP
                _axp = None
                try:
                    _axp = px_at((1-vcx)*r_sp) if (1-vcx) > 0 else 0.0
                except NameError:
                    _axp = None
                ax = _axp if _axp is not None else (1-vcx)*px_lim
        pyok = (all(shoot(p, yc + s*ay)[1] is None for s in (-1.0, 1.0))
                and all(shoot3(p, s*ax, yc)[1] is None for s in (-1.0, 1.0)))
        eb = ellipse_blocked(p, yc, ax, ay, 72)
        for i in range(24):
            th = 2*math.pi*i/24
            for g in (1.0, 0.7, 0.4):
                rr, b3, _iy, _ix = shoot3(p, g*ax*math.cos(th), yc + g*ay*math.sin(th))
                if b3 is not None: continue
                for k2, v2 in enumerate(rr): maxr[k2] = max(maxr[k2], v2)
        beams.append((Y, p, vdy, vcy, vcx, pyok, eb))
        if verbose:
            who = sorted(blk_cnt.items(), key=lambda x: -x[1])[:2]
            ang_deg = math.degrees(p) if obj is None else math.degrees(math.atan2(Y, 1e9))
            head = ('  %6.2f   %6.2f°' % (Y, math.degrees(p))) if obj is None \
                   else ('  %6.2f   物高%6.1f' % (Y, p))
            print('%s  %+7.4f  %7.4f  %7.4f   %5.1f%%   %s%s   %s'
                  % (head, vdy, vcy, vcx, 100*(1-vcy)*(1-vcx),
                     'OK' if pyok else '★Py/Px±1被挡', '' if eb == 0 else ' 斜%d/72' % eb,
                     ', '.join('面%d(%d)' % (k+1, v) for k, v in who) or '-'))
    return vig, maxr, blk_all, beams, axr


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('spec')
    ap.add_argument('--emb', type=int, default=0); ap.add_argument('--state', default=None)
    ap.add_argument('--margin', type=float, default=0.010,
                    help='光瞳两端内缩比例（×入瞳半径），默认 0.010；保证 Layout 的 Py=±1 不被切')
    ap.add_argument('--first-only', action='store_true', help='只解基准结构（旧行为）')
    ap.add_argument('--fit-ellipse', action='store_true',
                    help='把椭圆缩到斜方位也零遮挡（很费光，一般不要开）')
    ap.add_argument('--write', nargs='?', const='AUTO', default=None,
                    help='把 vignetting / vignetting_cfg / fix_semi_surfaces 写回 spec')
    a = ap.parse_args()
    spec = json.load(open(a.spec, encoding='utf-8')); emb = spec['embodiments'][a.emb]
    zx = spec['zmx']; state = a.state or emb['states'][0]
    cfgs = zx.get('configs') or []
    kb1 = (zx.get('focus') or {}).get('key_before')
    kb2 = (zx.get('focus2') or {}).get('key_before')
    if a.first_only or not cfgs:
        cfgs = [{'name': state, 'd0': 'INFINITY'}]

    stop_semi, aprows = aperture_cfg(spec, emb, state, cfgs)
    fpat = zx.get('fno_patent') or {}
    print('== 孔径：Paraxial Working F/# 逐结构（∞ 结构光阑近轴半径 %.4f）==' % stop_semi)
    print('  %-12s %8s %7s %9s %9s %9s %9s  %s'
          % ('结构', 'EFL', '|β|', '光阑不变', '专利印刷', '采用WFNO', '光阑半径', '来源'))
    for r in aprows:
        print('  %-12s %8.3f %7.4f %9.3f %9s %9.3f %9.4f  %s'
              % (r['name'], r['efl'], r['beta'], r['wfno_fixed'],
                 ('%.2f' % fpat[r['name']]) if r['name'] in fpat else '-',
                 r['wfno'], r['stop_semi'], r['src']))
    for t in getattr(aperture_cfg, 'notes', []):
        print('  ⚠ ' + t)
    if any(r['src'] != '固定光阑' for r in aprows):
        print('  注：专利近距 F 数比「光阑全开不变」暗 —— 这支镜头近距会收光圈（或专利按别的口径定义），'
              '按专利值走；中间结构是插值的，回话要说')
    vig_cfg, maxr, blk_all, axall = [], None, {}, None
    for c, ar in zip(cfgs, aprows):
        dmap = cfg_dmap(spec, emb, c)
        d0 = c.get('d0')
        obj = None if (d0 is None or str(d0).upper().startswith('INF')) else float(d0)
        S, zimg, _ = build(spec, emb, state, dmap)
        v, mr, blk, _b, axr = solve_state(S, zimg, zx, obj, a.margin, True,
                                          '\n== %s ==  ' % c.get('name', '?'), a.fit_ellipse,
                                          rEP_fixed=ar['rEP'])
        vig_cfg.append(v)
        maxr = mr if maxr is None else [max(x, y) for x, y in zip(maxr, mr)]
        axall = axr if axall is None else [max(x, y) for x, y in zip(axall, axr)]
        for k, n in blk.items(): blk_all[k] = blk_all.get(k, 0) + n

    S, zimg, _ = build(spec, emb, state)
    print('\n  逐面：全结构 3D 渐晕光束最大半径 vs 有効径/2（负余量 = 该面在切光）')
    bad = []
    for k, sf in enumerate(S):
        if sf.semi is None: continue
        m = maxr[k]; mark = ''
        if sf.stop and aprows:
            mark = '  （光阑：大小由孔径定义给出，不参与挡光）'
        elif m - sf.semi > 0.01: mark = '  ← ★ 光束超出口径'; bad.append(k+1)
        elif sf.semi - m < 0.30: mark = '  ← 贴边(<0.30)'
        print('    面%-3d 光束 %6.3f   有効径/2 %6.3f   余量 %+6.3f%s'
              % (k+1, m, sf.semi, sf.semi-m, mark))
    if bad:
        print('    ★ 面 %s 的口径比实际光束还小 —— 会把离轴光瞳边缘误切掉，先回去核口径' % bad)
    # ---- 轴上满光瞳体检：这一条比渐晕严重得多，切了就是把 F 数改掉 ----
    print('\n  逐面：全结构**轴上满光瞳**需求半径 vs 有効径/2（这一圈绝不许切，切了 F 数就变了）')
    axbad = []
    for k, sf in enumerate(S):
        if sf.semi is None or not axall: continue
        if sf.stop and aprows: continue    # 光阑大小由孔径定义给出，不是挡光口径
        need = axall[k]; mg = sf.semi - need
        mark = ''
        if mg < -0.001:
            mark = '  ← ★★ 切到轴上光瞳，F 数会变大'; axbad.append((k+1, need, sf.semi))
        elif mg < 0.10: mark = '  ← 贴边'
        if mark:
            print('    面%-3d 轴上需求 %6.3f   有効径/2 %6.3f   余量 %+6.3f%s'
                  % (k+1, need, sf.semi, mg, mark))
    if axbad:
        worst = min(axbad, key=lambda t: t[2] / t[1])
        print('    ★★ 共 %d 个面在切轴上光瞳，最严重的是面%d（%.3f < %.3f，切掉 %.1f%% 半径，'
              '等效 F 数约 %.3f）—— 口径必须放大到需求值，或者接受相对孔径变小'
              % (len(axbad), worst[0], worst[2], worst[1],
                 100*(1-worst[2]/worst[1]), zx['fno']*worst[1]/worst[2]))
    else:
        print('    全部通过：没有任何面切到轴上满光瞳')
    rank = sorted(blk_all.items(), key=lambda x: -x[1])
    print('\n全结构卡光面排名:', ', '.join('面%d×%d' % (k+1, v) for k, v in rank[:6]))
    tt = sum(v for _, v in rank) or 1
    defs = sorted(k+1 for k, v in rank if v/tt > 0.05)
    blocks, cur = [], None
    for k, sf in enumerate(S):
        if sf.n and abs(sf.n - 1.0) > 1e-9:
            cur = [k, k+1] if cur is None else [cur[0], k+1]
        elif cur is not None:
            blocks.append(tuple(cur)); cur = None
    if cur is not None: blocks.append(tuple(cur))
    grown = set(defs)
    for d in defs:
        for a0, b0 in blocks:
            if a0 <= d-1 <= b0: grown |= {k+1 for k in range(a0, b0+1)}
    defs = sorted(grown)
    print('渐晕定义面（同一镜片各面一并固定）:', defs)
    if a.write:
        out = a.spec[:-5] + '.vig.json' if a.write == 'AUTO' else a.write
        spec['zmx']['vignetting'] = vig_cfg[0]
        spec['zmx']['vignetting_cfg'] = vig_cfg
        spec['zmx']['semi_3d'] = [round(v, 3) for v in maxr]
        if axall: spec['zmx']['axial_3d'] = [round(v, 3) for v in axall]
        spec['zmx']['aperture'] = {
            'type': 'paraxial_working_fno',
            'stop_semi_paraxial': round(stop_semi, 6),
            'wfno_cfg': [round(r['wfno'], 4) for r in aprows],
            'wfno_fixed_stop': [round(r['wfno_fixed'], 4) for r in aprows],
            'stop_semi_cfg': [round(r['stop_semi'], 5) for r in aprows],
            'beta_cfg': [round(r['beta'], 6) for r in aprows],
            'source': [r['src'] for r in aprows],
            'epd_cfg': [round(2*r['rEP'], 4) for r in aprows],
            'configs': [r['name'] for r in aprows]}
        # 渐晕定义面单独记一份；fix_semi_surfaces 只并入、不覆盖 —— 以前直接覆盖成 defs，
        # 专利/断面图口径（aptrace 或手工全量固定）会退回自动，Zemax 按当前结构算口径，
        # 其它变焦位置轴上光线就显得被切（JP2023-039817A 70-200 GM II 踩过）。
        spec['zmx']['vig_def_surfaces'] = defs
        old = spec['zmx'].get('fix_semi_surfaces') or []
        spec['zmx']['fix_semi_surfaces'] = sorted(set(old) | set(defs), key=lambda v: (not isinstance(v, int), v))
        json.dump(spec, open(out, 'w', encoding='utf-8'), ensure_ascii=False, indent=1)
        print('已写出 %s（vignetting_cfg: %d 个结构 × %d 视场）' % (out, len(vig_cfg), len(vig_cfg[0])))


if __name__ == '__main__':
    main()
