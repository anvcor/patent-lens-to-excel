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


# --- 多波长折射率 -------------------------------------------------------------
# Zemax 的孔径（Paraxial Working F/#）、实光线瞄准、Set Vignetting 都按**主波长**（默认 e 线 0.5461）算，
# 轴上光束则是全部系统波长的包络。只拿 nd 追迹，光阑近轴半径就差 0.2%（TS-E24：9.7733 vs 9.7937）。
# 折射率来源（按优先级）：
#   ① 目录牌号：读 AGF 色散公式（与 lensmath 同一套 _n），再平移到 spec 的 nd；
#   ② Offset 玻璃解：基准牌号的色散曲线，按 (nd−1)/νd 与基准之比缩放 ——
#      对照 OpticStudio 2024R2 INDX（TS-E24 八个 Offset 面）：e 线 ≤1.5e-6、g 线 ≤5e-5；
#   ③ 模型玻璃 / 查不到目录：nd、νd、dPgF 定三项 Cauchy（n = A + B/λ² + C/λ⁴）。
#      对三家全目录：e 线 ≤2e-5，F/C/g 最差 5e-4（S-NPH3 这类超重火石）。
WAVES_DEFAULT = [(0.4861, 12.0), (0.5461, 30.0), (0.6563, 3.0), (0.5876, 22.0), (0.4358, 3.0)]
PWAV_DEFAULT = 2
LAM_D = 0.5875618
_LG, _LF, _LC = 0.4358343, 0.4861327, 0.6562725


def waves_of(zx):
    """(主波长 µm, [全部波长 µm])。与 make_zmx.waves_of / PWAV_DEFAULT 同一套约定。"""
    if zx.get('waves'):
        wv = [float(w[0]) for w in zx['waves']]
    elif zx.get('waves_um'):
        wv = [float(w) for w in zx['waves_um']]
    else:
        wv = [w for w, _ in WAVES_DEFAULT]
    k = int(zx.get('primary_wave', PWAV_DEFAULT))
    return wv[min(max(k, 1), len(wv)) - 1], wv


_CAT = {}
def _catalog(vendor, gcat=None):
    """{牌号: (公式号, CD 系数)}；按 PATENT_GLASS_DIR → 用户 Zemax Glasscat → skill 自带目录 找。"""
    stem = (gcat or {}).get(vendor) or vendor
    if (vendor, stem) in _CAT: return _CAT[(vendor, stem)]
    import os
    from lensmath import _decode
    here = os.path.dirname(os.path.abspath(__file__))
    dirs = [os.environ.get('PATENT_GLASS_DIR'),
            os.path.join(os.path.expanduser('~'), 'Documents', 'Zemax', 'Glasscat'),
            os.path.join(here, '..', 'assets', 'glass')]
    out = {}
    for d in dirs:
        if not d: continue
        for nm in (stem, vendor):
            p = os.path.join(d, nm + '.AGF')
            if not os.path.exists(p): continue
            cur = None
            for line in _decode(open(p, 'rb').read()).splitlines():
                if line.startswith('NM '):
                    q = line.split(); cur = q[1]
                    try: out[cur] = [int(float(q[2])), None]
                    except (IndexError, ValueError): cur = None
                elif line.startswith('CD ') and cur in out:
                    out[cur][1] = [float(x) for x in line.split()[1:]]
            break
        if out: break
    _CAT[(vendor, stem)] = out
    return out


def _n_cat(glass, lam, gcat=None):
    """'OHARA S-BAL42' → 目录色散公式在 lam 处的折射率（查不到返回 None）。"""
    from lensmath import _n
    q = (glass or '').split()
    if len(q) < 2: return None
    f, cd = _catalog(q[0], gcat).get(q[1], (None, None))
    return _n(cd, f, lam) if cd else None


def _n_cauchy(nd, vd, dpgf, lam):
    dFC = (nd - 1.0) / vd
    P = 0.6438 - 0.001682 * vd + (dpgf or 0.0)
    a1, b1 = _LF**-2 - _LC**-2, _LF**-4 - _LC**-4
    a2, b2 = _LG**-2 - _LF**-2, _LG**-4 - _LF**-4
    det = a1 * b2 - a2 * b1
    B = (dFC * b2 - P * dFC * b1) / det
    C = (a1 * P * dFC - a2 * dFC) / det
    return nd + B * (lam**-2 - LAM_D**-2) + C * (lam**-4 - LAM_D**-4)


def n_at(row, lam, gcat=None):
    """spec 一行的玻璃在波长 lam(µm) 处的折射率；lam=None 或 d 线直接返回 nd。"""
    nd = row.get('nd') or 1.0
    if lam is None or abs(lam - LAM_D) < 1e-5 or nd == 1.0: return nd
    go = row.get('glass_offset')
    if go:
        nb, nbd = _n_cat(go.get('base'), lam, gcat), _n_cat(go.get('base'), LAM_D, gcat)
        vd, vb = row.get('vd'), go.get('base_vd')
        if nb is not None and nbd is not None and vd and vb:
            return nd + (nb - nbd) * ((nd - 1.0) / vd) / ((nbd - 1.0) / vb)
    elif row.get('glass'):
        nc, ncd = _n_cat(row['glass'], lam, gcat), _n_cat(row['glass'], LAM_D, gcat)
        if nc is not None and ncd is not None and abs(ncd - nd) < 2e-4:
            return nc + (nd - ncd)
    if row.get('vd'):
        return _n_cauchy(nd, float(row['vd']), row.get('dpgf'), lam)
    return nd
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

def build(spec, emb, state, dmap=None, lam=None):
    """lam=None：按 nd（d 线）建系统（近轴孔径模型、口径体检都用这个）；
    给 lam(µm)：各面折射率换成该波长（n_at），DOE 偏折按 λ/λ0 缩放（Binary 2 相位固定，偏折 ∝ λ）。"""
    asph = {str(a['surface']).replace('面', ''): a for a in emb.get('aspheric', [])}
    gcat = (spec.get('zmx') or {}).get('gcat')
    S, z = [], 0.0
    rows = [s for s in emb['surfaces'] if s['i'] != 'IMG']
    img = [s for s in emb['surfaces'] if s['i'] == 'IMG']
    for s in rows:
        D = s['D']
        if isinstance(D, str): D = (dmap or {}).get(D, emb['variable'][D][state])
        a = asph.get(str(s['i']))
        A = acoef(a) if a else []
        phi = (s.get('extra') or {}).get('有効径 φi')
        doe = (s.get('doe') or {}).get('C')
        if doe and lam is not None:
            k_l = lam / (float(s['doe'].get('wl_nm', 587.56)) * 1e-3)
            doe = [float(c) * k_l for c in doe]
        S.append(Surf(0.0 if s['R'] in (None, 0) else 1.0/float(s['R']),
                      (a or {}).get('k', (a or {}).get('K', 0.0)) or 0.0, A, z,
                      n_at(s, lam, gcat) if s.get('nd') else 1.0, (phi/2.0) if phi else None,
                      bool(s.get('stop') or s['i'] == 'STO'), doe))
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


def trace(S, zimg, y0, z0, ang, apert=True, wide=False):
    """从 (y0,z0) 以角 ang(rad) 出发追子午实光线。返回 (各面高度, 被挡的面序号或None, 像高)。

    求交窗口默认按口径开（1.05×semi+0.2，面型多项式只在净口径附近才有意义），超出就当追失。
    wide=True 把窗口放到 1.5×semi+1：量「本来该走到多高」（轴上底线）时，口径恰恰可能已经切进光束 ——
    RF35 MFD 的 g 线轴上光线在面12 要到 16.994，窗口 16.985 直接追失，底线就量不出来。"""
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
        _yl = ((1.5 * s.semi + 1.0) if wide else (1.05 * s.semi + 0.2)) if s.semi is not None else 1e4
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

def trace3(S, zimg, P0, d, apert=True, xy=None):
    """3D 斜光线。旋转对称系统：面型只依赖 r=hypot(x,y)，法线 = (-sl*x/r, -sl*y/r, 1)。
    返回 (各面处的 r, 被挡面序号或 None, 像面 y, 像面 x)。给 xy=[] 时顺便记下各面交点 (x, y)。"""
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
        if xy is not None: xy.append((x, y))
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


def solve_state(S, zimg, zx, obj=None, margin=0.010, verbose=True, tag='', fit_ellipse=False, rEP_fixed=None,
                wfno=None, S_waves=()):
    """解一个结构（一组间隔 + 一个物距）的渐晕。

    obj=None 表示物在无限远，视场参数 p = 入射半角；否则 obj 是物距（mm，面1顶点起算），
    视场参数 p = 物高。两种情形都用「入瞳高度 ye」参数化光瞳（Py = (ye-yep)/rEP）。

    S 应按 .zmx 的**主波长**建（build(..., lam=主波长)）；wfno = 该结构的 Paraxial Working F/#
    （给了就在主波长上重算入瞳半径，与 OpticStudio 的 EPD 逐位一致）；
    S_waves = 其余系统波长的同一几何，只用于轴上满光瞳包络（axial_3d）。

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
        # 光阑面本身不参与挡光：它的大小已由孔径定义给出（Zemax 开 Real 瞄准时 Py=±1 打在
        # 「轴上实边缘光线在光阑上的高度」上，见下面 r_sp）。光阑再拿专利有効径去挡，
        # 近距结构会凭空在轴上写出 VCY≈0.017 的假渐晕（RF100：15.51 > 15.405）。
        for _s in S:
            if _s.stop: _s.semi = None
    ca, cb = -yB[ks], yA[ks]
    y1 = ca*yA[0] + cb*yB[0]; u1 = cb
    zEP = -y1/u1 if u1 else 0.0
    z_obj = None if obj is None else -abs(obj)
    if rEP_fixed and wfno:
        # 主波长上按 Paraxial Working F/# 重算入瞳半径（OpticStudio 的 EPD 就是这么来的）：
        # 近轴边缘光线像方斜率 = 1/(2·WFNO)。∞：平行光入射；有限共轭：轴上物点出发，
        # 入瞳半径 = 这条光线在近轴入瞳面上的高度。TS-E24 四个结构 EPD 与 EPDI 对到 1e-5。
        if obj is None:
            rEP = 0.5/(float(wfno)*abs(uA))
        else:
            _yo, _uo = par(abs(obj), 1.0)
            rEP = 0.5/(float(wfno)*abs(_uo)) * (zEP + abs(obj))

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

    r_par = rEP * abs(yA[ks])                  # 光阑近轴半径
    r_sp = r_par
    if rEP_fixed:
        # ★ OpticStudio「Ray Aiming = Real」的光阑半径（2024R2 ZOS-API 实测，TS-E24 / RF5.2 鱼眼逐位对上）：
        #   R = 主波长轴上实光线**瞄近轴入瞳边缘**（∞：入射高 = EPD/2 的平行光；有限共轭：
        #       轴上物点 → 近轴入瞳面上高 EPD/2 那一点的直线）在光阑面上的实际高度；
        #   之后每个视场、每个波长的 (Px,Py) 都线性落在光阑上：光阑坐标 = 主光线 + (Px,Py)·R。
        # 入瞳有球差时 R 比近轴半径大（TS-E24 10.061 vs 9.794，+2.7%；RF5.2 4.020 vs 3.977）——
        # 专利印的光阑有効径（TS-E24 20.09 → 10.045）对应的也是这个实半径。
        # 旧版反过来把入瞳缩到「实光线落在近轴半径上」（×0.975），轴上底线和光瞳都偏小一圈：
        # TS-E24 的 axial_3d 说全过，OpticStudio 里面14~22 切到轴上光束。
        hs0 = trace(S, zimg, *start(0.0, rEP), apert=False, wide=True)[0]
        if hs0 and len(hs0) > ks and hs0[ks]:
            r_sp = abs(hs0[ks])

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
    # 光束 = OpticStudio 的轴上光束：每个系统波长各自瞄到光阑实半径 r_sp（Real 瞄准对所有波长同一个 R），
    # 取全部波长、全部光瞳半径的逐面包络。轴上旋转对称，子午一条线扫到底即可；
    # 不只取边缘光线 —— 强光瞳像差时某些面的最高点在光瞳中间（焦散）。
    axr = [0.0]*len(S)
    for SS in [S] + list(S_waves or ()):
        def _hs(ye, _SS=SS):
            y0, z0, ang = start(0.0, ye)
            return trace(_SS, zimg, y0, z0, ang, apert=False, wide=True)[0]
        ye_e = rEP
        if SS is not S and rEP_fixed:
            # 从 rEP 出发做定点迭代（斜率取 rEP/r_sp，色差只差百分之几，几步就收敛）。
            # 不用大区间二分：1.15×rEP 那头在前片上超出追迹窗口直接追失，会静默退回 rEP。
            ok_ = None
            for _ in range(60):
                h = _hs(ye_e)
                if len(h) <= ks:                   # 追失：退回上一个好点和这里的中点，不许静默退回 rEP
                    if ok_ is None: ye_e = rEP; break
                    ye_e = 0.5*(ok_ + ye_e); continue
                ok_ = ye_e
                f = abs(h[ks]) - r_sp
                if abs(f) < 1e-9: break
                ye_e -= f * rEP / r_sp
        for _i in range(40, 0, -1):
            for _k, _v in enumerate(_hs(ye_e*_i/40.0)): axr[_k] = max(axr[_k], abs(_v))
    if verbose:
        print('%s入瞳半径 %.4f  入瞳位置 面1前 %.3f%s%s'
              % (tag, rEP, -zEP, '' if obj is None else '   物距 %.2f' % obj,
                 '' if abs(r_sp - r_par) < 1e-6 else
                 '   光阑实半径 %.4f（近轴 %.4f，%+.2f%%；Real 瞄准 Py=±1 落在这里）'
                 % (r_sp, r_par, 100*(r_sp/r_par - 1))))
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
        # Zemax RAIM Real：Py=±1 是「实光线打在光阑 主光线高度 ± 光阑实半径 r_sp」（r_sp 见上），
        # 不是入瞳上的 yep±rEP。两者差的是光瞳像差（轴上球差 + 离轴彗差/畸变），
        # 单一缩放补不齐：回归实测 RF100 离轴要 OpticStudio 补 4~7 步。
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
                def ye_dir(target, sign):
                    # 扫描范围逐级放宽：大视场光瞳畸变极大，Py=±1 在入瞳上可以离主光线 2~3 个 rEP
                    # （TS-E24 54°：Py=+0.6 已在 yep+1.7rEP）。只扫 1.6rEP、扫完没跨过就返回 None，
                    # stopref 会静默退回入瞳坐标 —— 把 Py≈−0.33 当成 −1，VCY 0.21 vs OpticStudio 0.73。
                    for span in (1.6, 4.0, 8.0):
                        r_ = ye_at(target, yep, yep + sign*span*rEP)
                        if r_ is not None: return r_
                    return None
                def ye_find(target):
                    # 按 target 在主光线哪一侧决定扫描方向（光阑高度随入瞳高度单调，朝向 gs）
                    up = 1.0 if gs*(target - s0) >= 0 else -1.0
                    r_ = ye_dir(target, up)
                    return r_ if r_ is not None else ye_dir(target, -up)
                def py_of(ye):
                    h = hstop(ye)
                    return None if h is None else gs*(h - s0)/r_sp
                ye_top = ye_dir(s0 + gs*r_sp, +1.0)
                ye_bot = ye_dir(s0 - gs*r_sp, -1.0)
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
            ye = bot_lim + (top_lim - bot_lim)*i/120
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
        aim_x = None
        if stopref:
            # ★ Zemax 的 Px 是**二维瞄准**：Px=t 的实光线精确落在光阑点 (t·R, 光瞳中心 y) 上。
            # 旧写法固定入瞳 y、只改入瞳 x，再用 sqrt(r²−y²) 估光阑 x —— 斜光线在光阑上的 y 会漂，
            # VCX 系统性偏小，回归里 6 格 Px±1 在 OpticStudio 被挡（TS-E24 MFD 视场4、RF50 MFD 视场1 等）。
            yc_t = s0 + gs*vdy*r_sp                 # 光瞳中心在光阑上的 y（vdy 此时仍是求解朝向）
            def _stop_xy(px_, ye_):
                xy_ = []
                trace3(S, zimg, *start3(p, px_, ye_), apert=False, xy=xy_)
                return xy_[ks] if len(xy_) > ks else None
            def aim_x(t, guess):
                """入瞳 (px, ye)，使实光线落在光阑 (t·R, yc_t)；解不出返回 None。"""
                px_, ye_ = guess
                h_ = 1e-4*rEP
                # 门槛按光阑半径取相对值：求交本身有 ~1e-6 mm 的数值噪声（RF800 实测），
                # 绝对 1e-8 永远收不住，瞄准「失败」会被当成挡光（VCX 假到 0.55）。1e-5·R 对光瞳坐标可忽略。
                tol_a = 1e-7*r_sp
                best = None
                for _ in range(30):
                    q0 = _stop_xy(px_, ye_)
                    if q0 is None: return None
                    fx, fy = q0[0] - t*r_sp, q0[1] - yc_t
                    e_ = max(abs(fx), abs(fy))
                    if best is None or e_ < best[0]: best = (e_, px_, ye_)
                    if e_ < tol_a: return px_, ye_
                    qa, qb = _stop_xy(px_ + h_, ye_), _stop_xy(px_, ye_ + h_)
                    if qa is None or qb is None: break
                    a11, a21 = (qa[0] - q0[0])/h_, (qa[1] - q0[1])/h_
                    a12, a22 = (qb[0] - q0[0])/h_, (qb[1] - q0[1])/h_
                    det = a11*a22 - a12*a21
                    if det == 0: break
                    px_ -= (fx*a22 - fy*a12)/det
                    ye_ -= (a11*fy - a21*fx)/det
                return (best[1], best[2]) if best[0] < 100*tol_a else None
            g0 = aim_x(0.0, (0.0, yc))
            if g0 is None: aim_x = None
        if aim_x is not None:
            # 沿光阑 x 从 0 扫到 R（续接上一个解当初值），第一个被挡 / 瞄不出来的就是边
            last_t, last_g, hit_t = 0.0, g0, None
            for i in range(1, 51):
                t = i/50.0; g_ = aim_x(t, last_g)
                if g_ is None or shoot3(p, g_[0], g_[1])[1] is not None: hit_t = t; break
                last_t, last_g = t, g_
            if hit_t is not None:  # X 向没挡住任何东西时 VCX 就是 0，绝不能再减余量
                lo2, hi2 = last_t, hit_t
                for _ in range(30):
                    m = 0.5*(lo2 + hi2); g_ = aim_x(m, last_g)
                    if g_ is not None and shoot3(p, g_[0], g_[1])[1] is None: lo2, last_g = m, g_
                    else: hi2 = m
                vcx = max(0.0, 1 - (lo2 - margin))
        elif shoot3(p, 0.0, yc)[1] is None:
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
        pxg = None
        if stopref and aim_x is not None:
            yc_t = s0 + gs*_vd*r_sp                 # 按取整后的 VDY 重瞄 Px=±1
            pxg = aim_x(1 - vcx, (0.0, yc))
            if pxg is not None: ax = abs(pxg[0])
        pyok = (all(shoot(p, yc + s*ay)[1] is None for s in (-1.0, 1.0))
                and (all(shoot3(p, s*pxg[0], pxg[1])[1] is None for s in (-1.0, 1.0)) if pxg else
                     all(shoot3(p, s*ax, yc)[1] is None for s in (-1.0, 1.0))))
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
    # 追迹用 .zmx 的主波长（OpticStudio 的 EPD / Real 瞄准 / Set Vignetting 都按它算），
    # 其余波长只进轴上满光瞳包络。近轴孔径模型（aperture_cfg 的 WFNO）仍按 nd 定 —— 它只给 F 数。
    lam_p, lams = waves_of(zx)
    print('  追迹波长：主 %.4f µm；轴上包络 %s' % (lam_p, ', '.join('%.4f' % l for l in lams)))
    vig_cfg, maxr, blk_all, axall = [], None, {}, None
    for c, ar in zip(cfgs, aprows):
        dmap = cfg_dmap(spec, emb, c)
        d0 = c.get('d0')
        obj = None if (d0 is None or str(d0).upper().startswith('INF')) else float(d0)
        S, zimg, _ = build(spec, emb, state, dmap, lam=lam_p)
        Sw = [build(spec, emb, state, dmap, lam=l)[0] for l in lams if abs(l - lam_p) > 1e-9]
        v, mr, blk, _b, axr = solve_state(S, zimg, zx, obj, a.margin, True,
                                          '\n== %s ==  ' % c.get('name', '?'), a.fit_ellipse,
                                          rEP_fixed=ar['rEP'], wfno=ar['wfno'], S_waves=Sw)
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
