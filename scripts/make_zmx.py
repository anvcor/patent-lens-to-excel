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


def fields_of(zx, emb):
    """zmx.fields_y 直接给就用；否则按最大像高取 6 等分点。
    最大像高优先级：zmx.max_y > 各种数据里的 Y > 像面有効径/2。"""
    if zx.get('fields_y'):
        return [float(y) for y in zx['fields_y']]
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

def build(spec, emb, catalog, asph_mode='auto'):
    zx = spec['zmx']; fc = zx['focus']
    kb, ka = fc['key_before'], fc.get('key_after')
    # 链式三段浮动（d_a + d_b + d_c = const，如本篇适马 105 微距：G2 与絞り各自移动）：
    # 第三段 key_last 的 DISZ = 守恒和 − 前两段，位置解随后会覆盖它。
    kl = fc.get('key_last')
    # 第二个对焦群（双浮动对焦，如索尼 135GM）。没有就是 None，行为与以前完全一致。
    fc2 = zx.get('focus2')
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
    atype, refit, rf_notes, xa_notes = {}, {}, [], []
    for s in surfs:
        key = str(s['i']); A = asp.get(key)
        if not A: continue
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
    a('FNUM %s 0' % num(zx['fno']))
    a('ENVD 20 1 0'); a('GFAC 0 0')
    if catalog and vendors: a('GCAT ' + ' '.join(vendors) + ' ')
    # RAIM 第 2 位 = Ray Aiming：0=Off / 1=Paraxial / 2=Real（官方样例实证，见 SKILL）。
    # 大孔径 / 大视场 / 光阑前有强负透镜时必须开 Real，否则光瞳定位是错的。
    ra = {'off': 0, 'none': 0, 'paraxial': 1, 'real': 2}.get(
        str(zx.get('ray_aiming', 'real')).lower(), 2)
    a('RAIM 0 %d 1 1 0 0 0 0 0 1' % ra); a('PUSH 0 0 0 0 0 0'); a('SDMA 0 1 0')
    a('OMMA 1 1')
    fy = fields_of(zx, emb); wv = waves_of(zx)
    a('FTYP 3 0 %d %d 0 0 0' % (len(fy), len(wv)))
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
    poslen = _poslen(fc)
    poslen2 = _poslen(fc2) if fc2 else None
    for k, s in enumerate(surfs, 1):
        i = _idx(s); key = str(s['i'])
        A = asp.get(key)
        a('SURF %d' % k)
        ty = atype.get(key, 'STANDARD') if A else 'STANDARD'
        a('  TYPE %s' % ty)
        c = 0.0 if s['R'] in (None, 0, 'inf') else 1.0/float(s['R'])
        a('  CURV %s 0 0 0 0 ""' % num(c, '%.12G'))
        a('  HIDE 0 0 0 0 0 0 0 0 0 0'); a('  MIRR 2 0')
        if s.get('stop') or s['i'] == 'STO': a('  STOP')
        if ty == 'EVENASPH':
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
        D = s['D']
        if isinstance(D, str):
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
        if i == fc['var_after']:
            a('  TOLE %d %s' % (fc['var_before'], num(poslen, '%.6G')))
        elif fc2 and i == fc2['var_after']:
            a('  TOLE %d %s' % (fc2['var_before'], num(poslen2, '%.6G')))
        a('  CONI %s' % num((A or {}).get('k', 0.0)))
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
    a('TOL TOFF   0   0              0              0   0 0 0 0')
    a('MNUM %d 1' % len(cfgs))
    for i, c0 in enumerate(cfgs, 1):
        a('LTTL   0   %d "%s" 0 0 0 1 1 1 0 0' % (i, ascii_(c0['name'])))
    for i, c0 in enumerate(cfgs, 1):
        d0 = c0['d0']
        v = '1.00000000E+10' if str(d0).upper().startswith('INF') else num(float(d0), '%.6G')
        a('THIC   0   %d %s 0 0 0 1 1 1 0 0' % (i, v))
    for i, c0 in enumerate(cfgs, 1):
        a('THIC  %2d   %d %s 0 0 0 1 1 1 0 0' % (fc['var_before'], i, num(c0[kb], '%.6G')))
    if fc2:
        for i, c0 in enumerate(cfgs, 1):
            a('THIC  %2d   %d %s 0 0 0 1 1 1 0 0' % (fc2['var_before'], i, num(c0[kb2], '%.6G')))
    # 链式三段浮动（如适马 85 Art / JP2018-5099 実施例4）：三个可变间隔互不成对，
    # 位置解只能钉住其中一个，剩下的第三个必须自己进 MCE。
    # zmx.mce_extra = [{"var": 27, "key": "d27"}, ...]
    for ex in zx.get('mce_extra', []):
        for i, c0 in enumerate(cfgs, 1):
            a('THIC  %2d   %d %s 0 0 0 1 1 1 0 0'
              % (ex['var'], i, num(float(c0[ex['key']]), '%.6G')))
    # ===== 逐结构渐晕：APER + FVCY/FVCX/FVDY/FVDX =====
    # Zemax 的 VDX/VDY/VCX/VCY 是**全局**量，只写在文件头里就等于所有结构共用一套。
    # 对焦镜头每个结构的渐晕差得很远（适马70微距 视场1 的 VDY 从 ∞ 的 +0.12 走到
    # 1:1 的 −0.34），不进多重结构的话近距结构的光瞳会整片被口径切掉 —— 用户实测打回过。
    # 操作数第 1 个参数是**视场号**（不是面号），第 2 个是结构号。
    vcfg = zx.get('vignetting_cfg')
    if vcfg and len(vcfg) == len(cfgs):
        for i, _c in enumerate(cfgs, 1):
            a('APER   0   %d %s 0 0 0 1 1 1 0 0' % (i, num(float(zx['fno']), '%.6G')))
        for op, ix in (('FVCY', 3), ('FVCX', 2), ('FVDY', 1), ('FVDX', 0)):
            for f in range(len(vcfg[0])):
                vals = [vcfg[ci][f][ix] for ci in range(len(cfgs))]
                if all(abs(v) < 1e-9 for v in vals):
                    continue          # 全 0 的视场行省掉：缺行时该视场直接用文件头的全局值
                for ci, v in enumerate(vals, 1):
                    a('%s  %2d   %d %s 0 0 0 1 1 1 0 0' % (op, f+1, ci, num(float(v), '%.6G')))
    a('CONF 1')
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

def verify(path, patent_f):
    txt = _read_zmx(path).splitlines()
    S=[]; cur=None; mce={}; names=[]
    for ln in txt:
        t = ln.strip()
        if t.startswith('SURF '): cur={'c':0.,'d':0.,'nd':None,'stop':False,'tole':None}; S.append(cur)
        elif cur is not None:
            if t.startswith('TYPE '): cur['type']=t.split()[1]
            elif t.startswith('XDAT '): cur['xdat']=cur.get('xdat',0)+1
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
    body=S[1:-1]
    n=1.; y=1.; u=0.
    for k,s in enumerate(body):
        n2=s['nd'] or 1.
        u=(n*u-y*s['c']*(n2-n))/n2
        if k<len(body)-1: y+=u*s['d']
        n=n2
    print('  自校验: 面数=%d 玻璃面=%d 光阑=面%s 位置解=%s' %
          (len(body), sum(1 for s in body if s['nd']),
           [i+1 for i,s in enumerate(body) if s['stop']],
           [(i+1,s['tole']) for i,s in enumerate(body) if s['tole']]))
    xa = [(i+1, s.get('xdat', 0)) for i, s in enumerate(body) if s.get('type') == 'XASPHERE']
    ev = [i+1 for i, s in enumerate(body) if s.get('type') == 'EVENASPH']
    if ev: print('  自校验: Even Asphere 面%s' % ev)
    if xa:
        bad = [i for i, n in xa if n != XDAT_NTERMS + 2]
        print('  自校验: Extended Asphere(XASPHERE) 面%s  XDAT 行数 %s%s'
              % ([i for i, _ in xa], sorted({n for _, n in xa}),
                 '  ★ 面%s 的 XDAT 条数不对' % bad if bad else '  (= 2 + 10 项，正确)'))
    print('  自校验: 反解析 EFL=%.4f (专利 %.2f)  结构=%s' % (-1/u, patent_f, names))

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('spec'); ap.add_argument('-o',required=True)
    ap.add_argument('--emb', type=int, default=0)
    ap.add_argument('--asph-type', choices=('auto', 'even', 'extended'), default='auto',
                    help='auto=有 A18/A20 的面用 Extended Asphere(XASPHERE)、其余 Even Asphere；'
                         'extended=所有非球面都用 XASPHERE；'
                         'even=全部强按 Even Asphere（高次项会被重拟合进 r^16，有残差）')
    a=ap.parse_args()
    spec=json.load(open(a.spec,encoding='utf-8'))
    emb=spec['embodiments'][a.emb]
    idx_of(emb['surfaces'])
    global _idx
    _idx = lambda s: s['_i']
    pf=dict(emb.get('general',[])).get('f (mm)', 0)
    for cat,tag in ((True,'catalog'),(False,'modelglass')):
        p='%s_%s.zmx'%(a.o,tag)
        open(p,'wb').write(build(spec,emb,cat,a.asph_type))
        print(p); verify(p, pf)

if __name__=='__main__':
    main()
