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
import json, math, argparse

class Surf:
    def __init__(s, c, k, A, z, n_after, semi, is_stop):
        s.c, s.k, s.A, s.z, s.n, s.semi, s.stop = c, k, A, z, n_after, semi, is_stop
    def sag(s, y):
        y2 = y * y
        # 定点迭代偶尔会发散（强弯月 + 负厚度的哑面序列上实测过），y 冲到 1e19 时
        # y**16 直接 OverflowError 把整个求解打断。超出任何真实镜头尺度就当追失。
        if not (y2 < 1.0e8): return None
        r = 1 - (1 + s.k) * s.c * s.c * y2
        if r < 0: return None
        z = s.c * y2 / (1 + math.sqrt(r))
        for i, a in enumerate(s.A): z += a * y ** (4 + 2 * i)
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
        A = [a.get(k) or 0.0 for k in ('A4','A6','A8','A10','A12','A14','A16')] if a else []
        phi = (s.get('extra') or {}).get('有効径 φi')
        S.append(Surf(0.0 if s['R'] in (None, 0) else 1.0/float(s['R']),
                      (a or {}).get('k', 0.0) or 0.0, A, z,
                      s.get('nd') or 1.0, (phi/2.0) if phi else None,
                      bool(s.get('stop') or s['i'] == 'STO')))
        z += float(D)
    return S, z, ((img[0].get('extra') or {}).get('有効径 φi') if img else None)

def trace(S, zimg, y0, z0, ang, apert=True):
    """从 (y0,z0) 以角 ang(rad) 出发追子午实光线。返回 (各面高度, 被挡的面序号或None, 像高)。"""
    dy, dz = math.sin(ang), math.cos(ang)
    y, z, n = y0, z0, 1.0
    hs = []
    for k, s in enumerate(S):
        t = (s.z - z) / dz                       # 先交到顶点平面
        tp = t
        ok = False
        for _ in range(120):                     # 定点迭代：t ← (顶点z + sag(y(t)) − z)/dz
            yy = y + t * dy
            sg = s.sag(yy)
            if sg is None:                       # 迭代跑出非球面有效域 → 阻尼回退
                t = 0.5 * (t + tp)
                if abs(t - tp) < 1e-12: break
                continue
            tn = (s.z + sg - z) / dz
            tp, t = t, tn
            if abs(tn - tp) < 1e-12:
                ok = True; break
            ok = True
        if not ok: return hs, k, None
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
        n = s.n
    t = (zimg - z) / dz
    return hs, None, y + t * dy

def trace3(S, zimg, P0, d, apert=True):
    """3D 斜光线。旋转对称系统：面型只依赖 r=hypot(x,y)，法线 = (-sl*x/r, -sl*y/r, 1)。
    返回 (各面处的 r, 被挡面序号或 None, 像面 y, 像面 x)。"""
    x, y, z = P0; dx, dy, dz = d
    n = 1.0; rs = []
    for k, s in enumerate(S):
        t = (s.z - z) / dz
        tp = t; ok = False
        for _ in range(120):
            r = math.hypot(x + t*dx, y + t*dy)
            sg = s.sag(r)
            if sg is None:
                t = 0.5 * (t + tp)
                if abs(t - tp) < 1e-12: break
                continue
            tn = (s.z + sg - z) / dz
            tp, t = t, tn
            ok = True
            if abs(tn - tp) < 1e-12: break
        if not ok: return rs, k, None, None
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
        n = s.n
    t = (zimg - z) / dz
    return rs, None, y + t*dy, x + t*dx




FRACS = (1.0, 0.8, 0.6, 0.4, 0.2, 0.0)
Z0 = -80.0                      # 无限远物：光线起始平面


def solve_state(S, zimg, zx, obj=None, margin=0.010, verbose=True, tag='', fit_ellipse=False):
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
            u = (n*u - y*s.c*(s.n - n)) / s.n
            ys.append(y)
            if k < len(S)-1: y += u * (S[k+1].z - S[k].z)
            n = s.n
        return ys, u
    yA, uA = par(1.0, 0.0)
    yB, _uB = par(0.0, 1.0)
    ks = [k for k, s in enumerate(S) if s.stop][0]
    rEP = (1.0/(2*zx['fno'])) / abs(uA)                 # 入瞳半径（随结构变：光阑前的 d12 在动）
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
        for span in (1.5, 4.0, 12.0):
            prev_m = prev_v = None; lo = hi = None
            N = 48
            for i in range(N+1):
                ye = -span*rEP + 2*span*rEP*i/N
                v, _ = val(ye)
                # ★ 变号只认**相邻两个有效采样**。旧写法把 prev 一直留着，于是能跨过中间
                # 一整段追失的 ye 去配一个反号点，解出一条根本不存在的"主光线"。
                # 实测 WO2021241230 Ex1：∞ 态在 35° 上这样配出 yep=3.78/像高 −3.23，
                # 外层解半视场的二分被骗着一路往大角爬，最后把 22.8° 的视场解成 35.24°，
                # 5 个离轴视场全退成 vig=[0,0,0,0]（＝满光瞳，最坏的兜底）。
                if v is None:
                    prev_m = prev_v = None; continue
                if prev_v is not None and (v == 0 or (v > 0) != (prev_v > 0)):
                    lo, hi = prev_m, ye; break
                prev_m, prev_v = ye, v
            if lo is not None: break
        if lo is None: return None, None
        for _ in range(44):
            m = 0.5*(lo+hi); v, _ = val(m)
            if v is None: break
            if v > 0: hi = m
            else: lo = m
        m = 0.5*(lo+hi); _v, yi = val(m)
        return m, yi

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
        print('%s入瞳半径 %.4f  入瞳位置 面1前 %.3f%s'
              % (tag, rEP, -zEP, '' if obj is None else '   物距 %.2f' % obj))
        print('  视场 Y′  实际半角      VDY      VCY      VCX   通过光瞳  Py/Px±1  主要卡光面')
    for fr in FRACS:
        Y = ymax*fr
        # --- 解视场参数 p：无限远解半角，有限共轭解物高 ---
        if Y == 0:
            p = 0.0
        elif obj is None:
            # 半视场角：从 0 起**逐步扫**，取第一次跨过 Y 的区间再二分。
            # 直接在 [0,70°] 上二分是错的 —— 大角上 chief 仍能解出「幽灵主光线」
            # （大部分光瞳追失、剩下一小段凑出的变号点），像高很小甚至反号，
            # 二分就被骗着一路往大角爬。实测 WO2021241230 Ex1：真值 22.8° 被解成 35.24°，
            # 于是 ∞ 态 5 个离轴视场全部退成 vig=[0,0,0,0]（＝满光瞳，最坏的兜底）。
            NP, PMAX = 140, math.radians(70)
            lo = hi = None; pm = pv = None
            for _i in range(1, NP + 1):
                q = PMAX * _i / NP
                _, yq = chief(q)
                if yq is None: pm = pv = None; continue
                yq = abs(yq)
                if pv is not None and pv < Y <= yq: lo, hi = pm, q; break
                pm, pv = q, yq
            if lo is None:
                vig.append([0.0, 0.0, 0.0, 0.0]); continue
            for _ in range(34):
                m = 0.5*(lo+hi); _, yi = chief(m)
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
                _, yi0 = chief(p_probe)
                if yi0: break
                p_probe *= 0.5
            if not yi0:
                vig.append([0.0, 0.0, 0.0, 0.0]); continue
            lo, hi = 0.0, 1.6*p_probe*Y/abs(yi0)
            for _ in range(34):
                m = 0.5*(lo+hi); _, yi = chief(m)
                if yi is None: hi = m
                elif abs(yi) < Y: lo = m
                else: hi = m
            p = 0.5*(lo+hi)
        yep, _yi = chief(p)
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
        def edge_seed(sign):
            ok_, bad_ = yseed, yep + sign*rEP
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
        hi_ok = min(hi_ok - m_ep, yep + rEP) if hi_ok < yep + rEP - 1e-9 else yep + rEP
        lo_ok = max(lo_ok + m_ep, yep - rEP) if lo_ok > yep - rEP + 1e-9 else yep - rEP
        vdy = ((hi_ok+lo_ok)/2 - yep)/rEP
        vcy = 1 - (hi_ok-lo_ok)/(2*rEP)
        vdy = 0.0 if abs(vdy) < 0.002 else vdy
        vcy = 0.0 if vcy < 0.002 else vcy
        yc = yep + vdy*rEP
        # VCX：光瞳中心高度上横向扫描（Zemax 的椭圆光瞳，X 半轴取在椭圆中心处）
        vcx = 0.0
        if shoot3(p, 0.0, yc)[1] is None:
            last, hit = 0.0, False
            for i in range(1, 201):
                px = i/200.0*rEP
                if shoot3(p, px, yc)[1] is None: last = px
                else: hit = True; break
            if hit:            # X 向没挡住任何东西时 VCX 就是 0，绝不能再减余量
                lo2, hi2 = last, min(rEP, last + rEP/200.0)
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

    vig_cfg, maxr, blk_all, axall = [], None, {}, None
    for c in cfgs:
        dmap = {}
        for k in (kb1, kb2):
            if k and k in c: dmap[k] = float(c[k])
        d0 = c.get('d0')
        obj = None if (d0 is None or str(d0).upper().startswith('INF')) else float(d0)
        S, zimg, _ = build(spec, emb, state, dmap)
        v, mr, blk, _b, axr = solve_state(S, zimg, zx, obj, a.margin, True,
                                          '\n== %s ==  ' % c.get('name', '?'), a.fit_ellipse)
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
        if m - sf.semi > 0.01: mark = '  ← ★ 光束超出口径'; bad.append(k+1)
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
        spec['zmx']['fix_semi_surfaces'] = defs
        json.dump(spec, open(out, 'w', encoding='utf-8'), ensure_ascii=False, indent=1)
        print('已写出 %s（vignetting_cfg: %d 个结构 × %d 视场）' % (out, len(vig_cfg), len(vig_cfg[0])))


if __name__ == '__main__':
    main()
