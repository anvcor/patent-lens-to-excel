#!/usr/bin/env python3
"""口径体检：算每片透镜的边缘厚度、每个空气隙的边缘间隙，找出「口径撑太大 → 前后打架 / 负边厚」。
专利没给有効径时，还能反解出「至少要多少渐晕才装得下」。

    python3 clearance.py spec.matched.cg.json                    # 体检（默认 Real 光线瞄准）
    python3 clearance.py spec.matched.cg.json --aim paraxial     # 对比：近轴瞄准
    python3 clearance.py spec.matched.cg.json --solve --write    # 反解最小渐晕并写回 spec

原理
----
Zemax 的自动口径 = 所有视场、整个光瞳追迹出来的最大光线高度。**不设渐晕**时，1.0 视场
按满瞳追，边缘那几片的自动口径会被撑到远超真实值，于是出现两种物理上不可能的东西：
  · 负边缘厚度：t + z2(y) − z1(y) < 0，透镜边上厚度是负的，磨不出来。
  · 负边缘间隙：相邻两片在边缘处相撞。
真实镜头一定是有渐晕的，渐晕正是把这些口径压回可造范围的东西。所以反过来：
**能装得下的最小渐晕，就是这颗镜头渐晕量的下界**，比干脆不设渐晕靠谱得多。

sag 用完整非球面式（含 k 与 A4..A16），在各面自己的 clear semi-dia 上取值。
仅子午面（VCX 留 0），单色（各面 nd）。
"""
import json, math, argparse

AK = ('A4','A6','A8','A10','A12','A14','A16','A18','A20')


class Surf:
    def __init__(s, c, k, A, z, n, name=''):
        s.c, s.k, s.A, s.z, s.n, s.name = c, k, A, z, n, name
        s.stop = False
        s.glass = None
    def sag(s, y):
        y2 = y * y
        # r^18 / r^20 项在追迹发散时会把 y 顶到 1e15，y**20 直接 OverflowError。
        # 任何真实镜片都不会有 |y| > 1000mm，越界一律当追失。
        if not (y2 < 1e6): return None
        r = 1 - (1 + s.k) * s.c * s.c * y2
        if r < 0: return None
        z = s.c * y2 / (1 + math.sqrt(r))
        for i, a in enumerate(s.A): z += a * y ** (4 + 2 * i)
        return z
    def dsag(s, y):
        h = 1e-7
        a, b = s.sag(y + h), s.sag(y - h)
        return None if a is None or b is None else (a - b) / (2 * h)


def build(spec, emb, state):
    asph = {str(a['surface']).replace('面', ''): a for a in emb.get('aspheric', [])}
    rows = [r for r in emb['surfaces'] if r['i'] != 'IMG']
    S, z = [], 0.0
    for r in rows:
        D = r['D']
        if isinstance(D, str): D = emb['variable'][D][state]
        a = asph.get(str(r['i']))
        A = [a.get(k) or 0.0 for k in AK] if a else []
        s = Surf(0.0 if r['R'] in (None, 0) else 1.0 / float(r['R']),
                 (a or {}).get('k', 0.0) or 0.0, A, z, r.get('nd') or 1.0, str(r['i']))
        s.stop = bool(r.get('stop') or r['i'] == 'STO')
        s.glass = r.get('glass') or (('n=%.5f' % r['nd']) if r.get('nd') else None)
        s.lens = r.get('lens', '')
        S.append(s); z += float(D)
    return S, z


def trace(S, zimg, y0, z0, ang, stop_at=None):
    """子午实光线。stop_at 给面序号时，追到该面就返回其高度。"""
    dy, dz = math.sin(ang), math.cos(ang)
    y, z, n = y0, z0, 1.0
    hs = []
    for k, s in enumerate(S):
        t = tp = (s.z - z) / dz
        ok = False
        for _ in range(120):
            sg = s.sag(y + t * dy)
            if sg is None:
                t = 0.5 * (t + tp)
                if abs(t - tp) < 1e-12: break
                continue
            tn = (s.z + sg - z) / dz
            tp, t = t, tn
            ok = True
            if abs(tn - tp) < 1e-13: break
        if not ok: return hs, None
        y, z = y + t * dy, z + t * dz
        hs.append(y)
        if stop_at is not None and k == stop_at: return hs, y
        sl = s.dsag(y)
        if sl is None: return hs, None
        nx, nz = -sl, 1.0
        L = math.hypot(nx, nz); nx, nz = nx / L, nz / L
        if dy * nx + dz * nz > 0: nx, nz = -nx, -nz
        mu, ci = n / s.n, -(dy * nx + dz * nz)
        d = 1 - mu * mu * (1 - ci * ci)
        if d < 0: return hs, None
        ct = math.sqrt(d)
        dy, dz = mu * dy + (mu * ci - ct) * nx, mu * dz + (mu * ci - ct) * nz
        n = s.n
    return hs, y + (zimg - z) / dz * dy


def paraxial(S):
    n, y, u, ys = 1.0, 1.0, 0.0, []
    for k, s in enumerate(S):
        u = (n * u - y * s.c * (s.n - n)) / s.n
        ys.append(y)
        if k < len(S) - 1: y += u * (S[k + 1].z - S[k].z)
        n = s.n
    return ys, u


def pupil_geom(S, fno):
    ya, ua = paraxial(S)
    n, y, u, yb = 1.0, 0.0, 1.0, []
    for k, s in enumerate(S):
        u = (n * u - y * s.c * (s.n - n)) / s.n
        yb.append(y)
        if k < len(S) - 1: y += u * (S[k + 1].z - S[k].z)
        n = s.n
    ks = [k for k, s in enumerate(S) if s.stop][0]
    rEP = (1.0 / (2 * fno)) / abs(ua)               # 入瞳半径
    ca, cb = -yb[ks], ya[ks]
    y1, u1 = ca * ya[0] + cb * yb[0], cb
    zEP = -y1 / u1 if u1 else 0.0                   # 入瞳位置（面1顶点前为负）
    rSTO = abs(ya[ks]) * rEP / abs(ya[0])           # 近轴光阑半径
    return ks, rEP, zEP, rSTO, -1.0 / ua


class Aimer:
    """把归一化光瞳坐标 Py∈[-1,1] 变成一条实光线。"""
    def __init__(s, S, zimg, fno, mode='real', Z0=-120.0):
        s.S, s.zimg, s.mode, s.Z0 = S, zimg, mode, Z0
        s.ks, s.rEP, s.zEP, s.rSTO, s.f = pupil_geom(S, fno)

    def _shoot_from(s, y0, ang, stop_at=None):
        return trace(s.S, s.zimg, y0, s.Z0, ang, stop_at)

    def y0_paraxial(s, Py, ang):
        ye = Py * s.rEP
        return ye + (s.Z0 - s.zEP) * math.tan(ang)

    def y0_real(s, Py, ang):
        """二分：让实光线在光阑面上的高度 = Py·rSTO。"""
        target = Py * s.rSTO
        lo, hi = s.y0_paraxial(-1.6, ang), s.y0_paraxial(1.6, ang)
        flo = s._val(lo, ang)
        if flo is None:
            lo, flo = s._bracket(lo, hi, ang, +1)
            if lo is None: return None        # 该侧整条都追不到，别把 None 传给下一次 _bracket
        fhi = s._val(hi, ang)
        if fhi is None:
            hi, fhi = s._bracket(hi, lo, ang, -1)
            if hi is None: return None
        if flo is None or fhi is None: return None
        if (flo - target) * (fhi - target) > 0: return None
        for _ in range(70):
            m = 0.5 * (lo + hi); v = s._val(m, ang)
            if v is None:                      # 追失：往有效侧收
                hi = m; continue
            if (flo - target) * (v - target) <= 0: hi, fhi = m, v
            else: lo, flo = m, v
        return 0.5 * (lo + hi)

    def _val(s, y0, ang):
        hs, v = s._shoot_from(y0, ang, s.ks)
        return v if (v is not None and len(hs) > s.ks) else None

    def _bracket(s, bad, good, ang, sgn):
        for i in range(40):
            m = bad + (good - bad) * (i + 1) / 40.0
            v = s._val(m, ang)
            if v is not None: return m, v
        return None, None

    def shoot(s, Py, ang):
        y0 = s.y0_real(Py, ang) if s.mode == 'real' else s.y0_paraxial(Py, ang)
        if y0 is None: return None, None
        return s._shoot_from(y0, ang)

    def field_angle(s, Y):
        """解出实际像高为 Y 的入射角。"""
        if Y == 0: return 0.0
        lo, hi = 1e-5, math.radians(75)
        for _ in range(60):
            m = 0.5 * (lo + hi)
            hs, yi = s.shoot(0.0, m)
            if hs is None or yi is None or abs(yi) >= Y: hi = m
            else: lo = m
        return 0.5 * (lo + hi)


def semidias(S, aim, angs, vig, npup=41):
    """按 Zemax 的规则算自动 clear semi-dia：全视场 × 全光瞳的最大光线高度。"""
    sd = [0.0] * len(S)
    lost = 0
    for ang, (vdy, vcy) in zip(angs, vig):
        for i in range(npup):
            Py = -1.0 + 2.0 * i / (npup - 1)
            Pv = vdy + Py * (1.0 - vcy)
            hs, yi = aim.shoot(Pv, ang)
            if hs is None: lost += 1; continue
            if yi is None: lost += 1
            for k, h in enumerate(hs): sd[k] = max(sd[k], abs(h))
    return sd, lost


def elements(S):
    """把面序列拆成 (透镜元件, 空气隙) 段。返回 [(kind, k1, k2), ...]"""
    out = []
    for k in range(len(S) - 1):
        out.append(('glass' if S[k].n > 1.0001 else 'air', k, k + 1))
    return out


def need(kind, t, et_min, ec_min, ratio):
    """这一段至少要留多少边缘量。玻璃段同时受两条约束：
       · 绝对最小边缘厚度（磨边、装配、崩边）
       · 中心厚/边缘厚 ≤ ratio（正透镜越厚，边就得留得越多，否则磨不出也装不牢）"""
    return max(et_min, t / ratio) if kind == 'glass' else ec_min


def report(S, sd, et_min, ec_min, ratio=6.0, quiet=False):
    """边缘厚度 / 边缘间隙。y 取两面 clear semi-dia 的较大者（Zemax 的机械口径规则）。"""
    bad = []
    if not quiet:
        print('  段        面        y(mm)   中心厚   边缘值   备注')
    for kind, k1, k2 in elements(S):
        y = max(sd[k1], sd[k2])
        z1, z2 = S[k1].sag(y), S[k2].sag(y)
        t = S[k2].z - S[k1].z
        if z1 is None or z2 is None:
            bad.append((kind, k1, k2, y, t, None)); 
            if not quiet: print('  %-5s 面%2s→%2s  %7.3f  %7.3f   面形超出有效域' % (kind, S[k1].name, S[k2].name, y, t))
            continue
        e = t + z2 - z1
        lim = need(kind, t, et_min, ec_min, ratio)
        flag = '' if e >= lim else ('  ← 负!' if e < 0
                                    else '  ← 过薄 (需 %.2f, CT/ET=%.1f)' % (lim, t / e) if e > 0
                                    else '  ← 过薄')
        if e < lim: bad.append((kind, k1, k2, y, t, e))
        if not quiet:
            print('  %-5s 面%2s→%2s  %7.3f  %7.3f  %7.3f%s'
                  % ('玻璃' if kind == 'glass' else '空气', S[k1].name, S[k2].name, y, t, e, flag))
    return bad


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('spec'); ap.add_argument('--emb', type=int, default=0)
    ap.add_argument('--state', default=None)
    ap.add_argument('--aim', choices=('real', 'paraxial'), default='real')
    ap.add_argument('--et-min', type=float, default=1.25,
                    help='玻璃最小边缘厚度 mm（全画幅摄影镜头 1.0~1.5，小模组按口径缩）')
    ap.add_argument('--ec-min', type=float, default=0.30, help='空气最小边缘间隙 mm')
    ap.add_argument('--ct-et-max', type=float, default=6.0,
                    help='中心厚/边缘厚 上限（摄影镜头工艺极限约 6）')
    ap.add_argument('--solve', action='store_true', help='反解装得下的最小渐晕')
    ap.add_argument('--compare-aim', action='store_true', help='对比近轴瞄准 vs 实际瞄准')
    ap.add_argument('--write', nargs='?', const='AUTO', default=None)
    a = ap.parse_args()

    spec = json.load(open(a.spec, encoding='utf-8'))
    emb = spec['embodiments'][a.emb]
    zx = spec['zmx']; state = a.state or emb['states'][0]
    S, zimg = build(spec, emb, state)
    aim = Aimer(S, zimg, zx['fno'], a.aim)
    ymax = zx['max_y']
    fr = [1.0, 0.8, 0.6, 0.4, 0.2, 0.0]
    angs = [aim.field_angle(ymax * f) for f in fr]
    print('EFL %.4f  F/%.2f  入瞳半径 %.4f  入瞳位置 面1前 %.3f  近轴光阑半径 %.4f  瞄准=%s'
          % (aim.f, zx['fno'], aim.rEP, -aim.zEP, aim.rSTO, a.aim))
    print('视场实际半角: ' + '  '.join('%.2f°' % math.degrees(x) for x in angs))

    if a.compare_aim:
        compare_aim(S, zimg, zx, ymax, fr); return

    vg = zx.get('vignetting')
    if vg and not a.solve:
        vig0 = [(float(r[1]), float(r[3])) for r in vg][:len(fr)]
        head = '== 按 spec 里的渐晕（VCY=%s）==' % ','.join('%.3f' % v[1] for v in vig0)
    else:
        vig0 = [(0.0, 0.0)] * len(fr)
        head = '== 不设渐晕 =='
    sd, lost = semidias(S, aim, angs, vig0)
    s3 = zx.get('semi_3d')
    if s3 and not a.solve and len(s3) == len(S):
        sd = [float(v) for v in s3]
        head += '  [口径改用 vignet.py 的 3D 斜光线结果，含 VCX]'
    print('\n' + head)
    if lost: print('  有 %d 条光线追失（近轴瞄准在大视场会瞄空）' % lost)
    print('  自动 clear semi-dia: ' + ' '.join('%s=%.3f' % (S[k].name, sd[k]) for k in range(len(S))))
    bad = report(S, sd, a.et_min, a.ec_min, a.ct_et_max)
    print('  → %d 处不合格' % len(bad))
    if not a.solve:
        # 体检模式也允许把追出来的 clear semi-dia 写回 spec（layout_check 要用），
        # 但不动 vignetting —— 那是 vignet.py / --solve 的活。
        if a.write:
            out = a.spec[:-5] + '.sd.json' if a.write == 'AUTO' else a.write
            spec['zmx']['semi_diameters'] = {S[k].name: round(sd[k], 3) for k in range(len(S))}
            json.dump(spec, open(out, 'w', encoding='utf-8'), ensure_ascii=False, indent=1)
            print('  已写出 %s（只加了 semi_diameters）' % out)
        return

    print('\n== 反解渐晕（逐视场解可行光瞳区间）==')
    caps = caps_of(S, a.et_min, a.ec_min, a.ct_et_max)
    print('  各面可造口径上限（由边缘厚度/间隙/薄厚比反解）:')
    print('    ' + ' '.join('%s=%.2f' % (S[k].name, caps[k]) for k in range(len(S))))

    # 不要用「整体缩光瞳」的单参数去扫 —— 那会把前组也一起缩小，前组余量白白浪费。
    # 真实渐晕是被具体口径切出来的，光瞳既被压窄也被偏心。所以逐视场做：
    # 对每个归一化光瞳坐标 Py 追一条实光线，只要它在任何一个面上超了该面的口径上限，
    # 就判这条光线不可行；再取**最长的一段连续可行区间**，
    # 中点 → VDY、半宽 → 1−VCY。这样每个面都自动只留下它真正需要的口径。
    N = 121
    vig, sd = [], [0.0] * len(S)
    print('\n   视场Y   实际半角    VDY      VCY   通过光瞳   最先卡住的面')
    for f, ang in zip(fr, angs):
        ok, who = [], {}
        for i in range(N):
            Py = -1.0 + 2.0 * i / (N - 1)
            hs, yi = aim.shoot(Py, ang)
            good = hs is not None and yi is not None
            if good:
                for k, hgt in enumerate(hs):
                    if abs(hgt) > caps[k]:
                        good = False; who[k] = who.get(k, 0) + 1; break
            ok.append(good)
        best, cur = (0, 0, -1), None
        for i, g in enumerate(ok + [False]):
            if g and cur is None: cur = i
            elif not g and cur is not None:
                if i - cur > best[0]: best = (i - cur, cur, i - 1)
                cur = None
        n, i0, i1 = best
        if n == 0:
            print('  %6.2f   %6.2f°   全被挡' % (ymax * f, math.degrees(ang)))
            vig.append((0.0, 0.999)); continue
        lo = -1.0 + 2.0 * i0 / (N - 1); hi = -1.0 + 2.0 * i1 / (N - 1)
        vdy = (lo + hi) / 2.0; vcy = 1.0 - (hi - lo) / 2.0
        if abs(vdy) < 0.005: vdy = 0.0
        if vcy < 0.005: vcy = 0.0
        vig.append((vdy, vcy))
        for i in range(N):
            Py = -1.0 + 2.0 * i / (N - 1)
            Pv = vdy + Py * (1 - vcy)
            hs, yi = aim.shoot(Pv, ang)
            if hs:
                for k, hgt in enumerate(hs): sd[k] = max(sd[k], abs(hgt))
        top = sorted(who.items(), key=lambda q: -q[1])[:2]
        print('  %6.2f   %6.2f°  %+7.4f  %7.4f   %5.1f%%   %s'
              % (ymax * f, math.degrees(ang), vdy, vcy, 100 * (1 - vcy),
                 ', '.join('\u9762%s' % S[k].name for k, _ in top) or '-'))

    print('\n  clear semi-dia: ' + ' '.join('%s=%.3f' % (S[k].name, sd[k]) for k in range(len(S))))
    report(S, sd, a.et_min, a.ec_min, a.ct_et_max)
    if a.write:
        out = a.spec[:-5] + '.clr.json' if a.write == 'AUTO' else a.write
        spec['zmx']['vignetting'] = [[0.0, round(v, 4), 0.0, round(c, 4)] for v, c in vig]
        spec['zmx']['ray_aiming'] = 'real'
        spec['zmx']['semi_diameters'] = {S[k].name: round(sd[k], 3) for k in range(len(S))}
        json.dump(spec, open(out, 'w', encoding='utf-8'), ensure_ascii=False, indent=1)
        print('  已写出 %s' % out)


def compare_aim(S, zimg, zx, ymax, fr):
    """近轴瞄准 vs 实际瞄准：同一个归一化光瞳坐标，实际落在光阑上的哪里。"""
    ap_ = Aimer(S, zimg, zx['fno'], 'paraxial')
    ar_ = Aimer(S, zimg, zx['fno'], 'real')
    print('\n== 光线瞄准对比：近轴(Off/Paraxial) vs 实际(Real) ==')
    print('近轴光阑半径 %.4f mm。下表 = 用近轴入瞳瞄准发出的光线，实际落在光阑上的高度 / 光阑半径。'
          % ar_.rSTO)
    print('  理想应当等于 Py 本身（±1.000 / 0.000）。偏差就是光瞳像差。\n')
    print('   视场Y    半角    Py=-1.0   Py=-0.5   Py= 0.0   Py=+0.5   Py=+1.0    主光线偏差')
    for f in fr:
        Y = ymax * f
        ang = ar_.field_angle(Y)
        row, chief = [], None
        for Py in (-1.0, -0.5, 0.0, 0.5, 1.0):
            y0 = ap_.y0_paraxial(Py, ang)
            hs, v = ap_._shoot_from(y0, ang, ap_.ks)
            r = (v / ar_.rSTO) if v is not None else None
            row.append(r)
            if Py == 0.0: chief = r
        print('  %6.2f  %6.2f°  ' % (Y, math.degrees(ang))
              + '  '.join(('  ----- ' if r is None else '%+8.3f' % r) for r in row)
              + ('    —' if chief is None else '   %+.3f (%.1f%% 光阑半径)' % (chief, 100*abs(chief))))
    print('\n  → 主光线本该穿过光阑中心。近轴瞄准下它偏了多少，就是 Ray Aiming 关掉时的错误量。')
    # 满瞳能不能瞄进光阑
    print('\n  实际瞄准能否解出（Real 模式二分是否收敛）:')
    for f in fr:
        Y = ymax * f; ang = ar_.field_angle(Y)
        ok = sum(1 for i in range(21) if ar_.y0_real(-1.0 + i / 10.0, ang) is not None)
        print('    Y=%6.2f  21 个光瞳点里解出 %2d 个' % (Y, ok))


def caps_of(S, et_min, ec_min, ratio=6.0):
    """每个面「还装得下」的最大半口径：由它参与的边缘厚度/间隙约束反解，取最小。"""
    caps = [1e9] * len(S)
    for kind, k1, k2 in elements(S):
        t = S[k2].z - S[k1].z
        lim = need(kind, t, et_min, ec_min, ratio)
        # 取**满足约束的最大 y**，而不是「第一次不满足就停」——
        # 负透镜边缘比中心厚，ET 随 y 单调增，y≈0 处必然不满足绝对最小边厚，
        # 但那是中心厚度的事实、不是口径问题，不该把它的口径上限压成 0。
        ymax_ok, y = 0.0, 0.0
        while y < 60.0:
            y += 0.05
            z1, z2 = S[k1].sag(y), S[k2].sag(y)
            if z1 is None or z2 is None: break
            if t + z2 - z1 >= lim: ymax_ok = y
        caps[k1] = min(caps[k1], ymax_ok)
        caps[k2] = min(caps[k2], ymax_ok)
    return [max(c, 1e-3) for c in caps]


if __name__ == '__main__':
    main()
