#!/usr/bin/env python3
"""apcap.py —— 固定口径的干涉体检与自动收口。

背景（用户 2026-09-08 打回）：把 clear semi-dia 按断面图量到的玻璃外径**全部固定**以后，
Zemax 不会再自己发明口径，但也不会再替你检查相邻两面会不会撞上。强弯月/强非球面的
相邻两面在大半径处可能反向弯，空气间隙的**边缘厚度变成负的**，前后透镜互相侵入。
实测 JP2025052870A：面31(R=+37.80) 与面32(R=-124.60,非球面) 轴上隔 7.53mm，
口径给到 19.45/19.55 时在 r≈18.95 处就已经穿插。

用法：  python3 apcap.py <spec.json> [-o out.json] [--ec 0.30] [--et 0.80] [--floor 0.10]
    ec    空气段要求的边缘间隙（自动放宽为 min(ec, 0.5*轴上间隔)，照顾 0.2mm 的设计薄隙）
    et    玻璃段要求的边缘厚度绝对下限
    floor 收口后必须仍比 semi_3d 光束半径大出这么多；顶不住就保留光束+floor 并打★
（CT/ET≤6 的薄厚比只报告、不用来收口 —— 大孔径厚正透镜本来就超，收了会切光束。）
"""
import json, argparse
from math import sqrt

def sag(s, r):
    R = s['R']
    c = 0.0 if R in (None, 0, 'inf') else 1.0/float(R)
    k = s.get('k', 0.0)
    t = 1 - (1+k)*c*c*r*r
    if t <= 0: return None
    z = c*r*r/(1+sqrt(t))
    for key, e in (('A4',4),('A6',6),('A8',8),('A10',10),('A12',12),('A14',14),('A16',16)):
        z += (s.get('asp') or {}).get(key, 0.0) * r**e
    return z

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('spec'); ap.add_argument('-o')
    ap.add_argument('--ec', type=float, default=0.30)
    ap.add_argument('--et', type=float, default=0.80)
    ap.add_argument('--floor', type=float, default=0.10)
    ap.add_argument('--trim', type=float, default=0.30,
                    help='把口径收到「实际光束半径 + 这个余量」。断面图量到的是**玻璃外径**'
                         '（含磨边/法兰），当净口径用普遍偏大 0.5~2mm —— 用户 2026-09 打回过。'
                         '收紧只会让口径变小，不会引入新的干涉，也不要在这之后重解渐晕'
                         '（那会越收越小地自激）。给 0 关掉。')
    ap.add_argument('--axm', type=float, default=0.05,
                    help='轴上满光瞳之上再留的余量(mm)。口径切到轴上大光瞳 = 直接把 F 数改大，'
                         '比渐晕严重得多，所以 zmx.axial_3d 是**硬底线**，收紧和收口都不许越过它。')
    ap.add_argument('--back', type=float, default=0.05,
                    help='收口后再退这么多 mm。上限有时正好落在「面不存在」的半球边上（面2 的 |R|=31.105），停在那里面型导数是无穷大，退 0.05 更稳')
    ap.add_argument('--no-fix-all', action='store_true',
                    help='不要重建 fix_semi_surfaces（默认会把所有有口径的面都固定住）')
    a = ap.parse_args()
    spec = json.load(open(a.spec, encoding='utf-8'))
    emb = spec['embodiments'][0]
    asp = {str(x['surface']).replace('面',''): x for x in emb.get('aspheric', [])}
    st = (emb.get('states') or [None])[0]
    S, z = [], 0.0
    for q in emb['surfaces']:
        if q['i'] == 'IMG': break
        D = q['D']
        if isinstance(D, str): D = float(emb['variable'][D][st])
        A = asp.get(str(q['i']))
        S.append({'i': q['i'], 'R': q['R'], 'D': float(D), 'z': z, 'nd': q.get('nd'),
                  'k': (A or {}).get('k', 0.0), 'asp': A,
                  'phi': (q.get('extra') or {}).get('有効径 φi')})
        z += float(D)
    beam = (spec.get('zmx') or {}).get('semi_3d') or []
    axl = (spec.get('zmx') or {}).get('axial_3d') or []
    for k, s in enumerate(S):
        s['beam'] = float(beam[k]) if k < len(beam) else 0.0
        s['axial'] = float(axl[k]) if k < len(axl) else 0.0
    if not axl:
        print('  ⚠ spec 里没有 zmx.axial_3d —— 先跑一次新版 vignet.py，否则查不了'
              '「口径有没有切到轴上满光瞳」这条硬底线。')
    cap = {s['i']: (s['phi']/2.0 if s['phi'] else None) for s in S}
    print('  段    面        轴上     现口径          上限     光束')
    for x, y in zip(S, S[1:]):
        if cap[x['i']] is None or cap[y['i']] is None: continue
        glass = x['nd'] is not None
        need = a.et if glass else min(a.ec, 0.5*x['D'])
        def ok(r):
            zx, zy = sag(x, r), sag(y, r)
            if zx is None or zy is None: return False
            return (y['z']+zy) - (x['z']+zx) >= need
        lo, hi = 0.5, max(cap[x['i']], cap[y['i']]) + 10
        if not ok(lo):
            print('  %s 面%-3s→%-3s 轴上就不满足（间隔 %.3f < %.3f）' %
                  ('玻璃' if glass else '空气', x['i'], y['i'], x['D'], need)); continue
        for _ in range(60):
            m = 0.5*(lo+hi)
            if ok(m): lo = m
            else: hi = m
        # 逐面比：上限要跟**这一面自己**的现口径比，不能跟两面的较小者比
        # （面1 现口径 31.75、面2 只有 26.25，拿 min 去比会漏掉面1 该收到 31.105）
        mark = ''
        for s in (x, y):
            if lo < cap[s['i']] - 1e-6:
                mark = mark or '  ← 收口'
                v = max(lo - a.back, s['beam'] + a.floor)
                if v > lo + 1e-9: mark = '  ★ 上限低于光束+%.2f，保留光束余量' % a.floor
                cap[s['i']] = min(cap[s['i']], v)
        print('  %s 面%-3s→%-3s %7.3f  %6.3f/%6.3f  %8.3f  %6.2f/%6.2f%s' %
              ('玻璃' if glass else '空气', x['i'], y['i'], x['D'],
               cap[x['i']], cap[y['i']], lo, x['beam'], y['beam'], mark))
    # ---- 按实际光束收紧 ----
    # 断面图画的是玻璃外径，含磨边和法兰；净口径应当只比实际通光光束大一点点。
    # 注意：**收紧之后绝不能再跑一次 vignet** —— 口径是渐晕解的输入，
    # 收小 → 光束变小 → 再收更小，会自激塌缩。收紧只是最后一道整形。
    if a.trim > 0:
        cut = []
        for s in S:
            if cap[s['i']] is None or not s['beam']: continue
            v = max(s['beam'] + a.trim, s['beam'] + 0.02)
            v = max(v, s['axial'] + a.axm)          # 轴上满光瞳是硬底线，收不到它下面
            if v < cap[s['i']] - 1e-6:
                cut.append((s['i'], cap[s['i']], v, s['beam'])); cap[s['i']] = v
        if cut:
            print('\n按光束收紧（口径 = 光束 + %.2fmm）：' % a.trim)
            for i, old, new, bm in cut:
                print('  面%-4s %6.3f → %6.3f  (光束 %6.3f，缩 %.3f)' % (i, old, new, bm, old-new))
            print('  共收紧 %d 个面；其余各面光束已经贴到玻璃外径，再收就要切光了' % len(cut))
    # ---- 轴上满光瞳硬底线：谁都不许切它，切了就是把相对孔径改掉 ----
    if axl:
        up = []
        for s in S:
            if cap[s['i']] is None or not s['axial']: continue
            need = s['axial'] + a.axm
            if cap[s['i']] < need - 1e-6:
                up.append((s['i'], cap[s['i']], need, s['axial'])); cap[s['i']] = need
        if up:
            print('\n★★ 下列面的口径切到了**轴上满光瞳**（切了 F 数就变大），已放大到需求值：')
            for i, old, new, ax in up:
                print('  面%-4s %6.3f → %6.3f  (轴上需求 %6.3f，原来切掉 %.1f%% 半径)'
                      % (i, old, new, ax, 100*(1-old/ax)))
            print('  注：断面图量到的玻璃外径在这几个面上偏小（图按 ∞ 态画，近距态轴上锥更粗），'
                  '以轴上需求为准；若几何干涉上限也顶不住，那是真冲突，要回去核数据。')
    n = 0
    for q in emb['surfaces']:
        if q['i'] in cap and cap[q['i']]:
            new = round(2*cap[q['i']], 3)
            if abs(new - (q['extra'] or {}).get('有効径 φi', new)) > 1e-6:
                q['extra']['有効径 φi'] = new; n += 1
    print('\n收口 %d 个面' % n)
    # 口径全固定是 apcap 的前提（不然 Zemax 自己算的口径会绕过这里的干涉检查），
    # 所以在这里**重建** fix_semi_surfaces —— vignet.py --write 会把它覆盖成
    # 「渐晕定义面」那一小串，apcap 必须排在它后面并把全部有口径的面重新列进去。
    if not a.no_fix_all:
        spec.setdefault('zmx', {})['fix_semi_surfaces'] = [
            q['i'] for q in emb['surfaces']
            if isinstance(q['i'], int) and (q.get('extra') or {}).get('有効径 φi')]
        print('已把 %d 个面列入 fix_semi_surfaces（DIAM 固定）'
              % len(spec['zmx']['fix_semi_surfaces']))
    if a.o:
        json.dump(spec, open(a.o, 'w', encoding='utf-8'), ensure_ascii=False, indent=1)
        print('已写出', a.o)

if __name__ == '__main__':
    main()
