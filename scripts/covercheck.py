#!/usr/bin/env python3
"""判定专利是否偷偷剥掉了传感器盖板（把它的空气换算长折进了 BF），以及该补多厚。

    python3 covercheck.py spec.matched.json            # 只报告
    python3 covercheck.py spec.matched.json --write    # 直接把盖板补回 spec

原理
----
日系专利（尤其索尼）常把「盖板 + 到像面的空气」整段折算成空气长度并入最后一个
BF 数字。近轴上完全等效，看不出来；但盖板会在会聚光束里引入正的球差，原设计是
含盖板一起优化的，所以**把盖板拿掉之后，轴上球差反而会变大**。

于是做一件事：保持空气换算长不变（d_last' = d_last − t/n），把厚 t 的平板扫一遍，
量轴上实光线的横向像差 RMS。
  · RMS 在某个 t>0 处出现明显极小 → 专利剥了盖板，该补，厚度就是那个 t。
  · RMS 在 t=0 单调最小 → 厂家是剥掉后重新优化过的，别乱补。
平板在最后一段空气里的轴向位置不影响球差（光束角度不变），所以只需扫厚度。
纯标准库。单色（按各面 nd），只看子午轴上光线。
"""
import json, math, argparse

NG_DEFAULT = 1.51680      # BSC7 / D263 / 一般盖板玻璃
VD_DEFAULT = 64.20


class Surf:
    def __init__(s, c, k, A, z, n_after):
        s.c, s.k, s.A, s.z, s.n = c, k, A, z, n_after
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


def build(emb, state, extra_tail=None):
    """extra_tail = [(厚度, nd), ...] 追加在最后一面之后；最后一段厚度到像面。"""
    asph = {str(a['surface']).replace('面', ''): a for a in emb.get('aspheric', [])}
    rows = [dict(s) for s in emb['surfaces'] if s['i'] != 'IMG']
    ds = []
    for s in rows:
        D = s['D']
        if isinstance(D, str): D = emb['variable'][D][state]
        ds.append(float(D))
    tail = list(extra_tail or [])
    if tail:
        # 最后一段空气按换算长扣掉平板，近轴焦点不动
        ds[-1] -= sum(t / n for t, n in tail if n > 1.0)
    S, z = [], 0.0
    for s, D in zip(rows, ds):
        a = asph.get(str(s['i']))
        A = [a.get(k) or 0.0 for k in ('A4','A6','A8','A10','A12','A14','A16','A18','A20')] if a else []
        S.append(Surf(0.0 if s['R'] in (None, 0) else 1.0 / float(s['R']),
                      (a or {}).get('k', 0.0) or 0.0, A, z, s.get('nd') or 1.0))
        z += D
    for t, n in tail:
        S.append(Surf(0.0, 0.0, [], z, n)); z += t
    return S, z


def trace_axial(S, zimg, h):
    """从无限远、入射高 h 的轴上平行光，返回像面上的高度。"""
    y, z, n, dy, dz = h, S[0].z - 10.0, 1.0, 0.0, 1.0
    for s in S:
        t, tp, ok = (s.z - z) / dz, (s.z - z) / dz, False
        for _ in range(120):
            sg = s.sag(y + t * dy)
            if sg is None:
                t = 0.5 * (t + tp)
                if abs(t - tp) < 1e-13: break
                continue
            tn = (s.z + sg - z) / dz
            tp, t = t, tn
            ok = True
            if abs(tn - tp) < 1e-13: break
        if not ok: return None
        y, z = y + t * dy, z + t * dz
        sl = s.dsag(y)
        if sl is None: return None
        nx, nz = -sl, 1.0
        L = math.hypot(nx, nz); nx, nz = nx / L, nz / L
        if dy * nx + dz * nz > 0: nx, nz = -nx, -nz
        mu, ci = n / s.n, -(dy * nx + dz * nz)
        disc = 1 - mu * mu * (1 - ci * ci)
        if disc < 0: return None
        ct = math.sqrt(disc)
        dy, dz = mu * dy + (mu * ci - ct) * nx, mu * dz + (mu * ci - ct) * nz
        n = s.n
    return y + (zimg - z) / dz * dy


def ta_rms(S, zimg, rmax, nz=24):
    """轴上横向像差：按光瞳面积加权的 RMS，外加全口径值。"""
    num = den = 0.0; edge = None
    for i in range(1, nz + 1):
        h = rmax * i / nz
        y = trace_axial(S, zimg, h)
        if y is None: return None, None
        num += (y * y) * h; den += h
        edge = y
    return math.sqrt(num / den), edge


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('spec'); ap.add_argument('--emb', type=int, default=0)
    ap.add_argument('--state', default=None)
    ap.add_argument('--ng', type=float, default=NG_DEFAULT, help='盖板折射率(默认 1.5168)')
    ap.add_argument('--vd', type=float, default=VD_DEFAULT)
    ap.add_argument('--air', type=float, default=1.0, help='盖板到像面的空气(默认 1.0)')
    ap.add_argument('--tmax', type=float, default=5.0)
    ap.add_argument('--t', type=float, default=2.5,
                    help='要补的盖板厚度(默认 2.5mm，全画幅/APS-C 可换镜惯例)')
    ap.add_argument('--force-add', action='store_true',
                    help='不管球差判据结论，强制补盖板。日系可换镜专利（尤其索尼）常常是'
                         '剥掉盖板后又重新优化过，球差判据在这种情况下是查不出来的')
    ap.add_argument('--write', nargs='?', const='AUTO', default=None,
                    help='把判定出的盖板补进 spec（默认 <spec>.cg.json）')
    a = ap.parse_args()

    spec = json.load(open(a.spec, encoding='utf-8'))
    emb = spec['embodiments'][a.emb]
    state = a.state or emb['states'][0]
    fno = spec['zmx']['fno']

    rows = [s for s in emb['surfaces'] if s['i'] != 'IMG']
    last = rows[-1]
    if last.get('nd'):
        print('最后一面是玻璃，不是空气 BF —— 这篇要么已经带盖板，要么结构特殊，跳过。')
        return
    # 末端已带盖板的样子：最后一面是平面空气段，它前面那面是平面且带玻璃。
    # （别去看 rows[-3] —— 那是盖板前面那片透镜，本来就不该是平面，
    #   照原样判会漏掉「专利自己给了盖板」这种情况。）
    if len(rows) >= 2 and rows[-2].get('nd') and rows[-2]['R'] in (None, 0) \
       and rows[-1]['R'] in (None, 0):
        print('末端已有平行平板：盖板 %.2fmm nd=%.5f + 空气 %.2fmm —— 专利自己给了，无需补。'
              % (float(rows[-2]['D']), float(rows[-2]['nd']), float(rows[-1]['D'])))
        return

    S0, zi0 = build(emb, state)
    # 近轴：入瞳半径 = f/(2·Fno)
    n, y, u = 1.0, 1.0, 0.0
    for k, s in enumerate(S0):
        u = (n * u - y * s.c * (s.n - n)) / s.n
        if k < len(S0) - 1: y += u * (S0[k + 1].z - S0[k].z)
        n = s.n
    f = -1.0 / u
    rmax = f / (2.0 * fno)
    print('EFL %.4f  F/%.2f  入瞳半径 %.4f mm  扫描盖板 nd=%.5f\n' % (f, fno, rmax, a.ng))

    print('   t(mm)  空气换算   轴上 TA-RMS   全口径 TA')
    tab, best = [], (0.0, 1e30)
    for i in range(0, int(a.tmax * 20) + 1):
        t = i * 0.05
        S, zi = build(emb, state, [(t, a.ng)] if t > 0 else None)
        r, e = ta_rms(S, zi, rmax)
        if r is None: continue
        tab.append((t, r, e))
        if r < best[1]: best = (t, r)
    m = dict((round(t, 2), (r, e)) for t, r, e in tab)
    r0 = m[0.0][0]
    show = sorted(set([0.0] + [x for x in (0.5, 1.0, 1.6, 2.0, 2.5, 3.0, 4.0)
                               if x <= a.tmax] + [round(best[0], 2)]))
    for t in show:
        if t not in m: continue
        r, e = m[t]
        print('  %5.2f    %6.4f    %9.5f %s   %9.5f'
              % (t, t / a.ng, r, '*' if abs(t - best[0]) < 1e-9 else ' ', e))

    t0 = a.t
    rt = m.get(round(t0, 2), (None, None))[0]
    gain = (r0 - rt) / r0 * 100 if (rt and r0) else 0.0
    print('\n  无盖板 TA-RMS %.5f → 补 %.1fmm 后 %.5f（%s %.1f%%）'
          % (r0, t0, rt or 0, '改善' if gain > 0 else '变差', abs(gain)))
    print('  球差最优厚度 t* = %.2f mm（只作方向性参考，见下）' % best[0])
    if gain < 8.0 and not a.force_add:
        print('\n结论：球差判据看不出该补 —— 数据要么本来就不含盖板，要么**厂家剥掉后重新优化过**。')
        print('  ★ 判据的盲区：剥掉盖板后重新优化过的数据，球差已经自洽，这个判据永远查不出来。')
        print('    实测索尼 JPWO2017138250A1（Planar T* FE 50mm F1.4 ZA）就是这种：判据说别补，'
              '但实物是带盖板的，索尼故意不在专利里出现。')
        print('    所以：可换镜产品对应的日系专利，末端没有平板时默认就该补，用 --force-add。')
        return
    if gain < 8.0:
        print('\n结论：球差判据没通过，但按 --force-add 强制补回。'
              '（剥掉盖板后重新优化过的数据，球差判据查不出来 —— 见 SKILL 的盲区说明）')
    else:
        print('\n结论：把盖板补回去轴上球差明显变好，说明原设计本来就含盖板，'
              '发表时被剥掉、折进了 BF。建议补回。')
    print('  厚度取值：t* 是纯轴上单色球差的极小点，会偏大（不含带区/色差/轴外平衡），'
          '别直接拿它当物理厚度。')
    print('  实际厚度按画幅惯例定，t* 只用来确认「该补」这个方向：')
    print('    全画幅 / APS-C 可换镜（索尼、尼康、佳能专利常规）  2.5 mm 玻璃 + 1.0 mm 空气')
    print('    1 型 / M4/3                                        1.5 ~ 2.0 mm')
    print('    手机、小型模组（IR-cut）                            0.21 ~ 0.40 mm')
    eq = t0 / a.ng + a.air
    dl = float(last['D'] if not isinstance(last['D'], str)
               else emb['variable'][last['D']][state])
    print('\n  按 %.1f mm + %.1f mm 空气补回:' % (t0, a.air))
    print('    d_last %.4f → %.4f   (减去空气换算长 %.1f/%.5f + %.2f = %.4f)'
          % (dl, dl - eq, t0, a.ng, a.air, eq))
    print('    近轴焦点不动，ΣD 增加 %.4f mm' % (t0 - t0 / a.ng))
    if not a.write: return
    if isinstance(last['D'], str):
        k = last['D']; vals = set(round(float(v), 6) for v in emb['variable'][k].values())
        if len(vals) != 1:
            print('    最后一段是随状态变化的变量，请手工改 variable 表。'); return
        # BF 各状态相同（常见：D22 只是「参考用」的一列）→ 直接落成定值，
        # 否则补进来的盖板会被可变间隔逻辑带偏
        last['D'] = float(vals.pop()); emb['variable'].pop(k, None)
        for st_ in emb.get('states', []): pass
        print('    末面 D 原是各状态恒定的变量 %s，已落成定值 %.4f' % (k, last['D']))
    last['D'] = round(dl - eq, 4)
    last['note'] = (last.get('note', '') + ' ★ 专利未给传感器盖板；印刷 BF 是空气换算长，'
                    '已减去盖板换算长 %.4f' % eq).strip()
    i0 = emb['surfaces'].index(last)
    emb['surfaces'][i0 + 1:i0 + 1] = [
        {'i': 'CG1', 'R': None, 'D': t0, 'nd': a.ng, 'vd': a.vd,
         'lens': 'CG 盖板', 'group': '传感器盖板（补回）', 'type': '球面', 'extra': {},
         'cover': True, 'glass': 'HOYA BSC7',
         'note': '★ 专利未列出，按 %.1f mm 补回' % t0},
        {'i': 'CG2', 'R': None, 'D': a.air, 'lens': 'CG 盖板', 'cover_air': True,
         'group': '传感器盖板（补回）', 'type': '球面'}]
    out = a.spec[:-5] + '.cg.json' if a.write == 'AUTO' else a.write
    json.dump(spec, open(out, 'w', encoding='utf-8'), ensure_ascii=False, indent=1)
    print('    已写出 %s' % out)


if __name__ == '__main__':
    main()
