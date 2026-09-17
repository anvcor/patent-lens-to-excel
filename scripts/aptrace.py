#!/usr/bin/env python3
"""aptrace.py —— 断面图上**逐面**量净口径（沿各面自己的面型曲线走，走到墨迹断掉为止）。

为什么不用 figmeas.py 的连通域：
  figmeas 是按**元件**给一个口径的。两个问题：
  1) 一页两图、单图只有 70~80mm 宽的小幅断面图上，连通域会整组串位 —— 实测
     JP2025052870A 図1 把前片 φ62.7 读成 φ21.6，多个元件读数还雷同。
  2) 更要命的是「一个元件一个口径」本身就是错的：强弯月片的深弯那一面会被推到
     半球边上出**锋利刃口**（该篇面2 R=31.105 被给成半径 30.48，r/|R|=0.98，用户当场打回）。
     同一元件前后两面差好几 mm 是正常的 —— 前面画到玻璃外径（带法兰），
     后面的弧止于一圈平环。

做法（确定性，不依赖连通域）：
  对每个面，从 r=--start 起以 --step 步长沿 z = z_vertex + sag(r)（含非球面项）外走，
  在图上 ±--rad 像素窗口里找墨迹；连续缺 --gap 毫米就停，最后一个有墨迹的 r 就是该面画到的半径。
  只在**光轴的一侧**扫（默认下半平面）：标注引出线基本都在上方，下半干净。

用法：
  python3 aptrace.py in.pdf 41 spec.matched.cg.json -o spec.fig.json \
      --crop 2080x760+400+950 --rotate 0 --dpi 600 [--pxmm 10.4989 --x0 96.1] \
      [--clip 1010:735,999999:606] [--edge 0]

  --clip  分区截断底边：'x阈:行上限' 逗号分隔，按 x 从小到大匹配第一个满足 x<x阈 的规则。
          断面图下方的群括号（AF / A / F1 / F2 / R）会被当成元件轮廓，必须按区截掉；
          左右两半的括号高度往往不同，所以要分区给。
  --edge  机械半径 → 有効径 的磨边余量(mm)，默认 0 = 直接用量到的玻璃外径。
          **滤镜/盖板一定要用 0**：减了余量 φ 会小于像高，反而自己切掉边缘视场。

配套：量完务必再跑一次 apcap.py —— 逐面口径全固定以后 Zemax 不再替你查相邻面撞不撞，
强弯月/强非球面的相邻两面在大半径处会反向弯、空气间隙边缘厚度变成负的。
"""
import argparse, json, os, subprocess, sys, tempfile
from math import sqrt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figmeas import read_pgm, geom, autocal, crop_rotate          # 复用渲染/标定，别重复造


def build(spec, emb):
    """逐面取 R / 累计 z / 非球面系数（这里不折盖板：图上画的就是含滤镜的实物）。"""
    st = (emb.get('states') or [None])[0]
    asp = {str(x['surface']).replace('面', ''): x for x in emb.get('aspheric', [])}
    out, z = [], 0.0
    for s in emb['surfaces']:
        if s['i'] == 'IMG': break
        d = s['D']
        if isinstance(d, str): d = float(emb['variable'][d][st])
        A = asp.get(str(s['i']))
        out.append({'i': s['i'], 'R': None if s['R'] in (None, 0) else float(s['R']),
                    'z': z, 'k': (A or {}).get('k', (A or {}).get('K', 0.0)), 'A': A,
                    'stop': bool(s.get('stop')) or s['i'] == 'STO'})
        z += float(d)
    return out


def sag(s, r):
    c = 0.0 if s['R'] is None else 1.0 / s['R']
    t = 1 - (1 + s['k']) * c * c * r * r
    if t <= 0: return None
    v = c * r * r / (1 + sqrt(t))
    if s['A']:
        for key, e in (('A4', 4), ('A6', 6), ('A8', 8), ('A10', 10),
                       ('A12', 12), ('A14', 14), ('A16', 16),
                       ('A18', 18), ('A20', 20)):
            v += (s['A'].get(key) or 0.0) * r ** e
    return v


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('pdf'); ap.add_argument('page', type=int); ap.add_argument('spec')
    ap.add_argument('-o', required=True); ap.add_argument('--emb', type=int, default=0)
    ap.add_argument('--crop'); ap.add_argument('--rotate', type=int, default=0)
    ap.add_argument('--dpi', type=int, default=600)
    ap.add_argument('--pxmm', type=float); ap.add_argument('--x0', type=float)
    ap.add_argument('--axis', type=int)
    ap.add_argument('--side', choices=('lower', 'upper'), default='lower')
    ap.add_argument('--clip', default='', help="分区截断底边，如 '1010:735,999999:606'")
    ap.add_argument('--rmax', type=float, default=60.0, help='最大扫描半径(mm)；超长焦前片要调大')
    ap.add_argument('--start', type=float, default=3.0, help='起扫半径(mm)，避开光轴线')
    ap.add_argument('--step', type=float, default=0.05)
    ap.add_argument('--gap', type=float, default=0.50, help='连续缺这么多 mm 就判定曲线到头')
    ap.add_argument('--rad', type=int, default=4, help='找墨迹的像素窗口半径')
    ap.add_argument('--thr', type=int, default=170,
                    help='灰度阈值，低于它算墨迹。专利图的细线抗锯齿后灰度不低，'
                         '用 figmeas 的 140 会把面6/8/12 这类细弧扫丢')
    ap.add_argument('--edge', type=float, default=0.0, help='磨边余量(mm)')
    ap.add_argument('--img-phi', type=float, help='像面有効径(直径)，默认 2*max_y')
    a = ap.parse_args()

    spec = json.load(open(a.spec, encoding='utf-8')); emb = spec['embodiments'][a.emb]
    S = build(spec, emb)
    R, D, N, ids, z, zimg = geom(emb, spec)

    td = tempfile.mkdtemp()
    subprocess.run('pdftoppm -r %d -gray -f %d -l %d %s %s/p'
                   % (a.dpi, a.page, a.page, a.pdf, td), shell=True, check=True)
    src = [f for f in os.listdir(td) if f.endswith('.pgm')][0]
    W, H, px = crop_rotate('%s/%s' % (td, src), a.crop, a.rotate, '%s/f.pgm' % td)
    dark = [[px[y * W + x] < a.thr for x in range(W)] for y in range(H)]
    axis = a.axis if a.axis is not None else max(range(H), key=lambda y: sum(dark[y]))

    if a.pxmm and a.x0 is not None:
        pxmm, x0 = a.pxmm, a.x0
        print('标定（手动指定）: %.4f px/mm  面1顶点 x=%.1f  光轴行 y=%d' % (pxmm, x0, axis))
    else:
        cal = autocal(dark, W, H, axis, R, z, zimg, (W * 0.75) / zimg)
        if not cal: sys.exit('标定失败：交点太少，检查 --crop / --rotate / --dpi\n'
                             '（小幅断面图在 150dpi 下 r=1.5mm 只有 4 个像素，务必上 600dpi）')
        pxmm, x0, _rmm, npair, res = cal
        print('标定（各面顶点最小二乘）: %.4f px/mm  面1顶点 x=%.1f  匹配 %d 个交点  最大残差 %.1f px'
              % (pxmm, x0, npair, res))

    rules = []
    for part in filter(None, a.clip.split(',')):
        xs, ys = part.split(':'); rules.append((float(xs), int(ys)))
    rules.sort()
    def ylim(x):
        for xs, yl in rules:
            if x < xs: return yl
        return H
    sgn = 1 if a.side == 'lower' else -1
    def ink(fx, fy):
        xi, yi = int(round(fx)), int(round(fy))
        for dy in range(-a.rad, a.rad + 1):
            for dx in range(-a.rad, a.rad + 1):
                x, y = xi + dx, yi + dy
                if 0 <= x < W and 0 <= y < H and y < ylim(x) and dark[y][x]: return True
        return False

    print('面     量到半径   φ(减余量后)   备注')
    out = {}
    for s in S:
        if s['stop']: continue
        rmax, gap, r = 0.0, 0.0, a.start
        while r < a.rmax:
            g = sag(s, r)
            if g is None: break                     # 球面走到 r=|R| 就没有面型了
            if ink(x0 + (s['z'] + g) * pxmm, axis + sgn * r * pxmm):
                rmax, gap = r, 0.0
            else:
                gap += a.step
                if gap > a.gap: break
            r += a.step
        if rmax <= 0: print('  面%-4s 没扫到墨迹（检查 --crop/--side/--clip）' % s['i']); continue
        note = ''
        if s['R'] and rmax > 0.90 * abs(s['R']):
            note = '★ r/|R|=%.2f 接近半球，实物不可能磨到这里' % (rmax / abs(s['R']))
        phi = round(2 * (rmax - a.edge), 3)
        out[str(s['i'])] = phi
        print('  面%-4s %8.2f   φ%-8.2f  %s' % (s['i'], rmax, phi, note))

    ip = a.img_phi or 2 * float((spec.get('zmx') or {}).get('max_y') or 0)
    for q in emb['surfaces']:
        k = str(q['i'])
        if q['i'] == 'IMG':
            if ip: q['extra'] = {'有効径 φi': round(ip, 3)}
        elif k in out:
            q.setdefault('extra', {})['有効径 φi'] = out[k]
    spec.setdefault('zmx', {})['fix_semi_surfaces'] = [
        q['i'] for q in emb['surfaces'] if isinstance(q['i'], int) and (q.get('extra') or {}).get('有効径 φi')]
    json.dump(spec, open(a.o, 'w', encoding='utf-8'), ensure_ascii=False, indent=1)
    print('\n已写出 %s（%d 个面填了有効径，并把它们全部列入 fix_semi_surfaces）' % (a.o, len(out)))
    print('下一步：vignet.py → clearance.py → **apcap.py（固定口径必跑）** → make_zmx.py')


if __name__ == '__main__':
    main()
