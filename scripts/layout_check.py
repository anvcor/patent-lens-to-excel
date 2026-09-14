#!/usr/bin/env python3
"""把算出的口径画到专利自己的断面图上，肉眼核对 —— 独立于面数据的第三方证据。

    python3 layout_check.py in.pdf 37 spec.matched.cg.clr.json -o ov.png \
        --crop 950x1400+200+230 --rotate 90

断面图是**按比例作的工程图**。用图上最没有歧义的那条线（像面线）同时定原点和比例尺：
横向位置 = z = ΣD，半长 = Ymax。独立校验用「末面顶点的 z」（没参与标定）。
然后把 spec 里的 clear semi-dia 画成红点。**红点应当都落在镜片轮廓里侧一点点**：
  · 红点跑到轮廓外 → 口径算大了（多半是渐晕设少了）
  · 红点离轮廓太远 → 口径算小了（渐晕设过头）
量到的轮廓是机械半径（含磨边余量、法兰），比 clear 大 0.5~1.5mm 才正常。

自动读数会被标注文字、引出线、群括号干扰，**以看图为准**，不要盲信数字。
"""
import os, json, argparse, subprocess, tempfile



def crop_rotate(src, crop, rotate, dst):
    """等价于 `convert src [-crop WxH+X+Y +repage] [-rotate deg] dst`（纯 Pillow 实现）。"""
    from PIL import Image
    im = Image.open(src).convert('L')
    if crop:
        wh, _, off = crop.partition('+')
        wv, hv = (int(v) for v in wh.split('x'))
        xv, yv = (int(v) for v in off.split('+'))
        im = im.crop((xv, yv, min(xv + wv, im.width), min(yv + hv, im.height)))
    if rotate:
        r = int(rotate) % 360
        if r == 90:    im = im.transpose(Image.ROTATE_270)
        elif r == 180: im = im.transpose(Image.ROTATE_180)
        elif r == 270: im = im.transpose(Image.ROTATE_90)
        elif r:        im = im.rotate(-r, expand=True, fillcolor=255)
    im.save(dst)
    return im.width, im.height, im.tobytes()


def im_draw(src, dr, out):
    """执行 `-draw` 子集（line / rectangle / text），红色描边，Pillow 实现。"""
    from PIL import Image, ImageDraw
    im = Image.open(src).convert('RGB')
    g = ImageDraw.Draw(im)
    RED = (255, 0, 0)
    k = 0
    while k < len(dr):
        if dr[k] != '-draw': k += 1; continue
        cmd = dr[k + 1]; k += 2
        op, _, rest = cmd.partition(' ')
        if op == 'line':
            p1, p2 = rest.split()
            g.line([tuple(float(v) for v in p1.split(',')),
                    tuple(float(v) for v in p2.split(','))], fill=RED, width=1)
        elif op == 'rectangle':
            p1, p2 = rest.split()
            x1, y1 = (float(v) for v in p1.split(','))
            x2, y2 = (float(v) for v in p2.split(','))
            g.rectangle([x1, y1, x2, y2], outline=RED, fill=RED)
        elif op == 'text':
            pos, _, s = rest.partition(' ')
            x1, y1 = (float(v) for v in pos.split(','))
            g.text((x1, y1 - 20), s.strip().strip("'"), fill=RED)
    im.save(out)


def read_pgm(path):
    d = open(path, 'rb').read(); pos, toks = 0, []
    while len(toks) < 4:
        while d[pos:pos + 1].isspace(): pos += 1
        if d[pos:pos + 1] == b'#':
            while d[pos] != 0x0A: pos += 1
            continue
        s = pos
        while not d[pos:pos + 1].isspace(): pos += 1
        toks.append(d[s:pos])
    pos += 1; w, h = int(toks[1]), int(toks[2])
    return w, h, d[pos:pos + w * h]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('pdf'); ap.add_argument('page', type=int); ap.add_argument('spec')
    ap.add_argument('-o', default='layout_check.png')
    ap.add_argument('--crop', default=None, help='WxH+X+Y，300dpi 像素，只留断面图')
    ap.add_argument('--rotate', type=int, default=90, help='日系断面图多是竖排，顺时针 90 摆正')
    ap.add_argument('--emb', type=int, default=0); ap.add_argument('--state', default=None)
    ap.add_argument('--dpi', type=int, default=300)
    ap.add_argument('--pxmm', type=float, default=None,
                    help='手动比例尺 px/mm（自动锚点被标注骗到时用）')
    ap.add_argument('--x0', type=float, default=None, help='手动给面1顶点的像素 x')
    ap.add_argument('--axis', type=int, default=None,
                    help='手动指定光轴所在行。默认取暗像素最多的一行，但群括号/'
                         '对焦箭头的长横线可能比点划线还长（实测适马 JP2018-5099 図16），'
                         '那时红点会整体偏到图外')
    ap.add_argument('--grid', type=float, default=None, metavar='MM',
                    help='画毫米刻度线（如 2）。量断面图上各片口径最可靠的办法：'
                         '别去自动抠轮廓，画好格子用眼睛读，一次到位')
    a = ap.parse_args()

    spec = json.load(open(a.spec, encoding='utf-8'))
    emb = spec['embodiments'][a.emb]; zx = spec['zmx']
    state = a.state or emb['states'][0]
    sd = zx.get('semi_diameters') or {}
    if not sd: raise SystemExit('spec 里没有 zmx.semi_diameters —— 先跑 clearance.py --write')

    # 断面图画的是**专利原始数据**（不含我们后补的盖板）。
    # 所以补回的盖板要按空气换算长 t/n 折回去，ΣD 才和图上一致。
    z, zs = 0.0, []
    for r in emb['surfaces']:
        if r['i'] == 'IMG': break
        D = r['D']
        if isinstance(D, str): D = emb['variable'][D][state]
        D = float(D)
        if str(r['i']).startswith('CG'):
            z += D / float(r['nd']) if r.get('nd') else D
            continue
        zs.append((str(r['i']), z))
        z += D
    ztot = z

    d = tempfile.mkdtemp()
    subprocess.run(['pdftoppm', '-r', str(a.dpi), '-gray', '-f', str(a.page), '-l', str(a.page),
                    a.pdf, os.path.join(d, 'p')], check=True)
    src = os.path.join(d, [f for f in os.listdir(d) if f.startswith('p')][0])
    # 裁切/旋转/画图改用 Pillow：本机无 ImageMagick，且 Windows 的 convert.exe 是
    # FAT->NTFS 转换器，同名撞车，绝对不能调。
    pgm = os.path.join(d, 'r.pgm')
    w, h, px = crop_rotate(src, a.crop, a.rotate, pgm)

    dark = [[x for x in range(w) if px[y * w + x] < 150] for y in range(h)]
    yax = a.axis if a.axis is not None else max(range(h), key=lambda y: len(dark[y]))
    col = [[] for _ in range(w)]
    for y in range(h):
        for x in dark[y]: col[x].append(y)

    # 群括号 / 下划线是横贯的长直线，先按行密度删掉
    for y in range(h):
        if abs(y - yax) > 3 and len(dark[y]) >= 0.12 * w: dark[y] = []
    col = [[] for _ in range(w)]
    for y in range(h):
        for x in dark[y]: col[x].append(y)

    def hh(x):
        up = [yax - y for y in col[x] if y < yax - 3]
        dn = [y - yax for y in col[x] if y > yax + 3]
        return 0 if not up or not dn else min(max(up), max(dn))

    prof = [hh(x) for x in range(w)]
    on = [x for x in range(w) if prof[x] > 25]
    gs, cur = [], [on[0]]
    for x in on[1:]:
        if x - cur[-1] > 25: gs.append(cur); cur = [x]
        else: cur.append(x)
    gs.append(cur)
    # 像面线取**最后一组的最右几列**，不要取组的形心 ——
    # 盖板离像面只有 1mm 时两者会并进同一组，形心会偏到盖板那边去。
    tail = gs[-1][-5:] if len(gs[-1]) >= 5 else gs[-1]
    x_img = sum(tail) / len(tail)
    hi = max(prof[g] for g in tail)
    s = hi / float(zx['max_y'])
    x0 = x_img - ztot * s
    if a.pxmm:                       # 手动标定优先
        s = a.pxmm
        x0 = a.x0 if a.x0 is not None else x_img - ztot * s
        hi = s * float(zx['max_y'])
    print('标定: 像面线 x=%.0f 半长 %dpx = Ymax %.2f  →  %.3f px/mm；面1顶点 x=%.0f'
          % (x_img, hi, zx['max_y'], s, x0))

    # 独立校验：把面数据预测的「这一段有没有玻璃」掩码和实测轮廓比重合度。
    # 比盯单个地标稳 —— 地标会被法兰、跨轴小记号、并组骗到，整体图案不会。
    rows2 = [r for r in emb['surfaces'] if r['i'] != 'IMG']
    zg, zz, spans = [], 0.0, []
    for r in rows2:
        D = r['D']
        if isinstance(D, str): D = emb['variable'][D][state]
        D = float(D)
        if r.get('nd'): spans.append((zz, zz + D))
        zz += D
    M = [prof[x] > 25 for x in range(w)]
    hit = tot = 0
    for x in range(int(max(0, x0)), min(w, int(x0 + ztot * s))):
        z = (x - x0) / s
        p_ = any(a - 0.15 <= z <= b + 0.15 for a, b in spans)
        tot += 1; hit += (M[x] == p_)
    ok = 100.0 * hit / max(tot, 1)
    print('独立校验: 玻璃段掩码重合度 %.1f%%  %s'
          % (ok, '✓' if ok >= 72 else '← 标定可疑，检查 --crop / --rotate'))

    dr = []
    if a.grid:
        # 自动抠轮廓在有标注的图上不可靠：标注文字、引出线、群括号常常正好落在
        # 和镜片边缘同一个半径上，min(上,下)、对称性、删长横线三种滤法都会漏。
        # 画刻度线人眼读是确定性的，一次读完全部元件。
        mm = a.grid
        k = mm
        while k * s < max(yax, h - yax):
            for sg in (-1, 1):
                y = yax + sg * k * s
                if 2 < y < h - 2:
                    dr += ['-draw', 'line %d,%.0f %d,%.0f' % (int(x0) - 20, y, w - 60, y)]
                    dr += ['-draw', "text %d,%.0f '%g'" % (w - 55, y + 8, k)]
            k += mm
        print('已画 %g mm 刻度线 —— 直接读每片外缘落在哪一格，减 0.5~1mm 磨边余量即为 clear 口径' % mm)
    for i, zz in zs:
        if i not in sd: continue
        x = x0 + zz * s
        for sg in (-1, 1):
            y = yax + sg * float(sd[i]) * s
            dr += ['-draw', 'rectangle %.0f,%.0f %.0f,%.0f' % (x - 2, y - 2, x + 2, y + 2)]
    im_draw(pgm, dr, a.o)
    print('已写出 %s —— 红点 = 算出的 clear semi-dia，应当都落在镜片轮廓里侧一点点' % a.o)


if __name__ == '__main__':
    main()
