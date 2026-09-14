#!/usr/bin/env python3
"""从专利断面图一步到位：自动标定 + 逐元件量机械口径 → 把有効径写回 spec。

    python3 figmeas.py in.pdf 30 spec.json -o spec.fig.json \
        --crop 900x1310+1215+1735 --rotate 90

替代「像面线倒推标定 + 画刻度线人眼读」那一整轮往返。两个关键点：

1) **标定不能拿像面线倒推**（那是循环论证：用假设的 Ymax 去定 px/mm，再用它量 Ymax）。
   这里改成在 r≈1.5 / 2.5 mm 两个高度上取各面曲线与该高度的交点，和已知的各面顶点 z
   做最小二乘拟合，直接解出 px/mm 与面1顶点 x。实测 JPWO2017138250A1：
   像面线法 8.091 px/mm（掩码重合度 73.3%），顶点拟合法 8.2058（77.9%），差 1.4%。

2) **量口径不能只看一侧、也不能取最大值**。标注文字、群括号、对焦箭头会在同一半径上留痕。
   对每个元件在光轴上下各扫一遍「最外侧暗点」的逐列剖面，各自取「众数平台」
   （没有平台时用连续性滤波后的极大值），再取两侧的较小者 —— 污染只会往外加东西，
   不会往内减，所以 min 是安全的。半径上限按该元件的近轴 |边缘光线|+|主光线| 的 1.2 倍卡。
"""
import os, sys, json, math, argparse, subprocess, tempfile


def read_pgm(path):
    d = open(path, 'rb').read(); pos, toks = 0, []
    while len(toks) < 4:
        while d[pos:pos+1].isspace(): pos += 1
        if d[pos:pos+1] == b'#':
            while d[pos] != 0x0A: pos += 1
            continue
        j = pos
        while not d[j:j+1].isspace(): j += 1
        toks.append(d[pos:j]); pos = j
    pos += 1
    w, h = int(toks[1]), int(toks[2])
    return w, h, d[pos:pos+w*h]



def crop_rotate(src, crop, rotate, dst):
    """等价于 `convert src [-crop WxH+X+Y +repage] [-rotate deg] dst`（纯 Pillow 实现）。"""
    from PIL import Image
    im = Image.open(src).convert('L')
    if crop:
        wh, _, off = crop.partition('+')
        w, h = (int(v) for v in wh.split('x'))
        x, y = (int(v) for v in off.split('+'))
        im = im.crop((x, y, min(x + w, im.width), min(y + h, im.height)))
    if rotate:
        r = int(rotate) % 360
        if r == 90:    im = im.transpose(Image.ROTATE_270)   # IM 的 +90 是顺时针
        elif r == 180: im = im.transpose(Image.ROTATE_180)
        elif r == 270: im = im.transpose(Image.ROTATE_90)
        elif r:        im = im.rotate(-r, expand=True, fillcolor=255)
    im.save(dst)
    return im.width, im.height, im.tobytes()


def geom(emb, spec):
    """展开基准状态的 R/D/nd；已补的盖板折回空气换算长（图画的是专利原始数据）。"""
    st = emb.get('states', [None])[0]
    R, D, N, ids = [], [], [], []
    for s in emb['surfaces']:
        if s['i'] == 'IMG': continue
        d = s['D']
        if isinstance(d, str): d = float(emb['variable'][d][st])
        if s.get('cover'):                    # 盖板玻璃面：折成空气换算长并入前一段
            D[-1] += float(d) / float(s['nd']); continue
        if s.get('cover_air'):
            D[-1] += float(d); continue
        R.append(None if s['R'] in (None, 0) else float(s['R']))
        D.append(float(d)); N.append(s.get('nd')); ids.append(s['i'])
    z, c = [], 0.0
    for d in D: z.append(c); c += d
    return R, D, N, ids, z, c


def clusters(dark, W, y):
    cols = [x for x in range(W) if dark[y][x]]
    if not cols: return []
    out, cur = [], [cols[0]]
    for x in cols[1:]:
        if x - cur[-1] <= 2: cur.append(x)
        else: out.append(sum(cur)/len(cur)); cur = [x]
    out.append(sum(cur)/len(cur))
    return out


def autocal(dark, W, H, axis, R, z, zimg, pxmm0):
    """在 r=1.5 / 2.5 mm 上取交点，与理论顶点 z 最小二乘拟合 (pxmm, x0)。"""
    best = None
    for rmm in (1.5, 2.5):
        rpx = int(round(rmm * pxmm0))
        if axis + rpx >= H or axis - rpx < 0: continue
        obs = sorted(set(clusters(dark, W, axis+rpx) + clusters(dark, W, axis-rpx)))
        obs = [o for i, o in enumerate(obs) if i == 0 or o - obs[i-1] > 2]
        th = [zz + (0.0 if r is None else rmm*rmm/(2*r)) for r, zz in zip(R, z)] + [zimg]
        if len(obs) < 4: continue
        # 初值：最左=面1、最右=像面
        p = (obs[-1]-obs[0]) / (th[-1]-th[0]); x0 = obs[0] - th[0]*p
        for _ in range(6):
            pair = []
            for t in th:
                xp = x0 + p*t
                o = min(obs, key=lambda v: abs(v-xp))
                if abs(o-xp) < 6: pair.append((t, o))
            if len(pair) < 4: break
            n = len(pair); sx = sum(t for t, _ in pair); sy = sum(o for _, o in pair)
            sxx = sum(t*t for t, _ in pair); sxy = sum(t*o for t, o in pair)
            den = n*sxx - sx*sx
            if abs(den) < 1e-9: break
            p = (n*sxy - sx*sy)/den; x0 = (sy - p*sx)/n
        res = max(abs(o - (x0+p*t)) for t, o in pair) if pair else 999
        if best is None or (len(pair), -res) > (best[3], -best[4]):
            best = (p, x0, rmm, len(pair), res)
    return best


def components(dark, W, H):
    """8-连通标号。返回 lab[y][x]（0=背景）与各标号的像素数。"""
    lab = [[0]*W for _ in range(H)]
    nxt, sizes = 0, {0: 0}
    from collections import deque
    for sy in range(H):
        row = dark[sy]
        for sx in range(W):
            if not row[sx] or lab[sy][sx]: continue
            nxt += 1; n = 0
            q = deque([(sy, sx)]); lab[sy][sx] = nxt
            while q:
                y, x = q.popleft(); n += 1
                for dy in (-1, 0, 1):
                    yy = y+dy
                    if yy < 0 or yy >= H: continue
                    dr = dark[yy]; lr = lab[yy]
                    for dx in (-1, 0, 1):
                        xx = x+dx
                        if 0 <= xx < W and dr[xx] and not lr[xx]:
                            lr[xx] = nxt; q.append((yy, xx))
            sizes[nxt] = n
    return lab, sizes


def seed_labels(lab, W, H, axis, x0, pxmm, R, z, ks):
    """用算出来的面型曲线当种子：在 r=2~10mm 上沿理论曲线找暗点，取其连通域号。
    光轴点划线是**断开的**，所以各元件轮廓并不连成一个域 —— 必须逐元件按曲线播种，
    不能指望「最大连通域」。标注、群括号、对焦箭头不在曲线上，自然被排除。"""
    out = set()
    for k in ks:
        r_ = R[k]
        rr = 2.0
        while rr <= 10.0:
            sag = 0.0 if r_ is None else rr*rr/(2*r_)
            x = int(round(x0 + (z[k]+sag)*pxmm))
            for side in (+1, -1):
                y = int(round(axis + side*rr*pxmm))
                for dx in range(-3, 4):
                    for dy in range(-3, 4):
                        yy, xx = y+dy, x+dx
                        if 0 <= yy < H and 0 <= xx < W and lab[yy][xx]:
                            out.add(lab[yy][xx])
            rr += 0.5
    return out


def edge_radius(lab, labs, W, H, axis, xa, xb, side, pxmm):
    best = 0
    for x in range(max(0, xa), min(W, xb+1)):
        for r in range(min(axis, H-axis-1)-1, 3, -1):
            y = axis + side*r
            if lab[y][x] in labs:
                if r > best: best = r
                break
    return best/pxmm


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('pdf'); ap.add_argument('page', type=int); ap.add_argument('spec')
    ap.add_argument('-o', required=True); ap.add_argument('--emb', type=int, default=0)
    ap.add_argument('--crop'); ap.add_argument('--rotate', type=int, default=90)
    ap.add_argument('--dpi', type=int, default=300)
    ap.add_argument('--pxmm', type=float, default=None,
                    help='手动标定：像素/mm。自动标定被标注线污染时用（配 --x0）')
    ap.add_argument('--x0', type=float, default=None, help='手动标定：面1顶点的像素 x')
    ap.add_argument('--axis', type=int, default=None,
                    help='手动指定光轴所在行。默认取暗像素最多的一行，但群括号/'
                         '箭头的长横线可能比点划线还长，那时必须手动给')
    ap.add_argument('--edge', type=float, default=None,
                    help='机械半径 → 有効径 的磨边余量(mm)。默认自动：按口径成比例 '
                         'clamp(0.02*r, 0.25, 0.80) —— 固定 0.6mm 对小口径件太狠，'
                         '实测 WO2019187633 面17(r=16.54) 被切到 15.94，比实际光束需要的 '
                         '16.06 还小，于是这个面被误判成渐晕面、把离轴上光瞳切掉了。')
    a = ap.parse_args()
    spec = json.load(open(a.spec, encoding='utf-8')); emb = spec['embodiments'][a.emb]
    R, D, N, ids, z, zimg = geom(emb, spec)

    td = tempfile.mkdtemp()
    subprocess.run('pdftoppm -r %d -gray -f %d -l %d %s %s/p' % (a.dpi, a.page, a.page, a.pdf, td),
                   shell=True, check=True)
    src = [f for f in os.listdir(td) if f.endswith('.pgm')][0]
    # 裁切/旋转改用 Pillow（本机没有 ImageMagick，而 Windows 的 convert.exe 是
    # FAT->NTFS 转换器，同名撞车，绝对不能调）。行为与 IM 一致：-crop WxH+X+Y、
    # -rotate 为顺时针角度。
    W, H, px = crop_rotate('%s/%s' % (td, src), a.crop, a.rotate, '%s/f.pgm' % td)
    dark = [[px[y*W+x] < 140 for x in range(W)] for y in range(H)]
    axis = a.axis if a.axis is not None else max(range(H), key=lambda y: sum(dark[y]))

    if a.pxmm and a.x0 is not None:
        pxmm, x0 = a.pxmm, a.x0
        print('标定（手动指定）: %.4f px/mm  面1顶点 x=%.1f  光轴行 y=%d' % (pxmm, x0, axis))
    else:
        pxmm0 = (W * 0.75) / zimg                # 只作交点扫描高度的粗尺度
        cal = autocal(dark, W, H, axis, R, z, zimg, pxmm0)
        if not cal: sys.exit('标定失败：交点太少，检查 --crop / --rotate')
        pxmm, x0, rmm, npair, res = cal
        print('标定（各面顶点最小二乘，非像面线倒推）: %.4f px/mm  面1顶点 x=%.1f  '
              '匹配 %d 个交点  最大残差 %.1f px' % (pxmm, x0, npair, res))

    lab, sizes = components(dark, W, H)

    # 元件分块（连续的玻璃段）；盖板已在 geom() 里折回空气，不会出现在这里
    blocks, cur = [], None
    for k in range(len(R)):
        if N[k]:
            cur = [k, k+1] if cur is None else [cur[0], k+1]
        elif cur is not None:
            blocks.append((cur[0], cur[1])); cur = None
    if cur is not None: blocks.append((cur[0], cur[1]))

    print('\n元件(面)        下半    上半    Δ     取值   余量   有効径 φ')
    res_map = {}; edges = {}; vals = []
    # 列窗必须按**矢高**展开，不能只取顶点到顶点：强弯月片（如索尼 135GM 的 G1，
    # R=67.6 的凸面在 r=38mm 处矢高已 11mm）边缘落在两个顶点之外，
    # 顶点窗会把整圈边缘漏掉 —— 实测面1 读成 33.5mm，真值 38.8mm。
    rmax_px = min(axis, H - axis - 1) - 4
    rmax_mm = rmax_px / pxmm
    def _sag(r_, rr):
        if not r_: return 0.0
        t = 1.0 - (rr*rr)/(r_*r_)
        return r_ * (1.0 - math.sqrt(t)) if t > 0 else r_
    for a0, b0 in blocks:
        zs = []
        for k in range(a0, b0+1):
            rl = rmax_mm if not R[k] else min(rmax_mm, abs(R[k])*0.999)
            zs += [z[k], z[k] + _sag(R[k], rl)]
        xa = int(x0 + min(zs)*pxmm) - 2; xb = int(x0 + max(zs)*pxmm) + 2
        labs = seed_labels(lab, W, H, axis, x0, pxmm, R, z, range(a0, b0+1))
        lo = edge_radius(lab, labs, W, H, axis, xa, xb, +1, pxmm)
        up = edge_radius(lab, labs, W, H, axis, xa, xb, -1, pxmm)
        v = min(lo, up)
        eg = a.edge if a.edge is not None else min(0.80, max(0.25, 0.02*v))
        for k in range(a0, b0+1): res_map[str(ids[k])] = round(v, 2)
        edges[str(ids[a0])] = eg
        vals.append((ids[a0], ids[b0], v))
        print('  面%-2s-%-3s   %6.2f  %6.2f  %+5.2f  %5.2f   %4.2f   φ%.1f'
              % (ids[a0], ids[b0], lo, up, lo-up, v, eg, 2*(v-eg)))
    print('  （两侧差 >0.4mm 说明有一侧被标注/箭头污染，已取较小者；差很大时去看一眼叠加图）')
    for (i1, j1, v1), (i2, j2, v2) in zip(vals, vals[1:]):
        if abs(v1 - v2) < 0.05:
            print('  ★ 面%s-%s 与 面%s-%s 读数相同 (%.2f) —— 两片轮廓在图上相接、'
                  '连通域被合并了，取到的是两者的外包络，小的那片会偏大' % (i1, j1, i2, j2, v1))

    n = 0
    blk_edge = {}
    for a0, b0 in blocks:
        for k in range(a0, b0+1): blk_edge[str(ids[k])] = edges[str(ids[a0])]
    for s in emb['surfaces']:
        v = res_map.get(str(s['i']))
        if v:
            s.setdefault('extra', {})['有効径 φi'] = round(2*(v-blk_edge[str(s['i'])]), 2)
            s['extra']['断面图机械半径 (mm)'] = v; n += 1
    ec = emb.setdefault('extra_columns', [])
    for c in ('有効径 φi', '断面图机械半径 (mm)'):
        if c not in ec: ec.append(c)
    emb.setdefault('general', []).append(
        ['断面图标定', '%.4f px/mm, 面1顶点 x=%.1f（各面顶点最小二乘拟合）' % (pxmm, x0)])
    json.dump(spec, open(a.o, 'w', encoding='utf-8'), ensure_ascii=False, indent=1)
    print('\n已写出 %s（%d 个面填了有効径）' % (a.o, n))


if __name__ == '__main__':
    main()
