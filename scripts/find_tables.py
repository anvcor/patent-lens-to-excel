#!/usr/bin/env python3
"""无文字层专利：不靠视觉、不靠联系表，直接按「贯通横线密度」定位表格页 / 附图页。

    python3 find_tables.py in.pdf                 # 表格候选页
    python3 find_tables.py in.pdf --figures       # 附图（像差图/断面图）候选页

只依赖 poppler 的 pdftoppm，纯标准库解析 PGM，无需 numpy/pillow。
119 页 30dpi 全本扫描约 5 秒，替代 pdftoppm -r40 + montage + 3~4 次视觉读联系表。
"""
import sys, os, glob, argparse, subprocess, tempfile

def read_pgm(path):
    """读 P5 二进制 PGM，返回 (w, h, bytes)。pdftoppm -gray 的输出就是它。"""
    with open(path, 'rb') as f:
        data = f.read()
    # 头部：P5 <w> <h> <maxval>，token 之间可有注释行
    pos, toks = 0, []
    while len(toks) < 4:
        while pos < len(data) and data[pos:pos + 1].isspace(): pos += 1
        if data[pos:pos + 1] == b'#':
            while pos < len(data) and data[pos] != 0x0A: pos += 1
            continue
        s = pos
        while pos < len(data) and not data[pos:pos + 1].isspace(): pos += 1
        toks.append(data[s:pos])
    pos += 1
    w, h = int(toks[1]), int(toks[2])
    return w, h, data[pos:pos + w * h]

def scan(pdf, dpi=30, thresh=128, fill=0.35):
    d = tempfile.mkdtemp()
    subprocess.run(['pdftoppm', '-r', str(dpi), '-gray', pdf, os.path.join(d, 'p')], check=True)
    rows = []
    for f in sorted(glob.glob(os.path.join(d, 'p-*.pgm'))):
        w, h, px = read_pgm(f)
        need = int(w * fill)
        hlines = ink = 0
        colcnt = [0] * w
        for y in range(h):
            line = px[y * w:(y + 1) * w]
            dark = 0
            for x, v in enumerate(line):
                if v < thresh:
                    dark += 1; colcnt[x] += 1
            ink += dark
            if dark > need: hlines += 1
        vneed = int(h * fill)
        rows.append({'page': int(os.path.basename(f).split('-')[1].split('.')[0]),
                     'hlines': hlines,
                     'vlines': sum(1 for c in colcnt if c > vneed),
                     'ink': ink / float(w * h)})
        os.remove(f)
    return rows

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('pdf'); ap.add_argument('--dpi', type=int, default=30)
    ap.add_argument('--top', type=int, default=15)
    ap.add_argument('--locate', default=None,
                    help='给表格标题正则（如 "TABLE 1$" / "【表１】"），一次算出页码/分栏/裁图参数')
    ap.add_argument('--figures', action='store_true')
    a = ap.parse_args()
    if a.locate:
        return locate(a.pdf, a.locate)
    rows = scan(a.pdf, a.dpi)
    if a.figures:
        rows.sort(key=lambda r: -(r['ink'] * 100 - r['hlines']))
        print('附图候选页:')
    else:
        rows.sort(key=lambda r: (-r['hlines'], -r['ink']))
        print('表格候选页:')
    for r in rows[:a.top]:
        print('  p%-4d 横线=%-4d 竖线=%-4d 墨色=%.3f' % (r['page'], r['hlines'], r['vlines'], r['ink']))
    print('\n下一步: python3 crop_sheet.py %s -o t --panels <页>:<y>:<高> ...' % a.pdf)


# ---------------------------------------------------------------- locate 模式
def locate(pdf, pat, dpi=300):
    """从表格标题一路算到裁图参数：定页 → 分栏 → 找内容带 → 打印 convert 命令。
    原来这一段要拆成 6~7 次 bash（pdftotext / 找页眉 / 切页 / 渲染 / 扫列 / 扫行 / 裁图），
    现在一次跑完。美国专利是双栏的，日系是单栏，靠列墨色自动判。"""
    import re, subprocess, tempfile, os
    d = tempfile.mkdtemp()
    txt = os.path.join(d, 'p.txt')
    subprocess.run(['pdftotext', '-layout', pdf, txt], check=True)
    pages = open(txt, encoding='utf-8', errors='replace').read().split('\f')
    rx = re.compile(pat, re.M)
    hits = [i + 1 for i, pg in enumerate(pages) if rx.search(pg)]
    if not hits:
        print('正文里没找到 %r —— 可能整本无文字层，改用横线密度排名' % pat); return
    pno = hits[0]
    print('「%s」出现在 PDF 第 %s 页（取第一处 %d）' % (pat, hits[:6], pno))
    subprocess.run(['pdftoppm', '-r', str(dpi), '-gray', '-f', str(pno), '-l', str(pno),
                    pdf, os.path.join(d, 'g')], check=True)
    f = os.path.join(d, [x for x in os.listdir(d) if x.startswith('g')][0])
    w, h, px = read_pgm(f)
    colink = [sum(1 for y in range(0, h, 6) if px[y * w + x] < 128) for x in range(w)]
    # 分栏：连续 >40px 的空白列即栏间隙
    # 栏间隙不会是「一个墨点都没有」（扫描噪点、页码、装订线），
    # 按最大列墨色的 2% 当空白阈值，否则双栏的美国专利会被并成一栏。
    thr = 0.02 * max(colink) if colink else 0
    cols, run, st = [], 0, None
    for x in range(w):
        if colink[x] <= thr:
            run += 1
            if st is not None and run > 30: cols.append((st, x - run)); st = None
        else:
            run = 0
            if st is None: st = x
    if st is not None: cols.append((st, w - 1))
    cols = [c for c in cols if c[1] - c[0] > 200]
    print('版面: %d 栏 %s' % (len(cols), cols))
    for ci, (x0, x1) in enumerate(cols, 1):
        bands, run, st = [], 0, None
        for y in range(h):
            ink = sum(1 for x in range(x0, x1, 3) if px[y * w + x] < 128)
            if ink < 2:
                run += 1
                if st is not None and run > 45: bands.append((st, y - run)); st = None
            else:
                run = 0
                if st is None: st = y
        if st is not None: bands.append((st, h - 1))
        bands = [b for b in bands if b[1] - b[0] > 60]
        print('  第%d栏 x %d-%d 内容带:' % (ci, x0, x1))
        for (ya, yb) in bands:
            hgt = yb - ya + 20
            # Read 不降采样的上限约 1500x2000；超高的带自动切两段并留 80px 重叠，
            # 否则接缝会把表格的一行撕成两半。
            segs = [(ya, yb)] if hgt <= 1900 else \
                   [(ya, ya + (yb - ya) // 2 + 40), (ya + (yb - ya) // 2 - 40, yb)]
            for k, (sa, sb) in enumerate(segs, 1):
                print('    convert <page.pgm> -crop %dx%d+%d+%d +repage -normalize out%s.png'
                      % (x1 - x0 + 20, sb - sa + 20, x0 - 10, sa - 10,
                         '' if len(segs) == 1 else str(k)))
    print('\n渲染该页: pdftoppm -r %d -gray -f %d -l %d %s pg' % (dpi, pno, pno, pdf))

if __name__ == '__main__':
    main()
