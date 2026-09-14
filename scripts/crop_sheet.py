#!/usr/bin/env python3
"""把若干页的若干条带拼成 1~2 张「刚好还看得清」的图，把 5 次视觉读压到 2 次。

    python3 crop_sheet.py in.pdf -o t --panels 20:250:1900 20:2000:1500 21:250:1750

页码:起始y:高度（y 按 300dpi 的像素）。每张输出宽 1500、最高 2000px —— 这是 Read 不降采样的上限，
超过就自动切成 t1 / t2 …，保证 5 位小数不被糊掉。只依赖 poppler + ImageMagick。
"""
import os, argparse, subprocess, tempfile
W, HMAX = 1500, 2000

def ident(path):
    out = subprocess.run(['identify', '-format', '%w %h', path],
                         capture_output=True, text=True, check=True).stdout.split()
    return int(out[0]), int(out[1])

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('pdf'); ap.add_argument('-o', default='t')
    ap.add_argument('--dpi', type=int, default=300)
    ap.add_argument('--panels', nargs='+', required=True, help='page:y:h')
    a = ap.parse_args()
    d = tempfile.mkdtemp()
    for p in sorted({int(s.split(':')[0]) for s in a.panels}):
        subprocess.run(['pdftoppm', '-r', str(a.dpi), '-png', '-f', str(p), '-l', str(p),
                        a.pdf, os.path.join(d, 'p%d' % p)], check=True)
    strips = []
    for i, spec in enumerate(a.panels):
        p, y, h = (int(x) for x in spec.split(':'))
        src = os.path.join(d, [f for f in os.listdir(d) if f.startswith('p%d-' % p)][0])
        out = os.path.join(d, 's%02d.png' % i)
        subprocess.run(['convert', src, '-crop', '%dx%d+0+%d' % (10 ** 5, h, y), '+repage',
                        '-resize', '%dx' % W, out], check=True)
        strips.append((spec, out, ident(out)[1]))
    # first-fit：条带按原顺序装箱，装不下就放进前面还有余量的那张，尽量少出图
    sheets, used = [], []
    for s in strips:
        for k in range(len(sheets)):
            if used[k] + s[2] <= HMAX:
                sheets[k].append(s); used[k] += s[2]; break
        else:
            sheets.append([s]); used.append(s[2])
    for i, grp in enumerate(sheets, 1):
        name = '%s%d.png' % (a.o, i)
        subprocess.run(['convert'] + [g[1] for g in grp] +
                       ['-background', 'white', '-append', name], check=True)
        w, h = ident(name)
        print('%s  %dx%d  条带: %s' % (name, w, h, ', '.join(g[0] for g in grp)))

if __name__ == '__main__':
    main()
