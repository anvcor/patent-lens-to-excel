#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""which_example.py —— 「这篇专利的哪个実施例是实物？」查表器

PhotonsToPhotos 的 Optical Bench Hub 已经替 1280+ 支量产镜头标好了
「产品 ↔ 专利号 ↔ 实施例号」。开工第一件事查这里，比自己从 f/Fno/片数
一个个排除快得多，也是独立于自己判断的旁证。

用法（本机 python 一律用 /e/Python/python 或 E:\\Python\\python.exe）：

    python which_example.py --refresh            # 抓一次官网，刷新本地缓存
    python which_example.py RF24mm               # 按产品名找
    python which_example.py US20220019061        # 按专利号找
    python which_example.py sigma 105            # 多个词 = AND

缓存文件 p2p_index.json 就放在本脚本旁边。没网时直接用缓存。

注意（skill 的原话）：Bill Claff 的处方**不要当数据源**（他用模型玻璃、
口径是他自己算的），但拿来交叉校验 R / d / 非球面系数极好用。
他的 .txt 处方可以直接取：
    https://www.photonstophotos.net/GeneralTopics/Lenses/OpticalBench/Data/<PatentID>_<ExampleID>.txt
里面有 [variable distances]（含 d0 / Magnification / Aperture Diameter）、
[lens data]（R / d / nd / 有効径 / vd）、[aspherical data]、[element data]、[group data]。
"""
import json
import os
import re
import sys
import urllib.request

HUB = ('https://www.photonstophotos.net/GeneralTopics/Lenses/'
       'OpticalBench/OpticalBenchHub.htm')
DATA = ('https://www.photonstophotos.net/GeneralTopics/Lenses/'
        'OpticalBench/Data/%s_%s.txt')
HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, 'p2p_index.json')


def refresh():
    req = urllib.request.Request(HUB, headers={'User-Agent': 'Mozilla/5.0'})
    html = urllib.request.urlopen(req, timeout=60).read().decode('utf-8', 'replace')
    pat = re.compile(r"PatentID:\s*'([^']*)'\s*,\s*ExampleID:\s*'([^']*)'"
                     r"\s*,\s*TaleID:\s*'([^']*)'\s*,\s*Description:\s*'([^']*)'")
    rows = [{'patent': a, 'example': b, 'tale': c, 'product': d}
            for a, b, c, d in pat.findall(html)]
    if not rows:
        sys.exit('解析失败：页面结构变了，去看一眼 OpticalBenchHub.htm 的 JS 数组')
    json.dump(rows, open(CACHE, 'w', encoding='utf-8'), ensure_ascii=False, indent=0)
    print('已缓存 %d 条 → %s' % (len(rows), CACHE))
    return rows


def load():
    if not os.path.exists(CACHE):
        print('本地无缓存，先抓一次…')
        return refresh()
    return json.load(open(CACHE, encoding='utf-8'))


def norm(s):
    return re.sub(r'[^a-z0-9]', '', s.lower())


def main():
    args = [a for a in sys.argv[1:] if a != '--refresh']
    if '--refresh' in sys.argv[1:]:
        rows = refresh()
        if not args:
            return
    else:
        rows = load()
    if not args:
        print(__doc__)
        return
    keys = [norm(a) for a in args]
    hit = [r for r in rows
           if all(k in norm(r['product'] + ' ' + r['patent']) for k in keys)]
    if not hit:
        print('没找到。换个写法试试（产品名去掉空格也行），'
              '或者 --refresh 再抓一次；Hub 只收录了有对应量产镜头的那部分专利。')
        return
    for r in hit:
        print('%-46s %-18s %s' % (r['product'], r['patent'], r['example']))
    if len(hit) == 1:
        r = hit[0]
        print('\n交叉校验用的处方（R/d/非球面可信，玻璃与口径不可信）：')
        print('  ' + DATA % (r['patent'], r['example']))


if __name__ == '__main__':
    main()
