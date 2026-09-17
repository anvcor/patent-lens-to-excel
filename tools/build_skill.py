#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""build_skill.py —— 把本仓库打成可上传到 claude.ai 的 .skill 包。

    python tools/build_skill.py                       # → E:\\Download\\patent-lens-to-excel_YYYYMMDD.skill
    python tools/build_skill.py -o D:\\x.skill
    python tools/build_skill.py --list                 # 只列出会打进去的文件

规则（都踩过坑）：
* .skill 就是 zip，**顶层必须是同名文件夹** `patent-lens-to-excel/`。
* 条目路径必须是**正斜杠**。PowerShell 的 Compress-Archive 写反斜杠，上传后解不开 —— 别用它。
* 只打 skill 运行需要的东西：SKILL.md、scripts/、references/、assets/。
  DEVELOPMENT.md / README.md / tools/ / docs/ / .git 都是开发用的，不进包。
* __pycache__、*.pyc、*.bak* 一律排除。
"""
import argparse, datetime, os, sys, zipfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
NAME = 'patent-lens-to-excel'
INCLUDE = ['SKILL.md', 'scripts', 'references', 'assets']


def collect():
    out = []
    for item in INCLUDE:
        full = os.path.join(ROOT, item)
        if os.path.isfile(full):
            out.append((item, full)); continue
        if not os.path.isdir(full):
            sys.exit('!! 缺少 %s' % item)
        for dp, dns, fns in os.walk(full):
            dns[:] = sorted(d for d in dns if d != '__pycache__')
            for fn in sorted(fns):
                if fn.endswith('.pyc') or '.bak' in fn:
                    continue
                p = os.path.join(dp, fn)
                out.append((os.path.relpath(p, ROOT).replace(os.sep, '/'), p))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('-o', default=os.path.join(r'E:\Download', '%s_%s.skill'
                                               % (NAME, datetime.date.today().strftime('%Y%m%d'))))
    ap.add_argument('--list', action='store_true')
    a = ap.parse_args()
    files = collect()
    if a.list:
        for arc, _ in files: print(arc)
        print('%d files' % len(files)); return
    dirs = sorted({'/'.join(arc.split('/')[:i]) for arc, _ in files for i in range(1, arc.count('/') + 1)})
    with zipfile.ZipFile(a.o, 'w', zipfile.ZIP_DEFLATED) as z:
        z.writestr(zipfile.ZipInfo(NAME + '/'), '')
        for d in dirs:
            z.writestr(zipfile.ZipInfo('%s/%s/' % (NAME, d)), '')
        for arc, src in files:
            z.write(src, '%s/%s' % (NAME, arc))
    names = zipfile.ZipFile(a.o).namelist()
    assert all('\\' not in n for n in names), '出现反斜杠路径'
    assert all(n.startswith(NAME + '/') for n in names), '顶层不是同名文件夹'
    npy = sum(1 for n in names if n.endswith('.py'))
    print('已写出 %s：%d 个条目（%d 个 .py），%d 字节' % (a.o, len(names), npy, os.path.getsize(a.o)))


if __name__ == '__main__':
    main()
