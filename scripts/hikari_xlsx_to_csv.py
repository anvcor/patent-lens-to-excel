#!/usr/bin/env python3
"""HIKARI 官方目录 xlsx → glasslib 用的 nd/vd CSV。

    python3 hikari_xlsx_to_csv.py HIKARI_ALL_Catalog_Data.xlsx -o ../assets/glass/HIKARI.csv

HIKARI（住田光学/光ガラス）只发 Excel 目录、不发 .AGF，所以单独转一道。
Data 工作表的表头分三行，真正的列名在第 3 行；按列名定位而不是写死列号，
这样官方调整列序时不至于静默取错数。
"""
import argparse, csv, sys
import openpyxl


def convert(xlsx, out):
    wb = openpyxl.load_workbook(xlsx, read_only=True, data_only=True)
    ws = wb['Data'] if 'Data' in wb.sheetnames else wb[wb.sheetnames[0]]
    rows = ws.iter_rows(values_only=True)
    head = [next(rows) for _ in range(3)]          # 第3行是列名行
    names = [str(x).split('\n')[0].strip() if x is not None else '' for x in head[2]]

    def find(*cands):
        for i, n in enumerate(names):
            if n in cands:
                return i
        raise SystemExit('找不到列: %s\n实际列名前30个: %s' % (cands, names[:30]))

    i_nd, i_vd = find('d'), find('νd', 'vd', 'nud')
    i_name = 0                                      # 第1列固定是 硝種 Glass type

    out_rows = []
    for r in rows:
        if not r or not r[i_name]:
            continue
        nm = str(r[i_name]).strip()
        try:
            nd, vd = float(r[i_nd]), float(r[i_vd])
        except (TypeError, ValueError):
            continue
        if not (1.0 < nd < 2.5 and 5 < vd < 120):   # 挡掉表尾说明行
            continue
        out_rows.append((nm, nd, vd))

    with open(out, 'w', newline='', encoding='utf-8') as f:
        w = csv.writer(f)
        w.writerow(['name', 'nd', 'vd'])
        w.writerows(out_rows)
    return len(out_rows)


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('xlsx')
    ap.add_argument('-o', '--out', default='HIKARI.csv')
    a = ap.parse_args()
    n = convert(a.xlsx, a.out)
    print('已写出 %d 种玻璃 -> %s' % (n, a.out))
