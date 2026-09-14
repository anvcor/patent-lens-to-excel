#!/usr/bin/env python3
"""玻璃目录加载与牌号匹配。

用法（命令行，快速查一批 nd/vd）:
    python3 glasslib.py --vendors HOYA OHARA --pairs 1.48749/70.44 1.883/40.81

作为库:
    from glasslib import load, match, grade, tolerances
"""
import os, re, sys, csv, json, argparse

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_DIR = os.path.join(HERE, '..', 'assets', 'glass')
MOLD = re.compile(r'^(D-|M-|MP-|MC-|Q-)')       # 模压/精密模造玻璃前缀
# CDGM: D- ／ HOYA: M- MP- MC- ／ HIKARI: Q-（J- 为常规抛光料）
# 注意用带连字符的前缀匹配，避免误伤 CDGM 的 QF50 / QK3L 等常规牌号。


def _decode(raw: bytes) -> str:
    """AGF 有 UTF-16LE(带BOM) 和 ASCII 两种，逐个试。"""
    if raw[:2] in (b'\xff\xfe', b'\xfe\xff'):
        return raw.decode('utf-16')
    for enc in ('utf-8', 'latin-1'):
        try:
            t = raw.decode(enc)
            if t.count('NM ') > 10:
                return t
        except UnicodeDecodeError:
            pass
    return raw.decode('utf-16', errors='ignore')


def _find_file(vendor: str, glass_dir: str = None):
    """按 目录优先级 × 扩展名优先级 找目录文件。"""
    dirs = [glass_dir, os.environ.get('PATENT_GLASS_DIR'), DEFAULT_DIR]
    for d in dirs:
        if not d:
            continue
        for ext in ('.AGF', '.agf', '.csv', '.CSV'):
            cand = os.path.join(d, vendor.upper() + ext)
            if os.path.exists(cand):
                return cand
    raise FileNotFoundError('找不到玻璃目录文件: %s(.AGF/.csv)，查找路径 %s'
                            % (vendor.upper(), [d for d in dirs if d]))


def _load_agf(path, vendor):
    """NM 行格式统一为: NM <名称> <色散公式> <MIL> <nd> <vd> ...
    OHARA / HOYA / CDGM 的字段位置一致。"""
    out = []
    for line in _decode(open(path, 'rb').read()).splitlines():
        if not line.startswith('NM '):
            continue
        p = line.split()
        try:
            out.append({'name': p[1], 'nd': float(p[4]), 'vd': float(p[5]),
                        'vendor': vendor})
        except (IndexError, ValueError):
            continue
    return out


def _load_csv(path, vendor):
    """三列 CSV: name,nd,vd（表头必须有这三个名字，多余列忽略）。
    给只发 Excel 目录、不发 .AGF 的厂家用（如 HIKARI）。"""
    out = []
    with open(path, encoding='utf-8-sig', newline='') as f:
        for row in csv.DictReader(f):
            try:
                out.append({'name': row['name'].strip(),
                            'nd': float(row['nd']), 'vd': float(row['vd']),
                            'vendor': vendor})
            except (KeyError, TypeError, ValueError, AttributeError):
                continue
    return out


def load(vendor: str, glass_dir: str = None):
    """读一家的目录，返回 [{'name','nd','vd','vendor'}, ...]。

    支持 Zemax .AGF（OHARA / HOYA / CDGM）和三列 CSV（HIKARI 等只发 Excel 的厂家）。
    """
    vendor = vendor.upper()
    path = _find_file(vendor, glass_dir)
    out = (_load_csv(path, vendor) if path.lower().endswith('.csv')
           else _load_agf(path, vendor))
    if not out:
        raise RuntimeError('未解析到玻璃: ' + path)
    return out


def tolerances(nd_values):
    """按专利给出的 nd 位数决定匹配容差。

    专利常把 nd 降精度印刷；容差必须跟着放宽，否则会把明显对不上的牌号
    评成"精确匹配"，反过来又会把真正的原设计牌号排到后面。
    """
    dec = 0
    for v in nd_values:
        s = ('%.6f' % v).rstrip('0')
        dec = max(dec, len(s.split('.')[1]) if '.' in s else 0)
    if dec >= 4:
        return 0.0005, 0.05
    if dec == 3:
        return 0.002, 0.15
    return 0.006, 0.25


def match(nd, vd, lib, tol=(0.0005, 0.05), k=3, aspheric=False):
    """加权最近邻。nd 权重取容差量级，vd 同理。

    球面件默认给模压牌号扣分（要抛光料）；非球面件不扣，因为模压料正是
    非球面镜的常规选择。
    """
    tn, tv = tol
    def key(g):
        s = ((g['nd'] - nd) / tn) ** 2 + ((g['vd'] - vd) / tv) ** 2
        if MOLD.match(g['name']) and not aspheric:
            s += 60.0
        return s
    return sorted(lib, key=key)[:k]


def grade(g, nd, vd, tol=(0.0005, 0.05)):
    """给替代品一个人能一眼看懂的评级。"""
    tn, tv = tol
    dn, dv = abs(g['nd'] - nd), abs(g['vd'] - vd)
    if dn <= tn * 0.1 and dv <= tv:
        return '优·完全一致'
    if dn <= tn and dv <= tv * 4:
        return '良'
    if dn <= tn * 10 and dv <= tv * 10:
        return '可'
    return '差·需重优化'


def exact_vendor(nd, vd, libs, tol=(0.0005, 0.05)):
    """哪几家能完全对上 —— 用来推定专利原设计用的是谁家的料。"""
    hits = []
    for lib in libs:
        g = match(nd, vd, lib, tol, k=1)[0]
        if grade(g, nd, vd, tol).startswith('优'):
            hits.append('%s %s' % (g['vendor'], g['name']))
    return ' / '.join(hits)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--vendors', nargs='+', default=['HOYA', 'OHARA'],
                    help='可选 HOYA / OHARA / CDGM / HIKARI')
    ap.add_argument('--pairs', nargs='+', required=True, help='形如 1.48749/70.44')
    ap.add_argument('--json', action='store_true')
    a = ap.parse_args()
    pairs = [tuple(float(x) for x in p.split('/')) for p in a.pairs]
    tol = tolerances([p[0] for p in pairs])
    libs = {v: load(v) for v in a.vendors}
    res = []
    for nd, vd in pairs:
        row = {'nd': nd, 'vd': vd,
               'origin': exact_vendor(nd, vd, list(libs.values()), tol)}
        for v, lib in libs.items():
            top = match(nd, vd, lib, tol)
            row[v] = [{'name': g['name'], 'nd': g['nd'], 'vd': g['vd'],
                       'dnd': round(g['nd'] - nd, 6), 'dvd': round(g['vd'] - vd, 3),
                       'grade': grade(g, nd, vd, tol)} for g in top]
        res.append(row)
    if a.json:
        print(json.dumps(res, ensure_ascii=False, indent=1))
    else:
        print('容差 nd±%.4f / vd±%.2f' % tol)
        for r in res:
            print('%.5f/%.2f  原设计推定: %s' % (r['nd'], r['vd'], r['origin'] or '(无完全一致)'))
            for v in a.vendors:
                g = r[v][0]
                alt = ', '.join('%s(%.5f/%.2f)' % (x['name'], x['nd'], x['vd']) for x in r[v][1:])
                print('   %-6s %-12s %.5f/%-6.2f Δ%+.5f/%+.2f  [%s]   备选: %s'
                      % (v, g['name'], g['nd'], g['vd'], g['dnd'], g['dvd'], g['grade'], alt))


if __name__ == '__main__':
    main()
