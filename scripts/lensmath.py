#!/usr/bin/env python3
"""一次跑完：d/e 线自动判定 → 玻璃匹配 → 近轴校验 → 对焦位置解。

    python3 lensmath.py spec.json --vendors HOYA OHARA CDGM --betas 0.02 0.06 --write

自动做的事:
  1. 先按 nd 匹配；若中位 |Δ| 明显偏大而按 ne(546.074nm) 匹配后 Δ 收敛，判定专利印的是 e 线，
     并把 spec 里的 nd 换成对应牌号的 d 线目录值（印刷值转存 extra）。
  2. ΣD vs 专利 L、近轴 EFL vs 专利 f、各群焦距，全部对一遍。
  3. 由「对焦组两侧间隔之和恒定 + 像面固定在 ∞ 近轴焦点」解出目标倍率的间隔。
--write 会把 glass / 修正后的 nd / zmx.configs 写回 spec.json。
"""
import os, re, sys, json, math, argparse
HERE = os.path.dirname(os.path.abspath(__file__))
GDIR = os.environ.get('PATENT_GLASS_DIR') or os.path.join(HERE, '..', 'assets', 'glass')
MOLD = re.compile(r'^(D-|M-|MP-|MC-|Q-|L-)')
# —— 牌号世代（只在**同一家目录内部**排序，不跨厂家）——
# HIKARI（尼康自家玻璃厂）：E- 是老款、P- 是模压料的旧名、完全没前缀的（SK15/LAF9/PK2）更老，
# 现行是 **J-（研磨）/ Q-（模压）**。用户 2026-09 明确：新尼康镜头一般就是 J- 开头的。
# 目录里 E-/J- 常常是同 nd 的新旧两代（E-LAK01↔J-LAK01、E-BK7↔J-BK7A、E-PSKH1↔J-PSKH1、
# E-SK15↔J-SK15），按 Δ 排序必然挑到 E-，所以必须显式让老牌号靠后。
CURRENT = {'HIKARI': re.compile(r'^[JQ]-')}
def _gen(g):
    r = CURRENT.get(g.get('vendor'))
    return 0 if (r is None or r.match(g['name'])) else 1
LD, LE = 0.5875618, 0.5460740
_BADF = {}   # 色散公式编号 -> 自校验没过的牌号数（加载后会打印告警）
LG, LF, LC = 0.4358343, 0.4861327, 0.6562725

def dpgf(g):
    """Zemax 模型玻璃的 dPgF：相对部分色散对「正常线」的偏离。
    Pg,F=(ng-nF)/(nF-nC)；正常线 Pg,F0 = 0.6438 - 0.001682*vd。
    没有色散公式系数的厂家（CSV）返回 None。"""
    c = g.get('cd')
    if not c: return None
    if not c: return None
    ng, nF, nC = (_n(c, g['f'], L) for L in (LG, LF, LC))
    if None in (ng, nF, nC) or nF == nC: return None
    return round((ng - nF) / (nF - nC) - (0.6438 - 0.001682 * g['vd']), 5)

def _decode(raw):
    if raw[:2] in (b'\xff\xfe', b'\xfe\xff'): return raw.decode('utf-16')
    for e in ('utf-8', 'latin-1'):
        try:
            t = raw.decode(e)
            if t.count('NM ') > 10: return t
        except UnicodeDecodeError: pass
    return raw.decode('utf-16', errors='ignore')

def _n(c, f, lam):
    """AGF 的 CD 行系数 → 折射率。f 是 NM 行第 2 个数（色散公式编号）。

    HIKARI 目录实测分布：formula 1 有 182 支、**13 有 175 支**、12 有 21 支、6 有 12 支。
    老版本把 13 当成 Sellmeier 型（1+ΣKλ²/(λ²−L)）算，那是错的（13 = Extended3），
    12 更是压根没实现 → 直接返回 None，dPgF 写成 0，正犯「模型玻璃 dPgF 不能留 0」那条红线。
    下面 1/2/6/12/13 都用 HIKARI/HOYA/OHARA 的目录 nd 逐一验过（能精确复现印刷值）。
    """
    L2 = lam * lam
    try:
        if f == 1:      # Schott
            n2 = c[0] + c[1]*L2 + c[2]/L2 + c[3]/L2**2 + c[4]/L2**3 + c[5]/L2**4
        elif f == 2:    # Sellmeier1
            n2 = 1 + sum(c[i] * L2 / (L2 - c[i + 1]) for i in (0, 2, 4))
        elif f == 6:    # Sellmeier3
            n2 = 1 + sum(c[i] * L2 / (L2 - c[i + 1]) for i in (0, 2, 4, 6))
        elif f == 11:   # Sellmeier5
            n2 = 1 + sum(c[i] * L2 / (L2 - c[i + 1]) for i in (0, 2, 4, 6, 8))
        elif f == 10:   # Extended
            n2 = (c[0] + c[1]*L2 + c[2]/L2 + c[3]/L2**2 + c[4]/L2**3
                  + c[5]/L2**4 + c[6]/L2**5 + c[7]/L2**6)
        elif f == 12:   # Extended2
            n2 = (c[0] + c[1]*L2 + c[2]/L2 + c[3]/L2**2 + c[4]/L2**3
                  + c[5]/L2**4 + c[6]*L2**2 + c[7]*L2**3)
        elif f == 13:   # Extended3
            n2 = (c[0] + c[1]*L2 + c[2]*L2**2 + c[3]/L2 + c[4]/L2**2
                  + c[5]/L2**3 + c[6]/L2**4 + c[7]/L2**5 + c[8]/L2**6)
        else: return None
        return math.sqrt(n2) if n2 > 0 else None
    except Exception:
        return None

def load(vendor):
    for ext in ('.AGF', '.agf'):
        p = os.path.join(GDIR, vendor.upper() + ext)
        if os.path.exists(p): break
    else:
        p = os.path.join(GDIR, vendor.upper() + '.csv')
        import csv
        out = []
        for r in csv.DictReader(open(p, encoding='utf-8-sig')):
            out.append({'name': r['name'].strip(), 'nd': float(r['nd']), 'vd': float(r['vd']),
                        'ne': None, 'vendor': vendor.upper()})
        return out
    out, cur = [], None
    for line in _decode(open(p, 'rb').read()).splitlines():
        if line.startswith('NM '):
            q = line.split()
            cur = {'name': q[1], 'f': int(float(q[2])), 'nd': float(q[4]), 'vd': float(q[5]),
                   'vendor': vendor.upper(), 'cd': None, 'ne': None}
            out.append(cur)
        elif line.startswith('CD ') and cur is not None:
            cur['cd'] = [float(x) for x in line.split()[1:]]
            # 自校验：色散公式没实现或实现错了，反算的 d 线折射率不会等于 NM 行的 nd。
            # 对不上就把系数作废（ne=None、dPgF 跳过），绝不拿错曲线往下算。
            chk = _n(cur['cd'], cur['f'], LD)
            if chk is None or abs(chk - cur['nd']) > 2e-4:
                _BADF[cur['f']] = _BADF.get(cur['f'], 0) + 1
                cur['cd'] = None
            else:
                cur['ne'] = _n(cur['cd'], cur['f'], LE)
    return out

def _ulp(v):
    """按印刷位数给出半个最小刻度，作为「印刷精度内即平局」的判据。"""
    t = ('%.8f' % v).rstrip('0')
    d = len(t.split('.')[1]) if '.' in t else 0
    return 0.5 * 10 ** (-d)

BRAND_SLACK = 0.0020      # Δnd 落在最优值 + 这个带宽内，视为「同一档」，改按厂家优先级选
BRAND_SLACK_VD = 1.0

def best(val, vd, libs, key, asph, brand_slack=BRAND_SLACK, match_tol=2e-5, alt_only=()):
    un, uv = _ulp(val), _ulp(vd)
    c = []
    for rank, lib in enumerate(libs.values()):
        for g in lib:
            gv = g[key]
            if gv is None: continue
            s = ((gv - val) / 2e-4) ** 2 + ((g['vd'] - vd) / 0.15) ** 2
            if MOLD.match(g['name']) and not asph: s += 60
            # 边界要含进来：印 νd=30.1 时半刻度正好 0.05，目录里 30.05 的牌号
            # 会因为浮点误差被判出局，把首选莫名其妙让给下一家。给 1e-9 的容差。
            tie = 0 if (abs(gv - val) <= un * (1 + 1e-9) and abs(g['vd'] - vd) <= uv * (1 + 1e-9)
                        and not (MOLD.match(g['name']) and not asph)) else 1
            c.append((s, abs(gv - val), g, rank, tie))
    # 印刷精度以内一律视为平局，按 --vendors 给出的厂家顺序优先
    c.sort(key=lambda x: (x[4], x[3] if x[4] == 0 else 0, x[0]))
    # ================= 厂家优先 = 硬约束（不是打分项）=================
    # 用户反复强调过：索尼/尼康/佳能/腾龙/适马等日系厂商基本不用国产玻璃。
    # 旧逻辑「精确命中优先于厂家顺序」在 νd 第 4 位小数上崩掉了：
    #   WO2019187633 面9 印 1.58313/59.38 —— CDGM D-ZK2 的 νd=59.3817（Δ0.0017，
    #   落在 2 位小数的半刻度 0.005 内 → tie=0），HOYA M-BACD12 νd=59.4600
    #   （Δ0.08 → tie=1），于是唯一的 tie=0 是 CDGM，厂家优先带又因为「已有精确命中」
    #   被跳过，结果给索尼专利选了国产料。
    # 各家 νd 的定义/舍入本来就有 0.1 量级的差异，用它当厂家裁决依据毫无道理。
    # 改成：按 --vendors 顺序找**第一家有可接受候选**的厂商，在这家内部再按分数选。
    #   可接受 = |Δnd| ≤ max(印刷半刻度, match_tol) 且 |Δvd| ≤ max(印刷半刻度, VD_BAND)
    # alt_only 里的厂家永远只进备选列，不做首选（日系专利传 --alt-only CDGM）。
    nd_band = max(un, match_tol) * (1 + 1e-9)
    vd_band_ = max(uv, VD_BAND) * (1 + 1e-9)
    names = list(libs.keys())
    for r in range(len(names)):
        if names[r] in (alt_only or ()): continue
        pool = [x for x in c if x[3] == r
                and abs(x[2][key] - val) <= nd_band
                and abs(x[2]['vd'] - vd) <= vd_band_
                and (asph or not MOLD.match(x[2]['name']))]
        if pool:
            pool.sort(key=lambda x: (_gen(x[2]), x[0]))
            if pool[0] is not c[0]:
                b = c[0]
                c = [pool[0]] + [x for x in c if x is not pool[0]]
                c[0] = c[0] + (b,)
            return c
        if any(x[3] == r and abs(x[2][key] - val) <= nd_band for x in c):
            # 这家在 nd 上够近但 νd 差太多 / 只有模压料配球面 —— 让给下一家，但别静默
            pass
    # 品牌优先带：日系专利用日系料。Δnd 差在 brand_slack 以内算「同一档」，
    # 这一档里改按 --vendors 顺序选，免得第 4~5 位小数把首选判给了厂家根本不会用的牌号。
    # 品牌带只在**没有任何牌号落在印刷精度内**时才启用。
    # 若已有精确命中（tie==0），说明专利印的就是那个牌号的目录值，
    # 这时再按厂家顺序挪走会把确凿的证据推翻 —— 精确命中之间的排序交给 tie 逻辑。
    if c and brand_slack > 0 and c[0][4] != 0:
        b = c[0]
        okb = [x for x in c
               if x[1] <= b[1] + brand_slack
               and abs(x[2]['vd'] - vd) <= abs(b[2]['vd'] - vd) + BRAND_SLACK_VD
               and (asph or not MOLD.match(x[2]['name']) or MOLD.match(b[2]['name']))]
        if okb:
            okb.sort(key=lambda x: (x[3], x[0]))
            if okb[0] is not b:
                c = [okb[0]] + [x for x in c if x is not okb[0]]
                c[0] = c[0] + (b,)          # 记下「本来更近的那个」，报告里说明
    # 兜底红线：谁都不「可接受」时上面两段都不生效，首选可能落到 alt_only 厂家头上
    # （实测 JP2019-144441 面25 盖板 1.5230/58.59 → CDGM H-K51）。日系专利里出现
    # 国产牌号做首选是用户明确打回过的，这里无条件换成最近的非 alt_only 牌号。
    if c and alt_only and c[0][2]['vendor'] in alt_only:
        na = [x for x in c if x[2]['vendor'] not in alt_only]
        if na:
            b = c[0]; na.sort(key=lambda x: x[0])
            c = [na[0]] + [x for x in c if x is not na[0]]
            c[0] = c[0] + (b,)
    return c

# 各厂 νd 的定义与舍入差异本来就有 0.1 量级，不能拿它当厂家裁决依据
VD_BAND = 0.12

def surfaces_with_glass(emb):
    return [s for s in emb['surfaces'] if s.get('nd')]

def detect_line(emb, libs):
    """返回 ('nd'|'ne', 中位Δ_nd, 中位Δ_ne)。"""
    med = {}
    for key in ('nd', 'ne'):
        ds = []
        for s in surfaces_with_glass(emb):
            c = best(s['nd'], s['vd'], libs, key, s.get('type') == '非球面')
            ds.append(c[0][1] if c else 9)
        ds.sort(); med[key] = ds[len(ds) // 2]
    return ('ne' if med['ne'] < med['nd'] * 0.2 else 'nd'), med['nd'], med['ne']

def paraxial(surfs, dmap=None, obj=None):
    """surfs: [{'R','D','nd'}...] 已展开可变间隔。返回 (EFL, BF, beta)。"""
    n, y = 1.0, 1.0
    u = (1.0 / obj) if obj else 0.0
    u1 = u
    last = len(surfs) - 1
    for k, s in enumerate(surfs):
        c = 0.0 if s['R'] in (None, 0) else 1.0 / s['R']
        n2 = s['nd'] or 1.0
        u = (n * u - y * c * (n2 - n)) / n2
        if k < last: y += u * s['D']
        n = n2
    return (None if obj else -1.0 / u), -y / u, (u1 / u if obj else None)

def brentq(f, lo, hi, xtol=1e-12, maxiter=200):
    """二分法，替代 scipy.optimize.brentq —— 本脚本只需要单调区间求根。"""
    flo, fhi = f(lo), f(hi)
    if flo * fhi > 0:
        raise ValueError('f(a) 与 f(b) 同号，无法括号求根')
    for _ in range(maxiter):
        mid = 0.5 * (lo + hi)
        fm = f(mid)
        if fm == 0 or (hi - lo) < xtol: return mid
        if flo * fm < 0: hi, fhi = mid, fm
        else: lo, flo = mid, fm
    return 0.5 * (lo + hi)


def solve_obj(g, f, ok=None):
    """解「像面固定」的物距 o：g(o)=0。解不出返回 None（＝物在无穷远之外）。

    旧写法固定拿括号 [1.02|f|, 1e4…1e9]，当 1.02|f| 落在**前焦点以内**时
    g(lo) 与 g(hi) 同号（中间跨了一个极点：物正好在前焦点上，像跑到无穷远），
    于是整条行程一个解都找不到 —— WO2024147268 的双浮动分支实测全军覆没
    （每个 D12 都报「超出对焦行程」）。

    现在：对数网格扫描收集**全部**变号区间，从**最大 o** 往小依次求根，
    第一个通过 ok() 体检的就是答案。
    ok 必须卡住倍率：极点两侧 g 也会真的穿过 0，那个根物理上是废的
    （实测 D12=1.749 解出 o=205.58、β=-2.6e9），只按「最大 o」挑仍会中招 ——
    因为 D12=1.749 本来就该无解（物在无穷远），此时极点根成了唯一候选。"""
    import math
    INF = float('inf')
    def gs(o):
        # 物正好落在前焦点上时 u=0（像跑到无穷远），paraxial 会 ZeroDivisionError。
        try: v = g(o)
        except Exception: return float('nan')
        return v if (v == v and abs(v) != INF) else float('nan')
    lo0 = 1.02 * abs(f)
    a0, b0, N = math.log(lo0), math.log(1e9), 240
    xs = [math.exp(a0 + (b0 - a0) * i / N) for i in range(N + 1)]
    vs = [gs(o) for o in xs]
    brs = []
    for i in range(N):
        v1, v2 = vs[i], vs[i + 1]
        if v1 != v1 or v2 != v2: continue
        if v1 == 0.0: brs.append((xs[i], xs[i])); continue
        if v1 * v2 < 0: brs.append((xs[i], xs[i + 1]))
    def bisect(lo, hi, flo):
        if lo == hi: return lo
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            fm = gs(mid)
            if fm != fm:                       # 正好踩到极点，微移一点
                mid = lo + 0.499 * (hi - lo); fm = gs(mid)
                if fm != fm: return None
            if fm == 0.0 or (hi - lo) < 1e-9 * max(1.0, abs(hi)): return mid
            if flo * fm < 0: hi = mid
            else: lo, flo = mid, fm
        return 0.5 * (lo + hi)
    for lo, hi in reversed(brs):               # 从最大 o 往小试
        r = bisect(lo, hi, gs(lo))
        if r is None: continue
        if ok is None or ok(r): return r
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('spec'); ap.add_argument('--vendors', nargs='+', default=['HOYA', 'OHARA', 'CDGM'])
    ap.add_argument('--brand-slack', type=float, default=BRAND_SLACK,
                    help='Δnd 落在最优值+该带宽内即按厂家优先级选（默认 0.0020，设 0 关闭）')
    ap.add_argument('--betas', nargs='*', type=float, default=[0.02, 0.06])
    ap.add_argument('--match-tol', type=float, default=2e-5,
                    help='d 线时认作「等效牌号」的 |Δnd| 上限（默认 2e-5）；超出则保留印刷 nd/vd 且不写 glass')
    ap.add_argument('--alt-only', nargs='*', default=[],
                    help='这些厂家只进备选列、永不做首选。日系申请人固定传 --alt-only CDGM')
    ap.add_argument('--mfd', type=float, default=None,
                    help='产品标称最短撮影距離(mm，物体→像面)。给了就额外解一个 MFD 状态：'
                         '物距 = MFD − ΣD。专利近距态离产品 MFD 很远时用它。')
    ap.add_argument('--emb', type=int, default=0)
    ap.add_argument('--write', nargs='?', const='AUTO', default=None,
                    help='写出匹配后的 spec（默认 <spec>.matched.json，不覆盖原文件）')
    a = ap.parse_args()
    spec = json.load(open(a.spec, encoding='utf-8'))
    emb = spec['embodiments'][a.emb]
    libs = {v: load(v) for v in a.vendors}
    if _BADF:
        print('  ⚠ 色散公式自校验未通过（这些牌号不参与 e 线判定与 dPgF）: '
              + ', '.join('formula %d × %d 支' % kv for kv in sorted(_BADF.items())))
    line, mnd, mne = detect_line(emb, libs)
    print('== 折射率列判定 ==')
    print('  按 nd 匹配中位偏差 %.5f ；按 ne 匹配中位偏差 %.5f  →  判定为 %s 线'
          % (mnd, mne, 'e' if line == 'ne' else 'd'))
    if line == 'ne':
        print('  （专利「ndi」列实为 ne(546.074nm)；下表 nd 已换成对应牌号 d 线目录值）')
    NOMATCH = []
    print('\n== 玻璃匹配 ==')
    ASPH_LENS = {q.get('lens') for q in emb['surfaces']
                 if q.get('type') == '非球面' and q.get('lens')}
    for s in surfaces_with_glass(emb):
        asph = s.get('type') == '非球面'
        asph_eff = asph or (s.get('lens') in ASPH_LENS)
        c = best(s['nd'], s['vd'], libs, line, asph, a.brand_slack, a.match_tol, tuple(a.alt_only))
        g = c[0][2]
        promoted = c[0][5] if len(c[0]) > 5 else None
        alts = []
        seen = {g['vendor']}
        for _, _, x, _r, _t in c[1:]:
            if x['vendor'] not in seen:
                seen.add(x['vendor']); alts.append('%s %s(Δ%+.5f)' % (x['vendor'], x['name'], x[line] - s['nd']))
            if len(alts) == 2: break
        tag = ''
        if promoted is not None:
            pg = promoted[2]
            tag = '   ← 品牌优先（%s %s 更近 Δ%+.5f，但按 --vendors 顺序不选它）' % (
                pg['vendor'], pg['name'], pg[line] - s['nd'])
        print('  面%-4s %.5f/%-5.1f %s → %-6s %-13s nd=%.5f νd=%.2f Δ%+.5f | %s%s'
              % (s['i'], s['nd'], s['vd'], 'ASP' if asph else '   ', g['vendor'], g['name'],
                 g['nd'], g['vd'], g[line] - s['nd'], ' ; '.join(alts), tag))
        # 无等效牌号的面要另挑一颗「贴合实物」的基准玻璃做 Offset 解：
        #   现行牌号优先；**非球面元件再优先取模压料**（Q-/M-/D-/L- …）——
        #   非球面落在模压料上本来就是匹配可信的标志，也正是实物用的料。
        base = g
        if abs(g[line] - s['nd']) > a.match_tol:
            # 候选必须同时靠近 nd **和** νd —— 只卡 nd 会挑出阿贝数差 30 的火石料
            # （实测面23 一度选到 J-F16：nd 只差 0.00024，νd 却差 32.6）。
            same = [x[2] for x in c if x[2]['vendor'] == g['vendor']
                    and abs(x[2][line] - s['nd']) <= 0.002
                    and abs(x[2]['vd'] - s['vd']) <= 1.0 and _gen(x[2]) == 0]
            # 非球面按**整块元件**判：胶合/单片里只要有一面是非球面，这块玻璃就是
            # 模压件，nd 挂在前表面而 * 号常常打在后表面（本篇 L31 = 面23+面24）。
            if asph_eff:
                mold = [x for x in same if MOLD.match(x['name'])]
                if mold:
                    mold.sort(key=lambda x: abs(x[line] - s['nd'])); base = mold[0]
            if base is g and same:
                same.sort(key=lambda x: abs(x[line] - s['nd'])); base = same[0]
        dn = g[line] - s['nd']
        if a.write:
            s.setdefault('extra', {})
            pg = dpgf(g)
            if pg is not None: s['dpgf'] = pg
            if line == 'ne':
                s['extra']['印刷 ndi (实为 e 线 ne)'] = s['nd']
                s['glass'] = '%s %s' % (g['vendor'], g['name'])
                s['nd'] = round(g['nd'], 5); s['vd'] = g['vd']
            elif abs(dn) <= a.match_tol:
                # 印刷精度以内 = 专利印的就是这个牌号的目录值，换成目录值无害
                s['glass'] = '%s %s' % (g['vendor'], g['name'])
                s['nd'] = round(g['nd'], 5); s['vd'] = g['vd']
            else:
                # d 线且无等效牌号：印刷值才是权威，绝不能被目录值顶掉（会改变 EFL）
                s['extra']['最接近牌号(非等效)'] = '%s %s Δnd%+.5f' % (g['vendor'], g['name'], dn)
                # Offset 玻璃解（Zemax 官方，solve code 4）：基准目录玻璃 + Nd/Vd 偏移。
                # 这样既精确还原专利印刷的 nd/vd，又保住真实玻璃的整条色散曲线 ——
                # 比模型玻璃好得多，是用户 2026-09 点名要的做法。
                s['glass_offset'] = {
                    'base': '%s %s' % (base['vendor'], base['name']),
                    'base_nd': base['nd'], 'base_vd': base['vd'],
                    'base_dpgf': dpgf(base), 'd_nd': round(s['nd'] - base['nd'], 6),
                    'd_vd': round(s['vd'] - base['vd'], 4)}
                pb = dpgf(base)
                if pb is not None: s['dpgf'] = pb
                nt = ('★ 四家目录无等效牌号（最近 %s %s Δnd%+.5f）；保留专利印刷值，'
                      'zmx 用 Offset 玻璃解：基准 %s + Nd%+.5f / Vd%+.4f'
                      % (g['vendor'], g['name'], dn, s['glass_offset']['base'],
                         s['glass_offset']['d_nd'], s['glass_offset']['d_vd']))
                s['note'] = (s.get('note', '') + '；' + nt).lstrip('；')
                NOMATCH.append((s['i'], s['nd'], s['vd'], g['vendor'], g['name'], dn,
                                s['glass_offset']))
    if NOMATCH:
        print('\n  ★ 无等效牌号（保留专利印刷值，zmx 该面走模型玻璃）:')
        for i, nd, vd, v, n, dn, go in NOMATCH:
            print('     面%-4s %.5f/%.1f  最近 %s %s Δnd%+.5f' % (i, nd, vd, v, n, dn))
            print('            → Offset 解：基准 %-18s nd=%.5f νd=%.4f  Nd offset %+.5f  Vd offset %+.4f'
                  % (go['base'], go['base_nd'], go['base_vd'], go['d_nd'], go['d_vd']))
    # ---- 近轴 ----
    zx = spec.get('zmx', {}); fc = zx.get('focus', {})
    st = emb.get('states', [None])[0]
    def expand(state, dmapover=None):
        out = []
        for s in emb['surfaces']:
            if s['i'] in ('IMG',): continue
            D = s['D']
            if isinstance(D, str):
                D = (dmapover or {}).get(D) or emb['variable'][D][state]
            out.append({'R': s['R'], 'D': float(D), 'nd': s.get('nd')})
        return out
    base = expand(st)
    SD = sum(s['D'] for s in base)
    f, bf, _ = paraxial(base)
    gen = dict((k, v) for k, v in emb.get('general', []))
    print('\n== 近轴校验 ==')
    print('  ΣD = %.2f   (专利 L = %s)' % (SD, gen.get('L 光学全长 (mm)', '?')))
    print('  EFL = %.4f  BF = %.4f   (专利 f = %s)' % (f, bf, gen.get('f (mm)', '?')))
    fc2 = zx.get('focus2')
    if not fc: print('\n(spec 无 zmx.focus，跳过对焦解)'); 
    elif fc2:
        # ===== 双浮动对焦群（两组独立移动，如索尼 135GM）=====
        # 一个像面共轭方程解不出两个未知量，所以把第 2 组的位置当成第 1 组的函数：
        # 用专利给出的各对焦状态拟合凸轮曲线 kb2 = q(kb1)（n 点 → n-1 次拉格朗日），
        # 再对每个 kb1 解物距使像面不动 → 得到 β。这是专利信息量的上限，
        # 中间态的两组分工是插值出来的，回话必须说明。
        kb1, ka1, tot1 = fc['key_before'], fc.get('key_after'), fc['sum']
        kb2, ka2, tot2 = fc2['key_before'], fc2.get('key_after'), fc2['sum']
        # 链式（三段浮动）：两个对焦群串在一起，只有 kb1+kb2+klast = tot1 守恒，
        # 两两配对都不守恒（适马 24 Art / JP2016-12034 就是这种）。
        klast = fc.get('key_last')
        sts = emb['states']
        xs = [float(emb['variable'][kb1][q]) for q in sts]
        ys = [float(emb['variable'][kb2][q]) for q in sts]
        def cam(x):
            r = 0.0
            for i in range(len(xs)):
                t = ys[i]
                for j in range(len(xs)):
                    if i != j: t *= (x - xs[j]) / (xs[i] - xs[j])
                r += t
            return r
        def dmap_of(x):
            y = cam(x)
            if klast:
                return {kb1: x, kb2: y, klast: tot1 - x - y}
            # key_after 为 null = 「整组前伸式对焦」：两个可变间隔各自独立，
            # 前面没有固定端可配对守恒（如适马 70 Art 这类全体繰り出し微距）。
            d = {kb1: x, kb2: y}
            if ka1: d[ka1] = tot1 - x
            if ka2: d[ka2] = tot2 - y
            return d
        def show(x):
            d = dmap_of(x)
            return '  '.join('%s=%.4f' % (k, v) for k, v in d.items())
        def bfd2(x, obj=None):
            return paraxial(expand(st, dmap_of(x)), obj=obj)
        # 链式对焦时 BF 本身是变量，像面固定的条件是「近轴后焦 = 该状态的 BF」，
        # 而不是「后焦恒等于 ∞ 态的后焦」。BFoff 吸收专利 BF 的印刷舍入。
        if klast:
            # 像面在空间中固定 ⇔「var_before 面到最后一面顶点的累计厚度 + 近轴后焦」恒定。
            # 直接拿 BF 变量当后焦是错的：补了盖板以后最后一面已经不是 BF 那一面。
            def zpos(x):
                ss = expand(st, dmap_of(x))
                return sum(q['D'] for q in ss[fc['var_before']-1:-1])
            CTOT = zpos(xs[0]) + bfd2(xs[0])[1]
            target = lambda x: CTOT - zpos(x)
        else:
            BF0 = bfd2(xs[0])[1]
            target = lambda x: BF0
        def solve_d0(x):
            # 只认「物在前焦点之外、成实像」的那一支：极点另一侧的根 |β| 会大到 1e9
            def ok(o):
                try: b = bfd2(x, o)[2]
                except Exception: return False
                return b is not None and -5.0 < b < 0.0
            return solve_obj(lambda o: bfd2(x, o)[1] - target(x), f, ok)
        def beta2(x):
            o = solve_d0(x)
            if o is None: return 0.0, None
            return bfd2(x, o)[2], o
        # 扫描范围按**机械行程**给，不是「专利记载的 ∞↔近距」那一段（与单组分支一致）。
        # 双浮动时还要保证凸轮映射出来的第 2 组间隔也留得住 0.10mm，否则解出的是非物理位置。
        p_lo, p_hi = min(xs), max(xs)
        d0v = xs[0]
        def ok_x(x):
            y = cam(x)
            if not (0.10 <= x <= tot1 - 0.10): return False
            if klast:
                return 0.10 <= y and 0.10 <= tot1 - x - y
            return 0.10 <= y <= tot2 - 0.10
        lo_d, hi_d = (d0v, tot1 - 0.10) if xs[-1] > d0v else (0.10, d0v)
        N = 800
        tab = []
        for i in range(N + 1):
            x = lo_d + (hi_d - lo_d) * i / N
            if not ok_x(x): continue
            b, o = beta2(x)
            if o: tab.append((x, abs(b)))
        print('\n== 对焦位置解（双浮动对焦群）==')
        if klast:
            print('  链式三段: %s+%s+%s = %.4f 恒定（G1 与像面固定，两组对焦透镜各自移动）' % (kb1, kb2, klast, tot1))
            print('  对焦群1 起于面%d，对焦群2 起于面%d' % (fc['var_before']+1, fc2['var_before']+1))
        elif ka1 and ka2:
            print('  对焦群1: 面%d〜面%d, %s+%s = %.4f 恒定' % (fc['var_before']+1, fc['var_after'], kb1, ka1, tot1))
            print('  对焦群2: 面%d〜面%d, %s+%s = %.4f 恒定' % (fc2['var_before']+1, fc2['var_after'], kb2, ka2, tot2))
        else:
            print('  全体繰り出し式：%s（面%d 后）与 %s（面%d 后）各自独立，无配对守恒；'
                  '像面固定条件 = 近轴后焦恒定' % (kb1, fc['var_before'], kb2, fc2['var_before']))
        print('  凸轮曲线 %s = q(%s) 由专利 %d 个状态拟合' % (kb2, kb1, len(sts)))
        print('  %-14s d0=∞           %s  β=0' % (sts[0], show(xs[0])))
        cfgs = [{'name': sts[0], 'd0': 'INFINITY', kb1: xs[0], kb2: ys[0]}]
        def x_for(m):
            pair = [t for t in zip(tab, tab[1:]) if (t[0][1] - m) * (t[1][1] - m) <= 0]
            if not pair: return None
            (x1, b1), (x2, b2) = pair[0]
            for _ in range(60):
                xm = 0.5 * (x1 + x2); bm = abs(beta2(xm)[0])
                if (b1 - m) * (bm - m) <= 0: x2, b2 = xm, bm
                else: x1, b1 = xm, bm
            return 0.5 * (x1 + x2)
        for m in a.betas:
            x = x_for(m)
            if x is None:
                print('  %.2fx: 该倍率超出对焦行程' % m); continue
            b, o = beta2(x); y = cam(x)
            ext = '' if p_lo - 1e-9 <= x <= p_hi + 1e-9 else '   ← ★外推：超出专利记载的对焦范围'
            print('  %-14s d0=%-10.2f  %s  β=%+.5f%s' % ('%.2fx' % m, o, show(x), b, ext))
            cfgs.append({'name': '%.2fx' % m, 'd0': round(o, 2),
                         kb1: round(x, 4), kb2: round(y, 4),
                         'extrapolated': bool(ext)})
        if a.mfd:
            tgt = a.mfd - SD
            xs_scan = [t[0] for t in tab]
            o_of = lambda x: (beta2(x)[1] or float('inf'))
            pr = [p for p in zip(xs_scan, xs_scan[1:])
                  if (o_of(p[0]) - tgt) * (o_of(p[1]) - tgt) <= 0]
            if not pr:
                print('  MFD(%gmm): 超出对焦行程' % a.mfd)
            else:
                x1, x2 = pr[0]
                for _ in range(60):
                    xm = 0.5 * (x1 + x2)
                    if (o_of(x1) - tgt) * (o_of(xm) - tgt) <= 0: x2 = xm
                    else: x1 = xm
                x = 0.5 * (x1 + x2); b, o = beta2(x); y = cam(x)
                ext = '' if p_lo - 1e-9 <= x <= p_hi + 1e-9 else '   ← ★外推：超出专利记载的对焦范围'
                nm = 'MFD(%gmm)' % a.mfd
                print('  %-14s d0=%-10.2f  %s  β=%+.5f  (1:%.1f)%s'
                      % (nm, o, show(x), b, 1 / abs(b) if b else 0, ext))
                cfgs.append({'name': nm, 'd0': round(o, 2), kb1: round(x, 4),
                             kb2: round(y, 4), 'extrapolated': bool(ext)})
        for extra in sts[1:]:
            x = float(emb['variable'][kb1][extra]); y = float(emb['variable'][kb2][extra])
            b, o = beta2(x)
            given = (emb.get('d0') or {}).get(extra)
            tagg = ''
            if given:
                tagg = '  (专利给定 d0=%.4f，解出 %.2f)' % (float(given), o or float('nan'))
                o = float(given)
            print('  %-14s d0=%-10.2f  %s  β=%+.5f%s'
                  % (extra, o or float('nan'), show(x), b, tagg))
            if extra == sts[-1]:
                cfgs.append({'name': extra, 'd0': round(o, 4) if o else None,
                             kb1: x, kb2: y})
        if a.write: spec.setdefault('zmx', {})['configs'] = cfgs
    else:
        vb, va, tot = fc['var_before'], fc['var_after'], fc['sum']
        kb, ka = fc['key_before'], fc['key_after']
        def bfd(d13, obj=None):
            s = expand(st, {kb: d13, ka: tot - d13})
            return paraxial(s, obj=obj)
        BF0 = bfd(emb['variable'][kb][st])[1]

        def solve_d0(d13):
            """像面固定在 ∞ 近轴焦点，解物距；解不出返回 None（等于物在无穷远之外）。"""
            def ok(o):
                try: b = bfd(d13, o)[2]
                except Exception: return False
                return b is not None and -5.0 < b < 0.0
            return solve_obj(lambda o: bfd(d13, o)[1] - BF0, f, ok)

        def beta(d13):
            o = solve_d0(d13)
            if o is None: return 0.0, None
            return bfd(d13, o)[2], o

        d0v = emb['variable'][kb][st]
        near = emb['variable'][kb][emb['states'][-1]] if len(emb.get('states', [])) > 1 else 0.2 * d0v
        # 扫描范围要按**机械行程**给，不能只给「专利记载的 ∞↔近距」那一段 ——
        # 有的专利（如 US20150092100 的近距态只到 0.033x）记载范围比镜头实际行程短得多，
        # 只扫记载范围会把 0.06x 误判成「超出对焦行程」。
        # 方向仍锁在对焦组实际移动的那一侧，避免扫到另一支非物理的解。
        p_lo, p_hi = (near, d0v) if near < d0v else (d0v, near)
        lo_d, hi_d = (d0v, tot - 0.10) if near > d0v else (0.10, d0v)
        N = 400
        tab = []
        for i in range(N + 1):
            x = lo_d + (hi_d - lo_d) * i / N
            b, o = beta(x)
            if o: tab.append((x, abs(b)))
        tab.sort(key=lambda t: t[1])
        print('\n== 对焦位置解 (对焦组 面%d〜面%d, %s+%s = %.4f 恒定) ==' % (vb + 1, va, kb, ka, tot))
        print('  %-12s d0=∞          %s=%.4f %s=%.4f  β=0' % (st, kb, d0v, ka, tot - d0v))
        cfgs = [{'name': st, 'd0': 'INFINITY', kb: d0v}]

        def d13_for(m):
            lo = max((x for x, b in tab if b <= m), key=lambda x: abs(x - d0v), default=None)
            pair = [t for t in zip(tab, tab[1:]) if (t[0][1] - m) * (t[1][1] - m) <= 0]
            if not pair: return None
            (x1, b1), (x2, b2) = pair[0]
            for _ in range(60):
                xm = 0.5 * (x1 + x2); bm = abs(beta(xm)[0])
                if (b1 - m) * (bm - m) <= 0: x2, b2 = xm, bm
                else: x1, b1 = xm, bm
            return 0.5 * (x1 + x2)

        for m in a.betas:
            x = d13_for(m)
            if x is None:
                print('  %.2fx: 该倍率超出对焦行程' % m); continue
            b, o = beta(x)
            ext = '' if p_lo - 1e-9 <= x <= p_hi + 1e-9 else '   ← ★外推：超出专利记载的对焦范围'
            print('  %-12s d0=%-10.2f %s=%.4f %s=%.4f  β=%+.5f%s'
                  % ('%.2fx' % m, o, kb, x, ka, tot - x, b, ext))
            cfgs.append({'name': '%.2fx' % m, 'd0': round(o, 2), kb: round(x, 4),
                         'extrapolated': bool(ext)})
        done = {c0['name'] for c0 in cfgs}
        for extra in emb.get('states', [])[1:]:
            if extra in done: continue
            x = emb['variable'][kb][extra]
            b, o = beta(x)
            given = (emb.get('d0') or {}).get(extra)   # 专利自己给了 d0 就用专利的
            if given: o = float(given)
            print('  %-12s d0=%-10.2f %s=%.4f %s=%.4f  β=%+.5f  (专利给定)'
                  % (extra, o or float('nan'), kb, x, ka, tot - x, b))
            cfgs.append({'name': extra, 'd0': round(o, 2) if o else None, kb: x})
        if a.write: spec.setdefault('zmx', {})['configs'] = cfgs
    if a.write:
        out = a.spec[:-5] + '.matched.json' if a.write == 'AUTO' else a.write
        json.dump(spec, open(out, 'w', encoding='utf-8'), ensure_ascii=False, indent=1)
        print('\n已写出 %s（原 spec 未改动）' % out)

if __name__ == '__main__':
    main()
