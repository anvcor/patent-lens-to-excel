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
# OHARA：**S- / L- 开头的才是 2000 年环保化以后的无铅(且无砷)牌号**。
# 目录里 210 个不带前缀的（PBM/PBH/BPH/BAL/BAM/BSL/BSM/TIM/TIH/LAL/LAH/NSL/FSL…）
# 是含铅老系列，只有老镜头（2000 年前）才会用。注意**光靠 AGF 的 Obsolete 位挡不住**：
# PBM2Y 在目录里仍是 Preferred，但它是含铅的 S-TIM2 的老款。用户 2026-09 打回过。
CURRENT = {'HIKARI': re.compile(r'^[JQ]-'), 'OHARA': re.compile(r'^(S-|L-)')}
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
                        'ne': None, 'vendor': vendor.upper(), 'status': 0})
        return out
    out, cur = [], None
    for line in _decode(open(p, 'rb').read()).splitlines():
        if line.startswith('NM '):
            q = line.split()
            cur = {'name': q[1], 'f': int(float(q[2])), 'nd': float(q[4]), 'vd': float(q[5]),
                   'vendor': vendor.upper(), 'cd': None, 'ne': None,
                   # NM 行第 8 个数 = status：0=Standard 1=Preferred 2=Obsolete 3=Special 4=Melt
                   'status': int(float(q[7])) if len(q) > 7 else 0}
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

def legacy(g):
    """停产牌号（AGF status=2）或该厂的老一代（含铅系列）—— 2000 年环保化以后的镜头不会用。"""
    return g.get('status') == 2 or _gen(g) != 0

def best(val, vd, libs, key, asph, brand_slack=BRAND_SLACK, match_tol=2e-5, alt_only=(), eco=True,
         mold_ok=True):
    """asph 是**本面**是否非球面（决定模压料的罚分与平局资格，沿用旧行为）；
    mold_ok 按**整块元件**给（元件任一面是非球面即 True），False 时模压料直接不进候选。
    默认 True 只为 detect_line 这类统计用途保持旧行为；逐面选料必须显式传。"""
    un, uv = _ulp(val), _ulp(vd)
    c = []
    for rank, lib in enumerate(libs.values()):
        for g in lib:
            gv = g[key]
            if gv is None: continue
            if eco and legacy(g): continue      # 含铅/停产的一律不参与（--allow-legacy 可关）
            # 球面元件绝不用模压料（M-/MP-/MC-/D-/Q-/L-）。原先只罚 +60 分，
            # 可 nd 一项动辄上千分，罚分形同虚设：TS-E24 II 面19（胶合球面）选出了 HOYA MP-NBFD10-20
            if MOLD.match(g['name']) and not mold_ok: continue
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

# ================= Offset 玻璃解的基准 =================
# Offset 解能把 nd 的差整块吃掉，吃不掉的是**色散曲线的形状**（νd、部分色散）。
# 所以基准要先贴 νd / dPgF，nd 只是次要项。旧逻辑的三个坑（US20100208366A1 TS-E24 II 实测）：
#   面7  1.61601/58.7 → HOYA BACD4：只在「最近牌号」那一家里找，OHARA S-BSM4（νd 58.72）根本没进候选；
#   面18 1.55400/52.2 → HOYA E-FEL1（νd 45.8，Vd offset +6.4）：窗口里没人就直接拿 nd 最近的；
#   面19 1.84175/37.2 → HOYA MP-NBFD10-20：球面胶合件用了模压料。
OFS_ND_WIN, OFS_VD_WIN = 0.002, 1.0        # 「窗口内」= 基本就是这颗料
OFS_ND_CAP = 0.03                          # 窗口外按色散挑时 nd 最多差这么多（没有就放开）
OFS_SIG_VD, OFS_SIG_PGF, OFS_SIG_ND = 0.2, 0.001, 0.02   # 打分：νd 0.2 ≈ dPgF 0.001 ≈ nd 0.02
OFS_OBSOLETE = 1.0                         # 停产（但无铅）牌号的罚分 ≈ νd 多差 0.2
OFS_VENDOR_VD = 0.5                        # 首选厂家 |Δvd| 不比全局最优差过这么多，也算「合理」
OFS_VENDOR_PGF = 0.002                     # 同上，dPgF（专利印了 θgF 时）

def offset_base(val, vd, libs, key, asph, alt_only=(), eco=True, pgf=None):
    """给无等效牌号的面挑 Offset 解的基准玻璃，返回 (glass, 说明)。
    asph   —— 按**整块元件**判；False 时模压料一律不要。
    pgf    —— 专利若印了 θgF(Pg,F)，传进来参与打分；没印就只看 νd。
    eco    —— 排除含铅老一代（OHARA 非 S-/L-、HIKARI 非 J-/Q-）。**停产但无铅**的
              （OHARA S-BSM4 这类）允许做基准、只小罚：老专利用的正是它们，色散数据仍是真的。"""
    tpg = None if pgf is None else pgf - (0.6438 - 0.001682 * vd)
    vendors = [v for v in libs if v not in alt_only] or list(libs)
    cand = []
    for v in vendors:
        for g in libs[v]:
            if g[key] is None: continue
            if eco and _gen(g): continue
            # CDGM/HOYA 没有世代前缀可判，停产的（CDGM F/ZF…）多半含铅 → 只放行有前缀可证无铅的厂家
            if eco and g.get('status') == 2 and v not in CURRENT: continue
            mold = bool(MOLD.match(g['name']))
            if mold and not asph: continue
            dn, dv = g[key] - val, g['vd'] - vd
            sc = (dv / OFS_SIG_VD) ** 2 + (dn / OFS_SIG_ND) ** 2
            p = dpgf(g) if tpg is not None else None
            if tpg is not None:
                sc += ((p - tpg) / OFS_SIG_PGF) ** 2 if p is not None else 25.0
            if g.get('status') == 2: sc += OFS_OBSOLETE
            cand.append({'g': g, 'v': v, 'dn': abs(dn), 'dv': abs(dv), 's': sc, 'mold': mold,
                         'dp': None if p is None else abs(p - tpg), 'old': g.get('status') == 2})
    if not cand: return None, ''
    eps = 1 + 1e-9
    # ① 窗口内（|Δnd|≤0.002 且 |Δvd|≤1）：按 --vendors 顺序取第一家；
    #    非球面元件优先模压料，其次现行牌号，再按 nd 最近（与旧逻辑一致，免得无谓改动）
    win = [x for x in cand if x['dn'] <= OFS_ND_WIN * eps and x['dv'] <= OFS_VD_WIN * eps]
    for v in vendors:
        w = [x for x in win if x['v'] == v]
        if w:
            w.sort(key=lambda x: (not (asph and x['mold']), x['old'], x['dn']))
            return w[0]['g'], '窗口内'
    # ② 窗口外：色散优先（νd、有 θgF 时加 dPgF），nd 只做次要项
    pool = [x for x in cand if x['dn'] <= OFS_ND_CAP * eps] or cand
    top = min(pool, key=lambda x: x['s'])
    for v in vendors:
        pv = [x for x in pool if x['v'] == v]
        if not pv: continue
        b = min(pv, key=lambda x: x['s'])
        # 「合理」按 νd 的绝对差距判，不按分数：各厂 νd 的定义/舍入本来就差 0.1 量级，
        # 分数一平方就把 2.16 vs 2.10 这种等价候选放大成「不合理」（RF28-70 面29 S-BSL7 vs BSC7）
        ok = b['dv'] <= max(OFS_VD_WIN, top['dv'] + OFS_VENDOR_VD) * eps and \
            (b['dp'] is None or b['dp'] <= max(0.003, (top['dp'] or 0) + OFS_VENDOR_PGF))
        if ok:
            return b['g'], '色散优先'
    return top['g'], '色散优先（首选厂家无合理候选）'

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
        # 衍射面（DOE）：φ_DOE = −2·C2（基准波长、m=1），见 vignet.doe_dw
        c2 = ((s.get('doe') or {}).get('C') or [0.0])[0]
        u = (n * u - y * c * (n2 - n) + 2.0 * c2 * y) / n2
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


def solve_obj(g, f, omax=1e9):
    """解 g(o)=0 求物距 o（g = 该物距下的近轴像距 − 目标像距）。

    ★ 不能用「lo = 1.02|f| 起、hi 逐级 ×10」那套固定括号 ——
    物在**前焦点**处像距发散，这个极点落在 lo 与真根之间时，g(lo) 与 g(hi) 同号，
    真根被整段跳过（实测 JP2022-61515A Ex2：真根 o≈6000，lo=89 恰在极点内侧，
    整个对焦解全部返回 None）。
    做法改成：在 (1.001|f|, omax] 上做对数扫描，物理支是极点右侧那一段 ——
    g 在该段从 +∞ 单调降到 bf_inf − BF0(<0)，所以**最右侧的「+ → −」变号**就是真根。
    极点本身是「− → +」的巨幅跳变，天然被排除。
    """
    lo0 = 1.001 * abs(f)
    N = 400
    xs = [lo0 * (omax / lo0) ** (i / N) for i in range(N + 1)]
    vs = [g(o) for o in xs]
    br = None
    for i in range(N):
        a, b = vs[i], vs[i + 1]
        if a == a and b == b and a > 0 >= b:
            br = (xs[i], xs[i + 1])
    if br is None: return None
    return brentq(g, br[0], br[1])


def zoom_configs(emb, zx, betas):
    """变焦镜头的结构表：每个变焦位置一个 ∞ 结构 + 每个 |β| 一个对焦结构。

    spec 写法（见 references/spec-schema.md「变焦」）：
        zmx.zoom = {"positions": [{"name":"W","label":"W","inf":"W-INF","near":["W-MFD"],"fno":2.91}, ...],
                    "betas": [0.06], "include_near": false}
    emb.states / emb.variable 里每个变焦位置的 ∞ 态与专利近距态各一列（所有可变间隔都要写全）。

    做法：
    - **结构里写全所有可变间隔**（变焦间隔 + 对焦间隔 + BF），不再只写 key_before ——
      变焦镜头每个位置的守恒和各不相同（本篇 D18+D21+D23 = 30.483 / 34.060 / 32.862），
      zmx 的位置解 TOLE 长度又进不了 MCE，所以两条出口都按「每个可变间隔逐结构定值」写。
    - 对焦路径：专利给了该位置的近距态 → 在「∞ 态 → 近距态₁ → 近距态₂ …」的间隔空间里**分段线性**走 t。
      两态时这正是 lensmath 定焦分支的凸轮直线（单组 / 双浮动 / 链式三段通吃），中间 0.06x 的分工是插值的。
      没给近距态、但 zmx.focus 是单组（key_before + key_after）→ 在该位置按 ±x 平移那一组。
    - 像面条件对任何对焦方式都成立：**近轴像距 = 该结构末面到像面的间隔**（末面间隔是变量时跟着走），
      再扣掉 ∞ 态「近轴后焦 − 印刷 BF」的舍入差。
    """
    var = emb['variable']
    keys = list(var)
    zz = zx['zoom']
    rows = [q for q in emb['surfaces'] if q['i'] != 'IMG']
    klast = rows[-1]['D'] if isinstance(rows[-1]['D'], str) else None

    def expand(dm):
        out = []
        for s in rows:
            D = s['D']
            if isinstance(D, str): D = dm[D]
            out.append({'R': s['R'], 'D': float(D), 'nd': s.get('nd'), 'doe': s.get('doe')})
        return out

    def dlast(dm):
        return float(dm[klast]) if klast else float(rows[-1]['D'])

    fcs = zx.get('focus') or {}
    inf_cfgs, near_cfgs, pat_cfgs = [], {m: [] for m in betas}, []
    print('\n== 变焦：逐变焦位置 ∞ + 对焦解（%d 个变焦位置 × |β| %s）=='
          % (len(zz['positions']), ', '.join('%.2fx' % m for m in betas)))
    for p in zz['positions']:
        st = p['inf']
        dm0 = {k: float(var[k][st]) for k in keys}
        S0 = expand(dm0)
        f0, bf0, _ = paraxial(S0)
        off = bf0 - dlast(dm0)
        label = p.get('label', p['name'])
        nm0 = '%s %.1fmm' % (label, f0)
        pf = p.get('f')
        print('\n  [%s] %s  EFL=%.4f%s  近轴后焦=%.4f（印刷末面间隔 %.4f，舍入差 %+.4f）  ΣD=%.3f'
              % (p['name'], st, f0, ('（专利 %.2f）' % pf) if pf else '（专利未印）', bf0, dlast(dm0), off,
                 sum(q['D'] for q in S0)))
        inf_cfgs.append(dict([('name', nm0 + ' INF'), ('zoom', p['name']), ('d0', 'INFINITY')]
                             + [(k, dm0[k]) for k in keys] + [('efl', round(f0, 4))]))
        near = list(p.get('near') or [])
        path = [dm0] + [{k: float(var[k][q]) for k in keys} for q in near]
        if len(path) > 1:
            fk = [k for k in keys if any(abs(d[k] - dm0[k]) > 1e-9 for d in path[1:])]
            nseg = len(path) - 1
            def dm_at(t, _path=path, _nseg=nseg):
                j = min(max(int(math.floor(t)), 0), _nseg - 1)
                u = t - j
                return {k: _path[j][k] + (_path[j + 1][k] - _path[j][k]) * u for k in keys}
            t_doc = float(nseg)
            print('       对焦间隔 %s：沿专利 %s 分段线性（凸轮是插值的）' % (fk, ' → '.join([st] + near)))
        elif fcs.get('key_before') and fcs.get('key_after'):
            if zx.get('focus2'):
                print('       ★ 变焦分支没有近距态时只支持单组平移兜底，zmx.focus2 已忽略（双浮动没有近距态定不了凸轮）')
            kb, ka = fcs['key_before'], fcs['key_after']
            fk = [kb, ka]
            def _mk(sg, _kb=kb, _ka=ka, _dm0=dm0):
                def f_(t):
                    d = dict(_dm0); d[_kb] = _dm0[_kb] + sg * t; d[_ka] = _dm0[_ka] - sg * t
                    return d
                return f_
            dm_at = None; t_doc = 0.0
            for sg in (+1.0, -1.0):
                cand = _mk(sg)
                try:
                    Sx = expand(cand(0.05)); fx = paraxial(Sx)[0]
                    o = solve_obj(lambda o_: paraxial(Sx, obj=o_)[1] - (dlast(cand(0.05)) + off),
                                  min(abs(fx), abs(f0)))
                except Exception:
                    o = None
                if o and o > 0:
                    dm_at = cand; break
            if dm_at is None:
                print('       ★ 单组 %s/%s 两个方向都解不出实物距，跳过该位置的对焦结构' % (kb, ka)); continue
            print('       专利没印该位置近距态 → 单组 %s+%s 平移（全部标 ★外推）' % (kb, ka))
        else:
            print('       ★ 专利没印该位置近距态，zmx.focus 也不是单组 → 只出 ∞ 结构'); continue

        def ok(dm, _fk=fk):
            return all(dm[k] >= 0.10 - 1e-9 for k in _fk)

        def solve(t, _dm_at=dm_at, _f0=f0, _off=off):
            dm = _dm_at(t); S = expand(dm)
            fx = paraxial(S)[0] or _f0
            try:
                o = solve_obj(lambda o_: paraxial(S, obj=o_)[1] - (dlast(dm) + _off),
                              min(abs(fx), abs(_f0)))
            except Exception:
                o = None
            if o is None: return None, None, dm
            return paraxial(S, obj=o)[2], o, dm

        # 可走的最远 t：对焦间隔全部 ≥ 0.10mm（专利记载段之外按最后一段外推，标 ★外推）
        # 整数步进（t = n·step，不累加浮点：50 次 +0.02 = 1.0000000000000004 会把正好 0.10 的间隔判出界），
        # 再在最后一步里二分到「最小对焦间隔 = 0.10」的精确位置，与定焦分支一致。
        step = 0.02 if t_doc else 0.05
        tmax = max(t_doc, 1.0) * 4 if t_doc else 200.0
        n = 0
        while (n + 1) * step <= tmax + 1e-12 and ok(dm_at((n + 1) * step)):
            n += 1
        t_hi = n * step
        if (n + 1) * step <= tmax + 1e-12:
            lo_, hi_ = t_hi, (n + 1) * step
            for _ in range(50):
                mid_ = 0.5 * (lo_ + hi_)
                if ok(dm_at(mid_)): lo_ = mid_
                else: hi_ = mid_
            t_hi = lo_
        N = 600
        tab = []
        for i in range(1, N + 1):
            t = t_hi * i / N
            b, o, _dm = solve(t)
            if o: tab.append((t, abs(b)))
        for j in range(1, int(t_doc) + 1):
            b, o, dm = solve(float(j))
            q = near[j - 1]
            # near_d0_printed 是专利表头那个「近距離(165mm)」—— 只属于最后一个近距态（最短撮影距離）
            given = (emb.get('d0') or {}).get(q) or (zz.get('near_d0_printed') if j == int(t_doc) else None)
            print('       %-8s 专利近距态  物距(面1起) %s  β=%+.5f  撮影距離(像面起) %.1fmm%s'
                  % (q, ('%.2f' % o) if o else '解不出', b or 0.0,
                     (o or 0) + sum(x['D'] for x in expand(dm)),
                     ('   ← 专利印 %s' % given) if given else ''))
            if zz.get('include_near') and o:
                # 最后一个近距态 = MFD；中间态用专利的状态名，避免结构重名（fno_patent / wfno_override 按名字查）
                nmq = ('%s MFD' % nm0) if j == int(t_doc) else ('%s %s' % (nm0, q))
                pat_cfgs.append(dict([('name', nmq), ('zoom', p['name']), ('d0', round(o, 4))]
                                     + [(k, round(dm[k], 5)) for k in keys]
                                     + [('beta', round(b, 5)), ('efl', round(f0, 4))]))
        for m in betas:
            pr = [(u, v) for u, v in zip(tab, tab[1:]) if (u[1] - m) * (v[1] - m) <= 0]
            if not pr:
                print('       %.2fx: 超出对焦行程（|β| 最大 %.4f）' % (m, max([x[1] for x in tab] or [0]))); continue
            (t1, b1), (t2, b2) = pr[0]
            for _ in range(60):
                tm = 0.5 * (t1 + t2); bm = abs(solve(tm)[0] or 0.0)
                if (b1 - m) * (bm - m) <= 0: t2, b2 = tm, bm
                else: t1, b1 = tm, bm
            t = 0.5 * (t1 + t2)
            b, o, dm = solve(t)
            ext = t > t_doc + 1e-9
            print('       %-8s t=%.4f  物距 %.2f  β=%+.5f  %s%s'
                  % ('%.2fx' % m, t, o, b, '  '.join('%s=%.4f' % (k, dm[k]) for k in fk),
                     '   ← ★外推：超出专利记载的对焦范围' if ext else ''))
            near_cfgs[m].append(dict([('name', '%s %.2fx' % (nm0, m)), ('zoom', p['name']), ('d0', round(o, 4))]
                                     + [(k, round(dm[k], 5)) for k in keys]
                                     + [('beta', round(b, 5)), ('efl', round(f0, 4)), ('extrapolated', ext)]))
    cfgs = inf_cfgs + [c for m in betas for c in near_cfgs[m]] + pat_cfgs
    print('\n  变焦交付结构（%d 个）：%s' % (len(cfgs), ' / '.join(c['name'] for c in cfgs)))
    return cfgs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('spec'); ap.add_argument('--vendors', nargs='+', default=['HOYA', 'OHARA', 'CDGM'])
    ap.add_argument('--brand-slack', type=float, default=BRAND_SLACK,
                    help='Δnd 落在最优值+该带宽内即按厂家优先级选（默认 0.0020，设 0 关闭）')
    ap.add_argument('--betas', nargs='*', type=float, default=None,
                    help='对焦结构的 |β|。默认：定焦 0.02 0.06（再加 MFD 态 = INF/0.02x/0.06x/MFD 四结构）；'
                         '变焦 zmx.zoom.betas，没写就 0.06（W/M/T 各 ∞ + 0.06x = 六结构）')
    ap.add_argument('--match-tol', type=float, default=2e-5,
                    help='d 线时认作「等效牌号」的 |Δnd| 上限（默认 2e-5）；超出则保留印刷 nd/vd 且不写 glass')
    ap.add_argument('--alt-only', nargs='*', default=[],
                    help='这些厂家只进备选列、永不做首选。日系申请人固定传 --alt-only CDGM')
    ap.add_argument('--mfd', type=float, default=None,
                    help='产品标称最短撮影距離(mm，物体→像面)。给了就额外解一个 MFD 状态：'
                         '物距 = MFD − ΣD。专利近距态离产品 MFD 很远时用它。')
    ap.add_argument('--allow-legacy', action='store_true',
                    help='允许含铅/停产牌号做候选。**只在 2000 年前的老专利上用** —— '
                         '默认排除 OHARA 的 PBM/PBH/BPH/BAL/BSM… 这些含铅老系列与 AGF 里 status=Obsolete 的牌号。')
    ap.add_argument('--emb', type=int, default=0)
    ap.add_argument('--write', nargs='?', const='AUTO', default=None,
                    help='写出匹配后的 spec（默认 <spec>.matched.json，不覆盖原文件）')
    a = ap.parse_args()
    spec = json.load(open(a.spec, encoding='utf-8'))
    emb = spec['embodiments'][a.emb]
    _zx0 = spec.get('zmx') or {}
    if a.betas is None:
        a.betas = list((_zx0.get('zoom') or {}).get('betas') or [0.06]) if _zx0.get('zoom') else [0.02, 0.06]
    if a.mfd is not None and _zx0.get('zoom'):
        print('★ --mfd 只对定焦生效；变焦的近距结构用 zmx.zoom.betas / zmx.zoom.include_near')
    if a.mfd is None and _zx0.get('mfd') and not _zx0.get('zoom'):
        # 产品标称 MFD 写在 spec 里（zmx.mfd，自像面起算 mm）就不必每次敲 --mfd
        a.mfd = float(_zx0['mfd'])
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
        c = best(s['nd'], s['vd'], libs, line, asph, a.brand_slack, a.match_tol,
                 tuple(a.alt_only), eco=not a.allow_legacy, mold_ok=asph_eff)
        if not c:        # 无铅池里一个都没有 —— 退回全目录并告警，绝不静默
            c = best(s['nd'], s['vd'], libs, line, asph, a.brand_slack, a.match_tol,
                     tuple(a.alt_only), eco=False, mold_ok=asph_eff)
            if c: print('  ★ 面%s 无铅目录里没有候选，回退到含铅/停产牌号 %s'
                        % (s['i'], c[0][2]['name']))
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
        # 无等效牌号的面要另挑一颗「贴合实物」的基准玻璃做 Offset 解（规则见 offset_base）：
        #   按 --vendors 顺序、球面元件不用模压料、窗口外色散优先。
        base, base_how = g, ''
        if abs(g[line] - s['nd']) > a.match_tol:
            bb, base_how = offset_base(s['nd'], s['vd'], libs, line, asph_eff, tuple(a.alt_only),
                                       eco=not a.allow_legacy, pgf=s.get('pgf'))
            if bb is not None: base = bb
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
                                s['glass_offset'], base_how))
    if NOMATCH:
        print('\n  ★ 无等效牌号（保留专利印刷值，zmx 该面走模型玻璃）:')
        for i, nd, vd, v, n, dn, go, how in NOMATCH:
            print('     面%-4s %.5f/%.1f  最近 %s %s Δnd%+.5f' % (i, nd, vd, v, n, dn))
            print('            → Offset 解：基准 %-18s nd=%.5f νd=%.4f  Nd offset %+.5f  Vd offset %+.4f  [%s]'
                  % (go['base'], go['base_nd'], go['base_vd'], go['d_nd'], go['d_vd'], how))
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
            out.append({'R': s['R'], 'D': float(D), 'nd': s.get('nd'), 'doe': s.get('doe')})
        return out
    base = expand(st)
    SD = sum(s['D'] for s in base)
    f, bf, _ = paraxial(base)
    gen = dict((k, v) for k, v in emb.get('general', []))
    print('\n== 近轴校验 ==')
    print('  ΣD = %.2f   (专利 L = %s)' % (SD, gen.get('L 光学全长 (mm)', '?')))
    print('  EFL = %.4f  BF = %.4f   (专利 f = %s)' % (f, bf, gen.get('f (mm)', '?')))
    def _var_on_last_surface(key):
        """可変間隔是否就是「最終レンズ面 → 像面」那一段。"""
        rows = [q for q in emb['surfaces'] if q['i'] != 'IMG']
        return bool(rows) and rows[-1]['D'] == key

    fc2 = zx.get('focus2')
    if zx.get('zoom'):
        # ===== 变焦镜头：W/M/T 各 ∞ + 各 |β| 对焦 =====
        cfgs = zoom_configs(emb, zx, a.betas)
        if a.write: spec.setdefault('zmx', {})['configs'] = cfgs
    elif not fc: print('\n(spec 无 zmx.focus，跳过对焦解)');
    elif not fc2 and fc.get('key_after') is None and _var_on_last_surface(fc['key_before']):
        # ===== 整組繰り出し：可変間隔 = 最終レンズ面〜像面 =====
        # 紧凑型定焦常见（本例 JP2023-140823A）：对焦时整个镜筒相对像面前伸，
        # 唯一在变的就是最后那一段空气。
        # ⚠ 它**不是内部间隔** —— 改它不会改变光学系统本身，所以通用分支那套
        # 「把像面钉在 ∞ 态后焦上、解物距」的方程恒等于 0=0，会把所有倍率
        # 都误报成「超出对焦行程」。这里反过来做：给物距 → 近轴像距**就是**该状态的间隔。
        kb = fc['key_before']
        d0v = float(emb['variable'][kb][st])
        base_s = expand(st, {kb: d0v})
        f0, bf0, _ = paraxial(base_s)
        SD_fix = sum(q['D'] for q in base_s) - d0v       # 面1 → 最終面 的轴上长度
        cap = float(fc.get('sum') or 0.0)                # 可変間隔的机械上限（只用来告警）
        bf_of = lambda o: paraxial(base_s, obj=o)[1]
        bet_of = lambda o: paraxial(base_s, obj=o)[2]

        def _root(g, N=2000):
            """在 (1.001|f|, 1e9] 上对数扫描，取最右侧的「+ -> -」变号。
            物在前焦点处像距发散，极点是「- -> +」的巨跳，天然排除。"""
            lo0 = 1.001 * abs(f0)
            xs = [lo0 * (1e9 / lo0) ** (i / float(N)) for i in range(N + 1)]
            gs = [g(x) for x in xs]
            br = None
            for i in range(N):
                if gs[i] is None or gs[i + 1] is None: continue
                if gs[i] > 0 >= gs[i + 1]: br = (xs[i], xs[i + 1])
            if not br: return None
            lo, hi = br
            for _ in range(200):
                mid = 0.5 * (lo + hi)
                if g(lo) * g(mid) <= 0: hi = mid
                else: lo = mid
            return 0.5 * (lo + hi)

        print('\n== 对焦位置解（整組繰り出し：面1〜面%d 全体前伸、%s = 最終面〜像面）=='
              % (len(base_s), kb))
        print('  ★ 专利只印 ∞ 一态，以下近距状态全部是近轴反解，一律标 ★外推 —— '
              '专利没有背书那里的像差校正。')
        print('  %-14s d0=∞          %s=%.4f  β=0' % (st, kb, d0v))
        cfgs = [{'name': st, 'd0': 'INFINITY', kb: round(d0v, 4)}]

        def _emit(nm, o):
            if o is None:
                print('  %-14s 超出对焦行程' % nm); return
            x = bf_of(o); b = bet_of(o)
            warn = '   ← ★ 超出设定的行程上限 %.2f' % cap if cap and x > cap else ''
            print('  %-14s d0=%-10.2f %s=%.4f  β=%+.5f  (1:%.1f)  撮影距離=%.1fmm  ★外推%s'
                  % (nm, o, kb, x, b, 1 / abs(b) if b else 0, o + SD_fix + x, warn))
            cfgs.append({'name': nm, 'd0': round(o, 2), kb: round(x, 4),
                         'extrapolated': True})

        for m in a.betas:
            _emit('%.2fx' % m, _root(lambda o, _m=m: abs(bet_of(o)) - _m))
        if a.mfd:
            # 撮影距離（像面起算）= 物距 + 面1→最終面 + 最終面→像面
            _emit('MFD(%gmm)' % a.mfd,
                  _root(lambda o: (a.mfd - (o + SD_fix + bf_of(o)))))
        if a.write: spec.setdefault('zmx', {})['configs'] = cfgs
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
            # ★ 括号起点必须用**该结构自己的 EFL**，不是 ∞ 态的 f：
            # 高倍率微距对焦后 EFL 会大幅下降（本篇 100.81 → 36.01），
            # 物距可以小于 ∞ 态的 f，用 1.001|f_inf| 起扫会把真根整段跳过（d0 变 nan）。
            fx = bfd2(x)[0] or f
            return solve_obj(lambda o: bfd2(x, o)[1] - target(x), min(abs(fx), abs(f)))
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
            sd_of = lambda x: sum(q['D'] for q in expand(st, dmap_of(x)))
            xs_scan = [t[0] for t in tab]
            o_of = lambda x: (beta2(x)[1] or float('inf'))
            gg = lambda x: o_of(x) - (a.mfd - sd_of(x))
            pr = [p for p in zip(xs_scan, xs_scan[1:]) if gg(p[0]) * gg(p[1]) <= 0]
            if not pr:
                print('  MFD(%gmm): 超出对焦行程' % a.mfd)
            else:
                x1, x2 = pr[0]
                for _ in range(60):
                    xm = 0.5 * (x1 + x2)
                    if gg(x1) * gg(xm) <= 0: x2 = xm
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
            # 专利记载的每个近距态都要进 configs（原来只收 sts[-1]，
            # 中间态如 0.5倍 会被丢掉，而那是专利原值、比插值更可信）
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
        # key_after 就是**最后一面到像面**的间隔（BF）时，对焦组 = 光阑以后整个后组前移、像面固定
        # （US20100208366A1 TS-E24 II：D12 减、BF 增，D12+BF 恒定）。这时像距目标不是 ∞ 的 BF0，
        # 而是跟着 BF 一起变：BF0 + (ka − ka_∞)。旧写法把像距钉在 BF0，任何倍率都「超出对焦行程」。
        _real = [q for q in emb['surfaces'] if q['i'] != 'IMG']
        ka_is_bf = bool(ka) and _real and _real[-1].get('D') == ka
        d0_inf = emb['variable'][kb][st]
        tgt = (lambda d13: BF0 + (d0_inf - d13)) if ka_is_bf else (lambda d13: BF0)

        def solve_d0(d13):
            """像面固定在 ∞ 近轴焦点，解物距；解不出返回 None（等于物在无穷远之外）。"""
            fx = bfd(d13)[0] or f      # 同上：用该结构自己的 EFL 定括号起点
            return solve_obj(lambda o: bfd(d13, o)[1] - tgt(d13), min(abs(fx), abs(f)))

        def beta(d13):
            o = solve_d0(d13)
            if o is None: return 0.0, None
            return bfd(d13, o)[2], o

        d0v = emb['variable'][kb][st]
        multi = len(emb.get('states', [])) > 1
        if multi:
            near = emb['variable'][kb][emb['states'][-1]]
        elif ka is None:
            # 整組繰り出し式（key_after 为 null）：对焦组只朝物侧走，可变间隔单调增大。
            # 专利只印 ∞ 一态时不能用 0.2*d0v 当近距端，那会把扫描方向整个弄反。
            near = tot - 0.10
        else:
            # 专利只印 ∞ 一态时方向不能预设：内对焦组可能朝物侧也可能朝像侧走
            # （EP4215968A1 的负 L2 朝像侧走，预设 0.2*d0v 会整段扫反）。
            # 两侧各探一步，取解出**正的实物距**的那一侧。
            near = 0.2 * d0v
            if not multi:
                def _ok(x):
                    try:
                        _b, _o = beta(x)
                        return bool(_o) and _o > 0
                    except Exception:
                        return False
                up = min(d0v + 0.05 * (tot - d0v), tot - 0.10)
                if _ok(up) and not _ok(0.5 * d0v):
                    near = tot - 0.10
        # 扫描范围要按**机械行程**给，不能只给「专利记载的 ∞↔近距」那一段 ——
        # 有的专利（如 US20150092100 的近距态只到 0.033x）记载范围比镜头实际行程短得多，
        # 只扫记载范围会把 0.06x 误判成「超出对焦行程」。
        # 方向仍锁在对焦组实际移动的那一侧，避免扫到另一支非物理的解。
        p_lo, p_hi = (near, d0v) if near < d0v else (d0v, near)
        lo_d, hi_d = (d0v, tot - 0.10) if near > d0v else (0.10, d0v)
        if not multi:
            p_lo = p_hi = d0v      # 专利只记载 ∞ 一态 → 其余解一律标 ★外推
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
        # --mfd: 按产品标称最短撮影距離反解一态（物距 = MFD - ΣD，均自像面起算）
        if a.mfd:
            # ΣD 不一定是常数：整組繰り出し式对焦（key_after 为 null）时两个可变间隔没有
            # 配对守恒，全长随对焦状态变化，拿 ∞ 态的 ΣD 当常数会把 MFD 解偏一整个行程量。
            sd_of = lambda x: sum(q['D'] for q in expand(st, {kb: x, ka: tot - x}))
            gg = lambda x: o_of_(x) - (a.mfd - sd_of(x))
            xs_scan = [x for x, _ in sorted(tab, key=lambda t: t[0])]
            o_of_ = lambda x: (beta(x)[1] or float('inf'))
            o_of = o_of_
            pr = [q for q in zip(xs_scan, xs_scan[1:]) if gg(q[0]) * gg(q[1]) <= 0]
            if not pr:
                print('  MFD(%gmm): 超出对焦行程' % a.mfd)
            else:
                x1, x2 = pr[0]
                for _ in range(60):
                    xm = 0.5 * (x1 + x2)
                    if gg(x1) * gg(xm) <= 0: x2 = xm
                    else: x1 = xm
                x = 0.5 * (x1 + x2); b, o = beta(x)
                ext = '' if p_lo - 1e-9 <= x <= p_hi + 1e-9 else '   ← ★外推：超出专利记载的对焦范围'
                nm = 'MFD(%gmm)' % a.mfd
                print('  %-12s d0=%-10.2f %s=%.4f %s=%.4f  β=%+.5f  (1:%.1f)%s'
                      % (nm, o, kb, x, ka, tot - x, b, 1 / abs(b) if b else 0, ext))
                cfgs.append({'name': nm, 'd0': round(o, 2), kb: round(x, 4),
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
    cf = (spec.get('zmx') or {}).get('configs') if a.write else None
    if cf and not zx.get('zoom'):
        # 定焦交付约定：INF / 0.02x / 0.06x / MFD 四个结构（用户 2026-09 定的）
        nms = [c.get('name') for c in cf]
        has_mfd = len(nms) >= 2 and not re.match(r'^\d+\.\d+x$', str(nms[-1]))
        print('\n  定焦结构（约定 INF / 0.02x / 0.06x / MFD）：%s%s'
              % (' / '.join(map(str, nms)),
                 '' if has_mfd else '   ★ 缺 MFD 态：专利没印近距态时在 spec 写 zmx.mfd（产品标称，像面起算 mm）或传 --mfd'))
    if a.write:
        out = a.spec[:-5] + '.matched.json' if a.write == 'AUTO' else a.write
        json.dump(spec, open(out, 'w', encoding='utf-8'), ensure_ascii=False, indent=1)
        print('\n已写出 %s（原 spec 未改动）' % out)

if __name__ == '__main__':
    main()
