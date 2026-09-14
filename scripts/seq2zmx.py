#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
seq2zmx.py —— CODE V 序列文件 .seq  →  Zemax .zmx（与 make_seq.py 相反的方向）

    python3 seq2zmx.py A2628.seq -o A2628.zmx \
        [--gcat "CDGM=CDGM2025011,HOYA=HOYA20260707"] [--glassdir DIR] \
        [--reverse-fields] [--raim 2] [--no-clap]

对照表见 references/codev-seq.md（那份是 zmx→seq，反过来读即可）。几条要点：

* 渐晕要反算，不是照抄：VCY=(VUY+VLY)/2、**VDY=(VLY−VUY)/2（带负号）**、VCX=VUX(=VLX)。
  CODE V 自带的「导出 Zemax」把 VDY 写成了 +(VUY−VLY)/2，符号是反的。
* `ZOO FNO` 那一串（2.9→3.45）是**有限共轭下的工作 F 数**，入瞳其实恒定。
  zmx 只写一行 `FNUM`，**绝不能**逐结构铺 APER —— 那会把光瞳二次缩小。
* CODE V 的 `CIR` 是**挡光的实口径** → `DIAM <v> 1`（固定）+ `CLAP 0 <v> 0`；没写 CIR 的面留自动。
* 玻璃名反查：`DQK3L_CDGM` → 去掉厂家后缀，按「去掉所有非字母数字」与目录牌号比对 → `D-QK3L`。
  查不到就报错退出，**不要**默默退化成模型玻璃（CODE V 自带导出就是全退成 nd=1.5/vd=0 的废文件）。
* 视场顺序默认**原样保留**（CODE V 轴上在前），这样两边逐项对得上；要 Zemax 习惯用 --reverse-fields。
* CCY/THC 变量标记、PIK 拾取解、CMP/DSX/DSY/BTX/BTY 公差数据**不转**，会在结束时逐条列出来。
"""
import re, os, sys, math, argparse

# ---------------- 玻璃目录 ----------------
def load_cat(path):
    b = open(path, 'rb').read()
    if b[:2] in (b'\xff\xfe', b'\xfe\xff'):
        txt = b.decode('utf-16')
    else:
        try:    txt = b.decode('utf-8')
        except UnicodeDecodeError: txt = b.decode('latin-1')
    out = {}
    if path.lower().endswith('.csv'):
        for ln in txt.splitlines()[1:]:
            f = ln.split(',')
            if len(f) >= 3:
                try: out[f[0].strip()] = (float(f[1]), float(f[2]))
                except ValueError: pass
        return out
    for ln in txt.splitlines():
        if ln.startswith('NM '):
            f = ln.split()
            try: out[f[1]] = (float(f[4]), float(f[5]))
            except (IndexError, ValueError): pass
    return out

def cat_dirs(extra):
    d = []
    if extra: d.append(extra)
    if os.environ.get('PATENT_GLASS_DIR'): d.append(os.environ['PATENT_GLASS_DIR'])
    d.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'assets', 'glass'))
    return d

CATS = {}
def vendor_cat(vend, extra):
    if vend in CATS: return CATS[vend]
    for d in cat_dirs(extra):
        if not os.path.isdir(d): continue
        for fn in sorted(os.listdir(d)):
            stem, ext = os.path.splitext(fn)
            if ext.lower() not in ('.agf', '.csv'): continue
            if stem.upper() == vend or stem.upper().startswith(vend):
                CATS[vend] = (load_cat(os.path.join(d, fn)), stem)
                return CATS[vend]
    CATS[vend] = ({}, None)
    return CATS[vend]

def key(s): return re.sub(r'[^A-Z0-9]', '', s.upper())

# ---------------- 读 .seq ----------------
def read_seq(path):
    raw = open(path, encoding='latin-1').read().replace('\r\n', '\n')
    lines = []
    for ln in raw.split('\n'):
        if lines and lines[-1].endswith('&'):
            lines[-1] = lines[-1][:-1] + ' ' + ln.strip()
        else:
            lines.append(ln.rstrip())
    return [l for l in lines if l.strip()]

def parse(lines):
    hdr, surfs, zoo, extras = {}, [], {}, []
    nz, cur = 1, None
    for ln in lines:
        t = ln.strip()
        if not ln.startswith(' '):                       # 顶格命令
            m = re.match(r'^(S|SO|SI)\s+(-?[\d.eE+-]+)\s+(-?[\d.eE+-]+)(.*)$', ln)
            if m:
                rest = m.group(4).split()
                s = {'kind': m.group(1), 'rdy': float(m.group(2)), 'thi': float(m.group(3)),
                     'glass': None, 'model': None, 'asp': None, 'cir': None,
                     'sto': False, 'oal': None}
                if len(rest) == 1:
                    s['glass'] = rest[0]
                elif len(rest) >= 2:
                    try:    s['model'] = (float(rest[0]), float(rest[1]))   # CODE V 模型玻璃
                    except ValueError: s['glass'] = rest[0]
                surfs.append(s); cur = s; continue
            m = re.match(r'^TIT\s+Z(\d+)\s+(.*)$', t)
            if m:
                zoo.setdefault('TITLES', {})[int(m.group(1))] = m.group(2).strip().strip('"').strip("'").strip('"')
                continue
            if t.startswith('ZOO '):
                p = t.split()
                if len(p) == 2 and p[1].isdigit(): nz = int(p[1]); continue
                if p[1] == 'TIT': continue
                if p[1] in ('THI', 'THC', 'VUY', 'VLY', 'VUX', 'VLX', 'WTF'):
                    zoo[' '.join(p[1:3])] = p[3:]
                else:
                    zoo[p[1]] = p[2:]
                continue
            m = re.match(r'^THI\s+S(\d+)\s+OAL\s+S(\d+)\.\.(\d+)\s+(-?[\d.eE+-]+)', t)
            if m:
                surfs[int(m.group(1))]['oal'] = (int(m.group(2)), float(m.group(4))); continue
            if t.split()[0] in ('PIK', 'CMP', 'DSX', 'DSY', 'BTX', 'BTY', 'DER', 'CLS', 'UID', 'RDM;LEN', 'GO', 'INI'):
                if t.split()[0] not in ('UID', 'DER', 'RDM;LEN', 'GO', 'INI'): extras.append(t)
                continue
            k = t.split()[0]
            if k in ('FNO','DIM','WL','REF','WTW','XRI','YRI','WTF','VUX','VLX','VUY','VLY','TITLE','DOR'):
                hdr[k] = t[len(k):].strip()
            continue
        if cur is None: continue                          # 面属性
        if   t.startswith('CIR'): cur['cir'] = float(t.split()[1])
        elif t == 'STO':          cur['sto'] = True
        elif t == 'ASP':          cur['asp'] = {'K': 0.0, 'c': [0.0]*9}
        elif t.startswith('K ') and cur['asp'] is not None:
            cur['asp']['K'] = float(t.split(';')[0].split()[1])
        elif cur['asp'] is not None:
            for part in t.split(';'):
                m = re.match(r'^\s*([A-HJ])\s+(-?[\d.]+(?:[eE]-?\d+)?)\s*$', part)
                if m:
                    k_ = m.group(1); i = 'ABCDEFGH'.find(k_)
                    cur['asp']['c'][i if i >= 0 else 8] = float(m.group(2))
    return hdr, surfs, zoo, nz, extras

# ---------------- 写 .zmx ----------------
def num(x, f='%.15G'): return f % x

def build(seq, args):
    hdr, S, zoo, nz, extras = parse(read_seq(seq))
    yfl = [float(x) for x in hdr['YRI'].split()]
    wtf = [float(x) for x in hdr.get('WTF', ' '.join(['1']*len(yfl))).split()]
    wl  = [float(x) for x in hdr['WL'].split()]
    wtw = [float(x) for x in hdr.get('WTW', ' '.join(['1']*len(wl))).split()]
    ref = int(float(hdr.get('REF', 1)))
    nf  = len(yfl)
    order = list(range(nf))[::-1] if args.reverse_fields else list(range(nf))

    def vg(cfg, j):   # j = CODE V 视场下标(0起)
        def g(tag):
            k = '%s F%d' % (tag, j+1)
            if k in zoo: return float(zoo[k][cfg-1])
            v = hdr.get(tag)
            return float(v.split()[j]) if v else 0.0
        z = lambda x: 0.0 if abs(x) < 1e-7 else x
        vuy, vly, vux = z(g('VUY')), z(g('VLY')), z(g('VUX'))
        return ((vly-vuy)/2.0, (vuy+vly)/2.0, vux)        # VDY, VCY, VCX

    gmap = dict(p.split('=') for p in args.gcat.split(',')) if args.gcat else {}
    gcat, notes, used = [], [], {}
    def zglass(s):
        cv = s['glass']; nm, vend = cv.rsplit('_', 1)
        cat, stem = vendor_cat(vend.upper(), args.glassdir)
        for k, v in cat.items():
            if key(k) == key(nm):
                name = gmap.get(vend.upper(), stem)
                if name and name not in gcat: gcat.append(name)
                used[cv] = (k, v[0], v[1])
                return k, v[0], v[1]
        sys.exit('!! %s 目录里找不到牌号 %s（先把该厂家的 .AGF 放进 assets/glass 或 PATENT_GLASS_DIR）' % (vend, nm))

    L = []; a = L.append
    a('VERS 240529 1 1 L000001')
    a('MODE SEQ')
    a('NAME %s' % (args.name or os.path.splitext(os.path.basename(seq))[0]))
    a('NOTE 0 Converted from CODE V sequence file %s by seq2zmx.py' % os.path.basename(seq))
    a('NOTE 1 ""')
    a('PFIL 0 0 0'); a('LANG 0'); a('UNIT MM X W X CM MR CPMM')
    a('FNUM %s 0' % num(float(hdr['FNO'])))
    a('ENVD 20 1 0'); a('GFAC 0 0')
    for s in S:
        if s['glass']: zglass(s)
    a('GCAT ' + ' '.join(gcat))
    a('RAIM 0 %d 1 1 0 0 0 0 0 1' % args.raim)
    a('PUSH 0 0 0 0 0 0'); a('SDMA 0 1 0'); a('OMMA 1 1')
    a('FTYP 3 0 %d %d 0 0 0 %d' % (nf, len(wl), nf))
    pad = lambda v, m=12: list(v) + [0.0]*(m-len(v))
    a('XFLN ' + ' '.join(num(x) for x in pad([0.0]*nf)))
    a('YFLN ' + ' '.join(num(yfl[i]) for i in order) + ' 0 0')
    a('FWGN ' + ' '.join(num(wtf[i]) for i in order) + ' 0 0')
    v1 = [vg(1, i) for i in order]
    a('VDXN ' + ' '.join(num(x) for x in pad([0.0]*nf)))
    a('VDYN ' + ' '.join(num(x[0]) for x in v1) + ' 0 0')
    a('VCXN ' + ' '.join(num(x[2]) for x in v1) + ' 0 0')
    a('VCYN ' + ' '.join(num(x[1]) for x in v1) + ' 0 0')
    a('VANN ' + ' '.join(num(x) for x in pad([0.0]*nf)))
    for i, (lam, w) in enumerate(zip(wl, wtw), 1):
        a('WAVM %d %s %s' % (i, num(lam/1000.0), num(w)))
    for i in range(len(wl)+1, 25): a('WAVM %d 0.55 0' % i)
    a('PWAV %d' % ref)
    a('POLS 1 0 1 0 0 1 0'); a('GLRS 1 0')
    a('GSTD 0 100.000 100.000 100.000 100.000 100.000 100.000 0 1 1 0 0 1 1 1 1 1 1')
    a('NSCD 100 500 0 0.001 10 9.9999999999999995e-07 0 0 0 0 0 0 1000000 0 2')
    a('COFN QF "COATING.DAT" "SCATTER_PROFILE.DAT" "ABG_DATA.DAT" "PROFILE.GRD"')
    a('COFN COATING.DAT SCATTER_PROFILE.DAT ABG_DATA.DAT PROFILE.GRD')

    for idx, s in enumerate(S):
        a('SURF %d' % idx)
        if s['sto']: a('  STOP')
        # CODE V 的 ASP 带到 r^20（A..D / E..H / J）。Zemax 的 Even Asphere 只到 r^16，
        # 所以 H(r^18)/J(r^20) 一旦非零就必须换 Extended Asphere（TYPE XASPHERE），
        # 系数改放 Extra Data：XDAT 1=项数(10)、2=归一化半径(1.0)、3=r^2(留 0)、4..12=r^4..r^20。
        # 格式实证自用户机器上 5 个真文件与 cv2zmx 宏，见 references/zmx-format.md。
        hi = bool(s['asp']) and (s['asp']['c'][7] or s['asp']['c'][8])
        ztyp = 'STANDARD' if not s['asp'] else ('XASPHERE' if hi else 'EVENASPH')
        a('  TYPE ' + ztyp)
        a('  CURV %s 0 0 0 0 ""' % num(0.0 if s['rdy'] == 0 else 1.0/s['rdy'], '%.16G'))
        a('  HIDE 0 0 0 0 0 0 0 0 0 0 0 0'); a('  MIRR 2 0'); a('  SLAB %d' % (idx+1))
        if ztyp == 'EVENASPH':
            a('  PARM 1 0')
            for p in range(2, 9): a('  PARM %d %s' % (p, num(s['asp']['c'][p-2])))
        elif ztyp == 'XASPHERE':
            xd = '  XDAT %d %.12E 0 0 1.000000000000E+00 0.000000000000E+00 0 ""'
            a(xd % (1, 10.0)); a(xd % (2, 1.0)); a(xd % (3, 0.0))
            for j in range(9): a(xd % (j + 4, float(s['asp']['c'][j])))
            notes.append('面%d 的 H(r^18)/J(r^20) 非零 → 已写成 Extended Asphere（XASPHERE），'
                         '系数在 Extra Data Editor 里' % idx)
        a('  DISZ ' + ('INFINITY' if s['kind'] == 'SO' and s['thi'] > 1e9 else
                       ('0' if s['kind'] == 'SI' else num(s['thi']))))
        if s['oal']:
            a('  TOLE %d %s' % (s['oal'][0], num(s['oal'][1])))
        if s['glass']:
            gn, nd, vd = zglass(s)
            a('  GLAS %s 0 0 %s %s 0 0 0 0 0 0' % (gn, num(nd), num(vd)))
        elif s['model']:
            a('  GLAS ___BLANK 1 0 %s %s 0 0 0 0 0 0' % (num(s['model'][0]), num(s['model'][1])))
            notes.append('面%d 是 CODE V 模型玻璃 → zmx 模型玻璃，dPgF 留 0（二级光谱不准，建议配真实牌号）' % idx)
        if s['asp']: a('  CONI %s' % num(s['asp']['K']))
        if s['cir'] is not None:
            a('  DIAM %s 1 0 0 1 ""' % num(s['cir'])); a('  MEMA %s 0 0 0 1 ""' % num(s['cir']))
        else:
            a('  DIAM 0 0 0 0 1 ""'); a('  MEMA 0 0 0 0 1 ""')
        a('  POPS 0 0 0 0 0 0 0 0 1 1 1 1 0 0 0 0')
        if s['cir'] is not None and s['kind'] == 'S' and not args.no_clap:
            a('  CLAP 0 %s 0' % num(s['cir']))

    if nz > 1:
        ml = lambda op, sn, cf, v: ('%-4s %3d %3d %.12E 0 0 0 1 1 1.000000000000E+00'
                                    ' 0.000000000000E+00 0 "" 0' % (op, sn, cf, v))
        a('TOL TOFF   0   0 0.0000000000000000E+00 0.0000000000000000E+00   0 0 0 0 0')
        a('MNUM %d 1' % nz)
        tit = zoo.get('TITLES', {})
        for c in range(1, nz+1):
            a('LTTL   0 %3d "%s" 0 0 0 1 1 0 0.0 "" 0' % (c, tit.get(c, 'CFG%d' % c)))
        for k in sorted([x for x in zoo if x.startswith('THI S')], key=lambda x: int(x[5:])):
            sn = int(k[5:])
            if S[sn].get('oal'):
                notes.append('面%d 带 OAL/位置解，已跳过它的 ZOO THI（解与逐结构定值不能同时写）' % sn)
                continue
            for c in range(1, nz+1):
                v = float(zoo[k][c-1])
                a(ml('THIC', sn, c, 1e10 if (sn == 0 and v > 1e9) else v))
        for n_, j in enumerate(order, 1):
            for op, ix in (('FVDY', 0), ('FVCY', 1), ('FVCX', 2)):
                vals = [vg(c, j)[ix] for c in range(1, nz+1)]
                if all(abs(x) < 1e-12 for x in vals): continue
                for c in range(1, nz+1): a(ml(op, n_, c, vals[c-1]))
        if 'FNO' in zoo:
            notes.append('ZOO FNO %s 未写成 APER：那是有限共轭的工作 F 数，入瞳恒定，逐结构改 APER 会二次缩小光瞳'
                         % ' '.join(zoo['FNO']))
        a('CONF 1')

    out = args.out or os.path.splitext(seq)[0] + '.zmx'
    open(out, 'w', newline='\r\n', encoding='latin-1').write('\n'.join(L) + '\n')
    print('写出 %s：%d 面, %d 结构, %d 视场, %d 波长, GCAT %s'
          % (out, len(S)-2, nz, nf, len(wl), ' '.join(gcat)))
    for cv, (k, nd, vd) in used.items():
        print('   玻璃 %-16s -> %-12s nd=%.5f vd=%.2f' % (cv, k, nd, vd))
    tol = [e for e in extras if e.split()[0] in ('CMP', 'DSX', 'DSY', 'BTX', 'BTY')]
    if tol:
        notes.append('未转换：%d 条公差/补偿记录（CMP/DSX/DSY/BTX/BTY）—— Zemax 用自己的公差编辑器' % len(tol))
    for e in extras:
        if e.split()[0] == 'PIK':
            notes.append('未转换的拾取解（数值已经是当前值，改上游不会跟随）：%s' % e)
    for n_ in dict.fromkeys(notes): print('   ★ ' + n_)
    return out

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('seq'); ap.add_argument('-o', '--out')
    ap.add_argument('--name'); ap.add_argument('--glassdir')
    ap.add_argument('--gcat', help='厂家=目录文件名，逗号分隔，如 "CDGM=CDGM2025011,HOYA=HOYA20260707"')
    ap.add_argument('--reverse-fields', action='store_true', help='视场倒成 Zemax 习惯（最大视场在前）')
    ap.add_argument('--raim', type=int, default=2, help='Ray Aiming 0=Off 1=Paraxial 2=Real（默认 2）')
    ap.add_argument('--no-clap', action='store_true', help='CIR 只写 DIAM，不写挡光的 CLAP')
    args = ap.parse_args()
    build(args.seq, args)
