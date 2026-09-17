#!/usr/bin/env python3
"""make_seq.py —— spec.json → CODE V 序列文件 (.seq)，与 make_zmx.py 平行的另一条出口。

省掉「Zemax 建好 → 导入 CODE V → 再导出 seq」那一趟往返。
输出格式照 CODE V 2026 自己导出的 .seq 一比一复刻（用户的 Z3512.seq 作参照）。

用法：
    python3 make_seq.py spec.final.json -o Z3512.seq [--glass catalog|exact] [--title 名字]

Zemax → CODE V 的对应关系（都验过）：
  波长      WAVM(升序,主波长 PWAV)  → WL 降序排列 + REF <主波长在降序里的序号> + WTW 同序权重
  视场      FTYP 3 (Real Image Height) → YRI（实际像高），XRI 全 0，WTF 全 1
  渐晕      VDY/VCY/VCX → VUY = VCY − VDY；VLY = VCY + VDY；VUX = VLX = VCX
            （实证 INF 视场1：Zemax VDY −0.1670 / VCY 0.4315 → CODE V VUY 0.5985 / VLY 0.2645）
  口径      DIAM(半口径,固定)      → CIR <半口径>（CODE V 的 CIR 就是实口径，会挡光）
  光阑      STOP                   → STO
  非球面    EVENASPH + PARM2..8    → ASP / K / CUF / A B C D（r⁴ r⁶ r⁸ r¹⁰）/ E F G H（r¹² r¹⁴ r¹⁶ r¹⁸）/ J（r²⁰）
            CODE V 的 K 与 Zemax 的 CONI 同义（都是 conic，不是日系专利的 κ）
  位置解    TOLE <vb> <len> 写在面 va → THI S<va> OAL S<vb>..<va+1> <len>
  多重结构  MNUM/LTTL/THIC        → ZOO n / ZOO TIT + TIT Zk / ZOO THI S<n> + ZOO THC S<n> 100…
  逐结构渐晕 FVCY/FVDY/FVCX       → ZOO VUY F<i> / ZOO VLY F<i> / ZOO VUX F<i> / ZOO VLX F<i>

玻璃：CODE V 的牌号是「去掉所有非字母数字 + _ + 厂家」，如 J-LAK01 → JLAK01_HIKARI。
  --glass catalog（默认）＝ 复刻 CODE V 自己导入时的做法：Offset 解的面直接用**基准目录玻璃**，
     **两个偏移量会丢掉**（实测 EFL 因此从 34.4045 变成 34.4369）。脚本会打印告警。
  --glass exact  ＝ Offset 解的面改写成 CODE V 的模型玻璃 `S <R> <THI> <nd> <vd>`，
     EFL 与专利一致，但丢掉真实玻璃的色散曲线。两者各有取舍，按用途选。
"""
import argparse, json, datetime, re

ASPKEYS = ['A4', 'A6', 'A8', 'A10', 'A12', 'A14', 'A16']
# CODE V 的 ASP 原生带到 r^20（A..D=r^4..r^10, E..H=r^12..r^18, J=r^20），
# 所以 A18/A20 在 .seq 里是**精确**的，不像 Zemax 的 Even Asphere 只到 r^16。
ASPKEYS_FULL = ASPKEYS + ['A18', 'A20']
CVLET = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'J']      # r⁴ … r²⁰


def _apow(A):
    """把非球面字典里的 A<n> 取成 {次数: 系数}，支持奇数次（佳能 A3..A15）。"""
    out = {}
    for k, v in (A or {}).items():
        m = re.fullmatch(r'[Aa]\s*(\d+)', str(k))
        if m and v: out[int(m.group(1))] = float(v)
    return out


def num(v):
    """CODE V 的写法：小数直接写，很小的数写成 0.xxxe±n（照它自己导出的样子）。"""
    v = float(v)
    if v == 0: return '0.0'
    if 1e-4 <= abs(v) < 1e6:
        s = repr(round(v, 12))
        return s if '.' in s or 'e' in s else s + '.0'
    m, e = ('%.10e' % v).split('e')
    m = float(m) / 10.0; e = int(e) + 1
    return '%se%d' % (repr(round(m, 12)), e)


def wrap(tag, vals, width=72):
    """CODE V 的续行：行尾 '&'，续行以一个空格开头。"""
    out, cur = [], tag
    for v in vals:
        if len(cur) + 1 + len(v) > width and cur != tag:
            out.append(cur + '&'); cur = ' ' + v
        else:
            cur += ' ' + v
    out.append(cur + ' ')
    return out


def cvglass(name):
    """'HIKARI J-LAK01' → 'JLAK01_HIKARI'"""
    ven, br = name.split()[0], name.split()[-1]
    return re.sub(r'[^A-Za-z0-9]', '', br).upper() + '_' + ven.upper()


def build(spec, emb, title, gmode):
    zx = spec['zmx']; cfgs = zx['configs']
    st = (emb.get('states') or [None])[0]
    surfs = [s for s in emb['surfaces'] if s['i'] != 'IMG']
    asp = {str(a['surface']).replace('面', ''): a for a in emb.get('aspheric', [])}
    # CODE V 的面号与 Zemax 的面号是同一套（都从物面后第 1 面起顺数，光阑也占一号），
    # 所以 zmx.focus 里的 var_before/var_after 可以直接当 S 号用，不需要再映射。
    L = []; a = L.append
    a('RDM;LEN       "VERSION: 2026       LENS VERSION: 92       Creation Date: %s"'
      % datetime.date.today().strftime('%e-%b-%Y').strip())
    # 结构名的引号写法：CODE V **导出**时写成 TITLE '"INF"' / TIT Z1 ""INF""（双层引号），
    # 但它自己**读**这种写法会读成空 —— 实测把脚本生成的文件在 CODE V 里过一遍，
    # 四个结构名全变成 ' '。所以写入要用单层引号，别照抄它的导出格式。
    a("TITLE '%s'" % title)
    # 孔径：CODE V 的 FNO = **所用共轭下的近轴工作 F 数** 1/(2·NA')（LensSystemSetupRM p.27-29：
    # "adjusts EPD to keep f_ratio satisfied"，NAO=n·sin(atan(EPD/2L))、NA=NAO/RED、FNO=1/(2NA)），
    # 与 Zemax 的 Paraxial Working F/#（FNUM <v> 1）是同一个量。所以两边写同一组逐结构值：
    # 头部 FNO = 结构 1，ZOO FNO 逐结构。值由 vignet.aperture_cfg 给出（专利各态 F 数 / 物理光阑固定）。
    # ★ 旧版把 ∞ 的 F 数原样铺满 ZOO FNO，CODE V 会在近距结构把入瞳**放大**去凑 F2.92
    #   （JP2021-148808A 1:1 结构 EPD 35.7 → 146.7），与 Zemax 旧版的「缩小」正好反向，一样是错的。
    # ★ 有限共轭下两边的定义差一个 sin/tan：CODE V 用物方 NAO = n·sin(atan(EPD/2L))
    #   （LensSystemSetupRM p.29；宏 fct_ABCD.seq：F/# used = |m|·√(r²+L²)/(n·EPD)），
    #   Zemax 的 Paraxial Working F/# 用近轴斜率 tan。要两边**入瞳完全相同**，就得
    #   FNO_CV² = PWFN² + (β/2)²（闭式，与 EPD、L 无关）。RF100 1.4x：6.640 → 6.677（+0.55%），∞ 结构不变。
    from make_zmx import wfno_cfg
    wf0, _ss, bts = wfno_cfg(spec, emb, with_beta=True)
    wf = [(w*w + (b/2.0)**2) ** 0.5 for w, b in zip(wf0, bts)]
    a('FNO   %s' % num(round(wf[0], 6)))
    a('DIM   M')
    # 波长：Zemax 那组固定加权（0.4861/12, 0.5461/30 主, 0.6563/3, 0.5876/22, 0.4358/3）
    wv = zx.get('waves') or [(0.4861, 12), (0.5461, 30), (0.6563, 3), (0.5876, 22), (0.4358, 3)]
    pw = zx.get('pwav', 2)                      # Zemax 的主波长序号（1 起）
    prim = wv[pw - 1][0]
    order = sorted(wv, key=lambda t: -t[0])     # CODE V 按波长降序列
    a('WL    ' + ' '.join('%.1f' % (w * 1000) for w, _ in order))
    a('REF   %d' % (1 + [w for w, _ in order].index(prim)))
    a('WTW   ' + ' '.join(str(int(g)) for _, g in order))
    a("INI   '   '")
    ymax = float(zx['max_y']); nf = 6
    # **CODE V 的视场顺序与 Zemax 相反**：第 1 行是轴上，逐行增大，第 6 行才是最大视场。
    # 用户 2026-09 明确要求。渐晕行（VUY/VLY/VUX/VLX 与 ZOO … F<i>）必须跟着一起倒。
    ys = [round(ymax * i / (nf - 1), 6) for i in range(nf)]
    FLIP = lambda seq: list(seq)[::-1]          # Zemax 是由大到小，CODE V 由小到大
    a('XRI   ' + ' '.join('0.0' for _ in ys))
    a('YRI   ' + ' '.join(num(y) for y in ys))
    a('WTF   ' + ' '.join('1.0' for _ in ys))
    vig = FLIP(zx.get('vignetting') or [[0, 0, 0, 0]] * nf)
    def vuy(f): return round(f[3] - f[1], 6)
    def vly(f): return round(f[3] + f[1], 6)
    for tag, fn in (('VUX', lambda f: round(f[2], 6)), ('VLX', lambda f: round(f[2], 6)),
                    ('VUY', vuy), ('VLY', vly)):
        L.extend(wrap('%-5s' % tag, [num(fn(f)) for f in vig]))
    a('DOR   1.15 1.05 1.1')
    d0 = cfgs[0]['d0']
    a('SO    0.0 %s' % ('0.1e11' if str(d0).upper().startswith('INF') else num(d0)))

    # 变焦镜头：结构里写全了所有可变间隔，逐个进 ZOO THI；不写 OAL 解（各变焦位置守恒和不同）
    zoom = bool(zx.get('zoom'))
    fc, fc2 = (None, None) if zoom else (zx.get('focus'), zx.get('focus2'))
    kb = fc and fc['key_before']; kb2 = fc2 and fc2['key_before']
    warn = []
    varsurf = {}
    _si_last[0] = 0
    for s in surfs:
        sidx = _si(s)                      # 每个面都要过一遍，非整数面号（STO/CG1）才编得对
        R = 0.0 if s['R'] in (None, 0, 'inf') else float(s['R'])
        D = s['D']
        if isinstance(D, str):
            varsurf[D] = sidx
        if isinstance(D, str) and zoom:
            D = float(cfgs[0][D]) if D in cfgs[0] else float(emb['variable'][D][st])
        elif isinstance(D, str):
            D = (cfgs[0].get(D) if D in (kb, kb2) else None) or float(emb['variable'][D][st])
        g = ''
        if s.get('nd'):
            # 面上写了 glass_codev 就原样用它（CODE V 目录里的真牌号）。
            # 两边目录并不一一对应：Zemax 的 ZEON.AGF 把 K22R 与 K26R 合并成
            # ZEONEX_K22R&K26R_2017，而 CODE V 的 ZEON.xml 是分开的 K22R / K26R（同一个
            # 指数 535557），机械地去掉符号拼出来的名字 CODE V 查不到。
            go = s.get('glass_offset')
            if s.get('glass_codev') and not gmode.startswith('exact'):
                g = ' ' + s['glass_codev']
                if go: warn.append((s['i'], s['glass_codev'], go['d_nd'], go['d_vd']))
            elif go and gmode.startswith('exact'):
                g = ' %s %s' % (num(s['nd']), num(s['vd']))
            elif go:
                g = ' ' + cvglass(go['base'])
                warn.append((s['i'], go['base'], go['d_nd'], go['d_vd']))
            elif s.get('glass'):
                g = ' ' + cvglass(s['glass'])
            else:
                g = ' %s %s' % (num(s['nd']), num(s['vd']))     # 模型玻璃
        a('S     %s %s%s' % (num(R), num(D), g))
        phi = (s.get('extra') or {}).get('有効径 φi')
        if phi: a('  CIR %s' % num(round(phi / 2.0, 6)))
        if s.get('stop') or s['i'] == 'STO': a('  STO')
        if s.get('doe'):
            # 衍射面：CODE V DOE，旋转对称相位多项式。HCO Cj = r^(2j) 的 OPD 系数(mm)，
            # 与专利 C_2j 同义（实证 E:/CODEV11.5/lens/bindoub.seq：DIF DOE / HOR / HWL nm / HCT R / HCO）
            d = s['doe']
            a('  DIF DOE'); a('  HOR 1.0')
            a('  HWL %s; HCT R' % num(float(d.get('wl_nm', 587.56))))
            a('  ' + '; '.join('HCO C%d %s' % (j, num(float(v))) for j, v in enumerate(d['C'], 1)))
        A = asp.get(str(s['i']))
        if A and any(n % 2 for n in _apow(A)):
            # 含奇数次项（佳能的 A3..A15）→ CODE V 的 Odd Polynomial 特殊面 SPS ODD。
            # 参数表见 CODE V 2026《Lens System Setup Reference Manual》p.353：
            #   C1 = K（円錐定数）、C2 = AR1(r^1)、C3 = AR2、… C(n+1) = AR n(r^n)，到 r^30；
            #   C33 = NRADIUS。SPS 后面那个数就是 NRADIUS，0.0 = 不归一化（等效 1.0），
            #   与 Zemax 侧 XOSPHERE 把 Rn 取 1.0 一致，两边系数都是专利印的原值。
            # ASP 面只有偶数次（A..J = r^4..r^20），装不下，所以这里不能走 ASP。
            pw = _apow(A)
            a('  SPS ODD 0.0')
            a('  SCO K %s' % num(A.get('k', A.get('K', 0.0))))
            for n in sorted(pw):
                a('  SCO AR%d %s' % (n, num(pw[n])))
        elif A:
            a('  ASP'); a('  K   %s' % num(A.get('k', A.get('K', 0.0)))); a('  CUF 0.0')
            co = [float(A.get(k) or 0.0) for k in ASPKEYS_FULL]   # 9 项 = A..D, E..H, J
            last = max([i for i, v in enumerate(co) if v] or [0])
            # CODE V 的习惯：只要写了 E..H 这一组，后面一定跟一行 J（哪怕是 0）；
            # 全部高次项为零时就只写 A..D 那一行。
            groups = ((0, 1, 2, 3),) if last < 4 else ((0, 1, 2, 3), (4, 5, 6, 7), (8,))
            for grp in groups:
                a('  ' + '; '.join('%s %s' % (CVLET[i], num(co[i])) for i in grp)
                  .replace('A ', 'A   ', 1).replace('E ', 'E   ', 1).replace('J ', 'J   ', 1))
    a('SI    0.0 0.0')
    noaal = '|no-oal' in gmode
    # 对焦间隔逐结构设成变量（ZOO THC 0），与 .zmx 那边 MCE THIC 的 Variable 同一组（make_zmx.focus_vars）。
    from make_zmx import focus_vars
    fvs = set() if '|no-vars' in gmode else {varsurf[k] for k in focus_vars(spec, emb) if k in varsurf}
    solved = []
    # 双浮动对焦群要写**两条** OAL 解。只写第一条时，第二组的 key_after（如 d22）
    # 既没有解也不进 ZOO THI，四个结构里会被一直钉在 ∞ 态的值上 —— 像面跟着跑掉。
    # 链式三段（key_last，如 US11768360B2 的 D14/D20/D22）也要写 OAL：对焦间隔设成变量以后，
    # 没有这条解 CODE V 优化时守恒和不保持、全长漂移（.zmx 那边是 TOLE 14 在面22，同一条约束）。
    for f in (fc, fc2):
        if (not f or noaal or 'var_after' not in f or 'sum' not in f
                or not (f.get('key_after') or f.get('key_last'))): continue
        va, vb = int(f['var_after']), int(f['var_before'])
        solved.append(va)
        tot = f['sum'] + sum(float(x['D']) for x in surfs
                             if isinstance(x['D'], (int, float))
                             and f['var_before'] < _si(x) < f['var_after'])
        a('THI   S%d OAL S%d..%d %s' % (va, vb, va + 1, num(round(tot, 6))))

    # ---- 多重结构 ----
    if len(cfgs) > 1:
        a('ZOO   %d' % len(cfgs)); a('ZOO   TIT')
        for i, c in enumerate(cfgs, 1): a('TIT   Z%d "%s"' % (i, c['name']))
        L.extend(wrap('ZOO   FNO', [num(round(v, 6)) for v in wf]))
        vc = [FLIP(v) for v in (zx.get('vignetting_cfg') or [])]
        if vc and len(vc) == len(cfgs):
            vx = lambda f: round(f[2], 6)
            for pair in ((('VUY', vuy), ('VLY', vly)), (('VUX', vx), ('VLX', vx))):
                for fi in range(nf):
                    for tag, fn in pair:
                        L.extend(wrap('ZOO   %s F%d' % (tag, fi + 1),
                                      [num(fn(v[fi])) for v in vc]))
        # ⚠ **带 OAL 解的那一面不能再写 ZOO THI** —— 文件是按顺序执行的命令流，
        # `ZOO THI S26 …` 排在 `THI S26 OAL …` 后面，会把解直接覆盖成一组定值，
        # 于是"求解不起作用"。CODE V 自己**导出**时两行都写（它只是在转储当前值），
        # 但那份导出并不是一份合法的输入 —— 别照抄。
        # 这条规则与 zmx 侧同源：用了位置解，MCE 里就不要再写该面的 THIC 行。
        rows = [(0, [c['d0'] for c in cfgs])]
        if zoom:
            for key, sn in sorted(varsurf.items(), key=lambda t: t[1]):
                vals = [float(c[key]) if key in c else float(emb['variable'][key][st]) for c in cfgs]
                if not all(abs(v - vals[0]) < 1e-9 for v in vals):
                    rows.append((sn, vals))
        if kb and fc: rows.append((int(fc['var_before']), [c[kb] for c in cfgs]))
        if fc and fc.get('key_last') and not solved:
            # 链式三段浮动：位置解那一面的厚度 = 守恒和 − 前两个可变间隔
            rows.append((int(fc['var_after']),
                         [round(fc['sum'] - c[kb] - c[kb2], 6) for c in cfgs]))
        elif fc and fc.get('key_after') and int(fc['var_after']) not in solved:
            rows.append((int(fc['var_after']),
                         [round(fc['sum'] - c[kb], 6) for c in cfgs]))
        if kb2 and fc2: rows.append((int(fc2['var_before']), [c[kb2] for c in cfgs]))
        if fc2 and fc2.get('key_after') and fc2.get('var_after') and int(fc2['var_after']) not in solved:
            # --no-oal 时第二对焦群的 key_after 也要逐结构给值（旧版漏了：RF100 的 D29 一直停在 ∞ 值）
            rows.append((int(fc2['var_after']), [round(fc2['sum'] - c[kb2], 6) for c in cfgs]))
        rows = [r for r in rows if r[0] not in solved]
        for sn, vals in rows:
            v = ['0.1e11' if str(x).upper().startswith('INF') else num(x) for x in vals]
            L.extend(wrap('ZOO   THI S%d' % sn, v))
            a('ZOO THC S%d %s' % (sn, ' '.join(('0' if sn in fvs else '100') for _ in vals)))
    a('GO')
    return L, warn


_si_last = [0]
def _si(s):
    try: _si_last[0] = int(s['i'])
    except (TypeError, ValueError): _si_last[0] += 1
    return _si_last[0]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('spec'); ap.add_argument('-o', required=True)
    ap.add_argument('--emb', type=int, default=0)
    ap.add_argument('--glass', choices=('catalog', 'exact'), default='catalog')
    ap.add_argument('--no-vars', action='store_true',
                    help='对焦间隔不设变量（默认 ZOO THC 0 = 变量，打开就能优化对焦位置）')
    ap.add_argument('--no-oal', action='store_true',
                    help='不写 THI S<va> OAL 位置解。CODE V 里对焦组位置一般是用**评价函数优化**出来的，'
                         '位置解会占掉一个自由度；想让三个间隔都自由参与优化时用这个开关。'
                         '默认写解（它编码的是「G1 与像面固定 → 三段间隔之和恒定」这条物理约束，'
                         '留着更省事，也不妨碍把 d18/d22 设成变量）。')
    ap.add_argument('--title')
    a = ap.parse_args()
    spec = json.load(open(a.spec, encoding='utf-8')); emb = spec['embodiments'][a.emb]
    title = a.title or (spec['zmx'].get('configs') or [{'name': 'LENS'}])[0]['name']
    L, warn = build(spec, emb, title, a.glass + ('|no-oal' if a.no_oal else '')
                    + ('|no-vars' if a.no_vars else ''))
    open(a.o, 'wb').write(('\n'.join(L) + '\n').encode('latin-1'))
    print('已写出 %s（%d 行）' % (a.o, len(L)))
    if warn:
        print('  ⚠ 下列面在 zmx 里是 Offset 玻璃解，CODE V 没有对应功能，本文件按 --glass catalog '
              '只写了基准牌号，**Nd/Vd 偏移丢失**：')
        for i, b, dn, dv in warn:
            print('     面%-3s 基准 %-20s 丢掉 Nd%+.5f / Vd%+.4f' % (i, b, dn, dv))
        print('  → 想让 EFL 与专利一致就改用 --glass exact（该面写成模型玻璃 nd/vd）。')


if __name__ == '__main__':
    main()
