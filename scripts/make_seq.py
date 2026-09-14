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
CVLET = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'J']      # r⁴ … r²⁰


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
    a('FNO   %s' % num(zx['fno']))
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

    fc, fc2 = zx.get('focus'), zx.get('focus2')
    kb = fc and fc['key_before']; kb2 = fc2 and fc2['key_before']
    warn = []
    for s in surfs:
        R = 0.0 if s['R'] in (None, 0, 'inf') else float(s['R'])
        D = s['D']
        if isinstance(D, str):
            D = (cfgs[0].get(D) if D in (kb, kb2) else None) or float(emb['variable'][D][st])
        g = ''
        if s.get('nd'):
            go = s.get('glass_offset')
            if go and gmode.startswith('exact'):
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
        A = asp.get(str(s['i']))
        if A:
            a('  ASP'); a('  K   %s' % num(A.get('k', 0.0))); a('  CUF 0.0')
            co = [float(A.get(k) or 0.0) for k in ASPKEYS] + [0.0, 0.0]   # 补到 9 项(H、J)
            last = max([i for i, v in enumerate(co) if v] or [0])
            # CODE V 的习惯：只要写了 E..H 这一组，后面一定跟一行 J（哪怕是 0）；
            # 全部高次项为零时就只写 A..D 那一行。
            groups = ((0, 1, 2, 3),) if last < 4 else ((0, 1, 2, 3), (4, 5, 6, 7), (8,))
            for grp in groups:
                a('  ' + '; '.join('%s %s' % (CVLET[i], num(co[i])) for i in grp)
                  .replace('A ', 'A   ', 1).replace('E ', 'E   ', 1).replace('J ', 'J   ', 1))
    a('SI    0.0 0.0')
    noaal = gmode.endswith('|no-oal')
    solved = set()
    if fc and not noaal:
        va, vb = int(fc['var_after']), int(fc['var_before'])
        solved.add(va)
        tot = fc['sum'] + sum(float(x['D']) for x in surfs
                              if isinstance(x['D'], (int, float))
                              and fc['var_before'] < _si(x) < fc['var_after'])
        a('THI   S%d OAL S%d..%d %s' % (va, vb, va + 1, num(round(tot, 6))))
        # 双浮动对焦（两对各自守恒，如 WO2024147268 尼康 135/1.8）：**第 2 组也要有自己的
        # OAL 解**。少写这一条的话 ka2 那个间隔在所有结构里被冻结成基准值 —— 第 2 对焦群
        # 只前进不后退，全长跟着变，像面在近距结构上整体跑掉。
        # 链式三段（key_last）只有一个守恒和，不能再加第二条。
        if fc2 and fc2.get('key_after') and not fc.get('key_last'):
            va2, vb2 = int(fc2['var_after']), int(fc2['var_before'])
            solved.add(va2)
            tot2 = fc2['sum'] + sum(float(x['D']) for x in surfs
                                    if isinstance(x['D'], (int, float))
                                    and fc2['var_before'] < _si(x) < fc2['var_after'])
            a('THI   S%d OAL S%d..%d %s' % (va2, vb2, va2 + 1, num(round(tot2, 6))))

    # ---- 多重结构 ----
    if len(cfgs) > 1:
        a('ZOO   %d' % len(cfgs)); a('ZOO   TIT')
        for i, c in enumerate(cfgs, 1): a('TIT   Z%d "%s"' % (i, c['name']))
        a('ZOO   FNO ' + ' '.join(num(zx['fno']) for _ in cfgs))
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
        if kb and fc: rows.append((int(fc['var_before']), [c[kb] for c in cfgs]))
        if fc and fc.get('key_last') and not solved:
            # 链式三段浮动：位置解那一面的厚度 = 守恒和 − 前两个可变间隔
            rows.append((int(fc['var_after']),
                         [round(fc['sum'] - c[kb] - c[kb2], 6) for c in cfgs]))
        elif fc and fc.get('key_after') and not solved:
            rows.append((int(fc['var_after']),
                         [round(fc['sum'] - c[kb], 6) for c in cfgs]))
        if kb2 and fc2: rows.append((int(fc2['var_before']), [c[kb2] for c in cfgs]))
        if fc2 and fc2.get('key_after') and not solved:
            rows.append((int(fc2['var_after']),
                         [round(fc2['sum'] - c[kb2], 6) for c in cfgs]))
        rows = [r for r in rows if r[0] not in solved]
        for sn, vals in rows:
            v = ['0.1e11' if str(x).upper().startswith('INF') else num(x) for x in vals]
            L.extend(wrap('ZOO   THI S%d' % sn, v))
            a('ZOO THC S%d %s' % (sn, ' '.join('100' for _ in vals)))
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
    ap.add_argument('--no-oal', action='store_true',
                    help='不写 THI S<va> OAL 位置解。CODE V 里对焦组位置一般是用**评价函数优化**出来的，'
                         '位置解会占掉一个自由度；想让三个间隔都自由参与优化时用这个开关。'
                         '默认写解（它编码的是「G1 与像面固定 → 三段间隔之和恒定」这条物理约束，'
                         '留着更省事，也不妨碍把 d18/d22 设成变量）。')
    ap.add_argument('--title')
    a = ap.parse_args()
    spec = json.load(open(a.spec, encoding='utf-8')); emb = spec['embodiments'][a.emb]
    title = a.title or (spec['zmx'].get('configs') or [{'name': 'LENS'}])[0]['name']
    L, warn = build(spec, emb, title, a.glass + ('|no-oal' if a.no_oal else ''))
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
