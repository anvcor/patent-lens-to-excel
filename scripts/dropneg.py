# -*- coding: utf-8 -*-
"""删掉「负厚度的无光焦度面」（遮光平面 / FS / 虚拟基准面），把负厚度并进前一个间隔，其后面号 −1。

    python dropneg.py spec.json [-o out.json] [--emb N] [--dry]

用户 2026-09-23 要求交付前去掉（JP2020-118807A RF 28-70 F2 L 的面22 d=−2.59、WO2023181666A1 24-70 GM II 的面21 d=−1.10）。
判定条件（全部满足才删）：R 为平面（None/0/inf）、不是玻璃（没有 nd）、不是光阑、不是盖板、厚度 < 0
（数值负，或可变间隔在所有状态都 < 0）。光学完全不变：它没有光焦度，前后两个间隔相加 = 原来的总间隔。

并入规则（设被删面为 k，前一面为 k−1）：
  - 前间隔与 k 的厚度都是数值 → 前间隔 += d_k
  - 有一边是可变间隔 → 合成可变间隔，**沿用前间隔的名字**（前间隔是数值时沿用 k 的名字、改名为 D<k-1>），
    各状态值 = 两边逐状态相加；合出来有 ≤0 的状态会报错（说明遮光面后面的间隔撑不住，不能直接删）。
改号：面号 > k 的全部 −1；可变间隔名 `D<n>`（大小写均可）n>k 的改成 D<n−1>，
  同步改 variable / d0 / zmx.configs / zmx.focus_vars / zmx.focus* 的 key_* 与 var_* / zmx.mce_extra /
  zmx.fix_semi_surfaces / zmx.vig_def_surfaces / 非球面「面N」/ 逐面列表 semi_3d / axial_3d（删掉该面那一项）。
被删面的口径（遮光作用）随之消失 —— 下游 vignet / apcap / zapi_vigfit 要重跑。

建议位置：mkspec.py 生成 spec 之后、lensmath.py 之前（此时下游字段还没有，最干净）。
"""
import argparse, io, json, re

VRE = re.compile(r'^([Dd])(\d+)$')


def _flat(R):
    return R is None or R == 0 or (isinstance(R, str) and R.lower().startswith('inf'))


def _ren(name, k):
    """可变间隔名改号：D<n>，n>k → D<n-1>"""
    if not isinstance(name, str):
        return name
    m = VRE.match(name)
    if m and int(m.group(2)) > k:
        return '%s%d' % (m.group(1), int(m.group(2)) - 1)
    return name


def _renum(i, k):
    return i - 1 if isinstance(i, int) and i > k else i


def candidates(emb):
    var = emb.get('variable') or {}
    out = []
    for n, s in enumerate(emb['surfaces']):
        if s['i'] == 'IMG' or n == 0:
            continue
        if not _flat(s.get('R')) or s.get('nd') or s.get('stop') or s.get('cover'):
            continue
        D = s['D']
        if isinstance(D, str):
            vs = list((var.get(D) or {}).values())
            neg = bool(vs) and all(float(v) < 0 for v in vs)
        else:
            neg = float(D) < 0
        if neg:
            out.append(n)
    return out


def drop_one(spec, emb, n):
    """删 emb['surfaces'][n]，返回说明文字。"""
    surfs = emb['surfaces']; var = emb.setdefault('variable', {})
    s, p = surfs[n], surfs[n - 1]
    k = s['i']
    if not isinstance(k, int):
        raise SystemExit('面号 %r 不是整数，脚本不处理' % k)
    dk, dp = s['D'], p['D']
    if not isinstance(dk, str) and not isinstance(dp, str):
        newD = round(float(dp) + float(dk), 6)
        if newD <= 0:
            raise SystemExit('面%s 并入后间隔 %.4f ≤ 0，不能直接删' % (p['i'], newD))
        p['D'] = newD
        msg = '面%s（d=%s）并入面%s：%s → %s' % (k, dk, p['i'], dp, newD)
    else:
        states = list(var[dk].keys()) if isinstance(dk, str) else list(var[dp].keys())
        a = var[dp] if isinstance(dp, str) else {st: float(dp) for st in states}
        b = var[dk] if isinstance(dk, str) else {st: float(dk) for st in states}
        merged = {st: round(float(a[st]) + float(b[st]), 6) for st in states}
        bad = [st for st, v in merged.items() if v <= 0]
        if bad:
            raise SystemExit('面%s 并入后可变间隔在 %s 状态 ≤ 0，不能直接删' % (p['i'], bad))
        name = dp if isinstance(dp, str) else 'D%d' % p['i']
        if isinstance(dk, str) and dk != name:
            var.pop(dk, None)
        var[name] = merged
        p['D'] = name
        msg = '面%s（d=%s）并入面%s 的可变间隔 %s：%s' % (k, dk, p['i'], name,
                                                   ' / '.join('%g' % v for v in merged.values()))
        # zmx.configs（lensmath 之后才有）里逐结构的间隔值也要相加；被删面的可变间隔名并掉
        for c in (spec.get('zmx') or {}).get('configs') or []:
            pv = c.get(dp) if isinstance(dp, str) else float(dp)
            kv = c.get(dk) if isinstance(dk, str) else float(dk)
            if pv is not None and kv is not None:
                c[name] = round(float(pv) + float(kv), 6)
            if isinstance(dk, str) and dk != name:
                c.pop(dk, None)
        if isinstance(dk, str) and dk != name:
            _rename_key(spec, emb, dk, name)
    note = '专利面%s（无光焦度平面，d=%s）已删，厚度并入本间隔' % (k, dk)
    p['note'] = (p.get('note') + '；' if p.get('note') else '') + note
    del surfs[n]
    # ---- 改号 ----
    for q in surfs:
        q['i'] = _renum(q['i'], k)
        q['D'] = _ren(q['D'], k)
    ren = {v: _ren(v, k) for v in list(var)}
    emb['variable'] = {ren[v]: var[v] for v in var}
    for a in emb.get('aspheric', []):
        m = re.match(r'^面(\d+)$', str(a.get('surface', '')))
        if m and int(m.group(1)) > k:
            a['surface'] = '面%d' % (int(m.group(1)) - 1)
        elif m and int(m.group(1)) == k:
            raise SystemExit('被删面%s 带非球面系数，不该删' % k)
    zx = spec.get('zmx') or {}
    for c in zx.get('configs') or []:
        for key in list(c):
            nk = _ren(key, k)
            if nk != key:
                c[nk] = c.pop(key)
    if zx.get('focus_vars') is not None:
        zx['focus_vars'] = [_ren(v, k) for v in zx['focus_vars']]
    for fk in ('focus', 'focus2'):
        f = zx.get(fk)
        if isinstance(f, dict):
            for key in list(f):
                if key.startswith('key_'):
                    f[key] = _ren(f[key], k)
                elif key.startswith('var_') and isinstance(f[key], int):
                    f[key] = _renum(f[key], k)
    for ex in zx.get('mce_extra') or []:
        ex['var'] = _renum(ex.get('var'), k)
        ex['key'] = _ren(ex.get('key'), k)
    for key in ('fix_semi_surfaces', 'vig_def_surfaces'):
        if zx.get(key) is not None:
            zx[key] = [_renum(v, k) for v in zx[key] if v != k]
    for key in ('semi_3d', 'axial_3d'):
        if isinstance(zx.get(key), list) and len(zx[key]) > n:
            del zx[key][n]
    return msg


def _rename_key(spec, emb, old, new):
    zx = spec.get('zmx') or {}
    for c in zx.get('configs') or []:
        if old in c:
            c[new] = c.pop(old)
    if zx.get('focus_vars'):
        zx['focus_vars'] = [new if v == old else v for v in zx['focus_vars']]


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('spec')
    ap.add_argument('-o', help='输出（默认覆盖原文件）')
    ap.add_argument('--emb', type=int, default=0)
    ap.add_argument('--dry', action='store_true', help='只列出会删的面，不写文件')
    a = ap.parse_args()
    spec = json.load(io.open(a.spec, encoding='utf-8'))
    emb = spec['embodiments'][a.emb]
    idx = candidates(emb)
    if not idx:
        print('没有负厚度的无光焦度面')
        return
    print('会删：', ', '.join('面%s（d=%s）' % (emb['surfaces'][n]['i'], emb['surfaces'][n]['D']) for n in idx))
    if a.dry:
        return
    for n in sorted(idx, reverse=True):          # 从后往前删，前面的下标不受影响
        print(' ', drop_one(spec, emb, n))
    zx = spec.get('zmx') or {}
    if any(zx.get(k) for k in ('semi_3d', 'vignetting_cfg')):
        print('  ⚠ spec 里已有渐晕/口径结果：被删面的遮光作用没了，vignet → clearance → apcap → zapi_vigfit 要重跑')
    out = a.o or a.spec
    json.dump(spec, io.open(out, 'w', encoding='utf-8'), ensure_ascii=False, indent=1)
    print('已写出', out, '（%d 面）' % len([s for s in emb['surfaces'] if s['i'] != 'IMG']))


if __name__ == '__main__':
    main()
