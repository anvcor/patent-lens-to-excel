#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""refocus_merge.py —— 把 zapi_refocus.ps1（轴上离焦 MTF 峰值居中）的逐结构对焦间隔写回 spec。

    python refocus_merge.py spec.final.json X.refocus.json [-o spec.final.json]

zapi_refocus 给的是 MCE 里对焦间隔那一面（focus.var_before 或 focus2.var_before）的新值，
按结构名对回 zmx.configs，覆盖该结构的 key_before；key_after / key_last 由位置解（TOLE / OAL）跟随，
守恒和不变。覆盖后 make_zmx.py / make_seq.py 两条出口用的是同一组重新对焦后的间隔。
"""
import json, argparse, io


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('spec'); ap.add_argument('fit'); ap.add_argument('-o')
    a = ap.parse_args()
    spec = json.load(io.open(a.spec, encoding='utf-8'))
    fit = json.loads(open(a.fit, 'rb').read().decode('utf-8-sig'))
    zx = spec['zmx']
    sn = int(fit['surface'])
    key = None
    for f in (zx.get('focus'), zx.get('focus2')):
        if f and int(f.get('var_before', -1)) == sn:
            key = f['key_before']
    if key is None:
        raise SystemExit('面%d 不是 focus / focus2 的 var_before —— 不知道该写哪个间隔' % sn)
    cfgs = zx.get('configs') or []
    byname = {c['name']: c for c in cfgs}
    for r in fit['configs']:
        c = byname.get(r['name'])
        if c is None:
            ci = int(r['config']) - 1
            c = cfgs[ci] if ci < len(cfgs) else None
        if c is None:
            raise SystemExit('结构 %s 在 spec 里找不到' % r['name'])
        print('  %-14s %s %.4f → %.4f  (%+.4f)   轴上 %g lp/mm MTF %.3f → %.3f   峰位 %+.4f → %+.4f mm'
              % (c['name'], key, float(c.get(key, r['old'])), r['new'], r['new'] - r['old'],
                 fit.get('freq', 50), r['mtf_old'], r['mtf_new'], r['peak_old'], r['peak_new']))
        c[key] = round(float(r['new']), 4)
    zx['refocus_source'] = ('OpticStudio ZOS-API 轴上离焦 MTF 峰值居中（zapi_refocus.ps1, %g lp/mm, 多色, 面%d）'
                            % (fit.get('freq', 50), sn))
    out = a.o or a.spec
    json.dump(spec, io.open(out, 'w', encoding='utf-8'), ensure_ascii=False, indent=1)
    print('已写回 %s' % out)


if __name__ == '__main__':
    main()
