#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""vigfit_merge.py —— 把 zapi_vigfit.ps1（OpticStudio 真追迹微调）的逐结构渐晕写回 spec。

    python vigfit_merge.py spec.final.json X.vigfit.json [-o spec.final.json]

zapi_vigfit 输出的是 Zemax 视场顺序（最大视场在前，与 vignet.py 的 FRACS 同序）、
每项 [VDX, VDY, VCX, VCY] —— 与 zmx.vignetting_cfg 的 [vdx, vdy, vcx, vcy] 同一排布，直接覆盖。
覆盖后 make_zmx.py / make_seq.py 两条出口用的就是同一套被 OpticStudio 验过的渐晕。
"""
import json, argparse, io


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('spec'); ap.add_argument('fit'); ap.add_argument('-o')
    a = ap.parse_args()
    spec = json.load(io.open(a.spec, encoding='utf-8'))
    raw = open(a.fit, 'rb').read()
    fit = json.loads(raw.decode('utf-8-sig'))          # Out-File -Encoding utf8 带 BOM
    zx = spec['zmx']
    old = zx.get('vignetting_cfg') or []
    new = [[[float(v) for v in f] for f in cfg] for cfg in fit['vig']]
    if old and (len(old) != len(new) or len(old[0]) != len(new[0])):
        raise SystemExit('结构数/视场数对不上：spec %dx%d vs vigfit %dx%d'
                         % (len(old), len(old[0]), len(new), len(new[0])))
    names = [c.get('name') for c in zx.get('configs') or []] or ['cfg%d' % (i+1) for i in range(len(new))]
    for ci, (cn, cfg) in enumerate(zip(names, new)):
        for fi, v in enumerate(cfg):
            o = old[ci][fi] if old else [0, 0, 0, 0]
            d = [round(v[k] - o[k], 4) for k in range(4)]
            if any(abs(x) > 1e-4 for x in d):
                print('  %-12s 视场%d  VDY %+.4f→%+.4f  VCX %.4f→%.4f  VCY %.4f→%.4f'
                      % (cn, fi+1, o[1], v[1], o[2], v[2], o[3], v[3]))
    zx['vignetting_cfg'] = new
    zx['vignetting'] = new[0]
    zx['vignetting_source'] = 'OpticStudio ZOS-API 真追迹微调（zapi_vigfit.ps1, step %s）' % fit.get('step')
    out = a.o or a.spec
    json.dump(spec, io.open(out, 'w', encoding='utf-8'), ensure_ascii=False, indent=1)
    print('已写回 %s' % out)


if __name__ == '__main__':
    main()
