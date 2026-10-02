"""Coordinate-descent re-tuning of build_best.PARAMS with the corrected metric. usage: retune.py key [key...]"""
import itertools
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import build_best as bb  # noqa
from batch import batch  # noqa


def tune(key, steps):
    base = bb.PARAMS[key]
    vals = base if isinstance(base, tuple) else (base,)
    cands = set()
    for idx, st in enumerate(steps):
        for d in st:
            v = list(vals); v[idx] += d; cands.add(tuple(v))
    cands.add(tuple(vals))
    t = {}
    for c in cands:
        bb.PARAMS[key] = c if isinstance(base, tuple) else c[0]
        t['rt_%s_%s' % (key, '_'.join(map(str, c)))] = bb.build()
    bb.PARAMS[key] = base
    r = batch(t, None, timeout=900, jobs=4)
    r = sorted(r, key=lambda q: (q[2] if q[2] is not None else 10**9, q[1]))
    best = r[0][0].split('_')[2:]
    v = tuple(int(x) for x in best)
    bb.PARAMS[key] = v if isinstance(base, tuple) else v[0]
    print('BEST', key, bb.PARAMS[key], r[0][2], flush=True)
    return r[0][2]


STEPS = {2: [(-4, -2, 2, 4), (-4, -2, 2, 4)], 3: [(-2, 2), (-4, -2, 2, 4), (-1, 1)], 1: [(-2, -1, 1, 2)]}
if __name__ == '__main__':
    for key in sys.argv[1:]:
        base = bb.PARAMS[key]
        n = len(base) if isinstance(base, tuple) else 1
        tune(key, STEPS[n])
    print('PARAMS =', bb.PARAMS)
