"""Ridge polylines of the parabola/sine hills from an RNA file: {bucket-signature: [(x, y), ...]}"""
import sys
from collections import defaultdict
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rnalines import lines  # noqa


def ridges(rna):
    byb = defaultdict(list)
    for q in lines(rna):
        byb[q[3]].append(q)
    out = {}
    for b, v in byb.items():
        if len(v) > 200 and all(abs(q[2][0] - q[1][0]) <= 4 for q in v[:50]) and v[0][1][1] > 150:
            pts = [v[0][1]] + [q[2] for q in v]
            out[''.join(sorted(b))] = pts
    return out


if __name__ == '__main__':
    for b, p in ridges(sys.argv[1]).items():
        print(b[:20], len(p), p[:12])
