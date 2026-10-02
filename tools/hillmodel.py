"""Model of surfaceTransform hill ridges (verified on our render): the drawing turtle integrates
f(i) = a*i + b + sine[((sa*i) >> 8) + sb] * sc  (sine table of 256 entries, analysis/sine_table.json)
in 1/4096 px and moves by a whole pixel only when the exact value is a full pixel away (dead band)."""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from polyfit import line  # noqa

T = json.load(open(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'analysis', 'sine_table.json')))


def ridge(y0, x0, n, p=(0, 0), s=(0, 0, 0)):
    a, b = p
    S = 0; kk = 0; out = []
    for i in range(n):
        q = S / 4096.0
        if q - kk >= 1:
            kk = int(q // 1)
        elif kk - q > 1:
            kk = -int((-q) // 1)
        out.append((x0 + i, y0 - kk))
        S += a * i + b + T[(((s[0] * i) >> 8) + s[1]) & 255] * s[2]
    return out


def tops(pts, W=600):
    """topmost drawn pixel per column of the polyline"""
    top = {}
    for p0, p1 in zip(pts, pts[1:]):
        for x, y in line(p0, p1):
            if 0 <= x < W and (x not in top or y < top[x]):
                top[x] = y
    return top
