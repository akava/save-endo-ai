"""Exact emulator of `weeds` (G+6608208): 9 stochastic L-system weeds (lsystem-weed-F) driven by randomInt.
randomInt(n): seed = (seed mod 3513381) * 3067221 mod 2^22, returns seed mod n (verified on all 141 draws of the
original scene). Turtle: 8-bit fixed point, sineArray from the green zone, pixel = floor."""
import json
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SIN = json.load(open(os.path.join(ROOT, 'analysis', 'sine_table.json')))
C1, C2 = 3513381, 3067221
RULES = ['F[+F]F[-F]F', 'F[+F]F', 'F[-F]F']


class Rng:
    def __init__(self, seed):
        self.s = seed

    def __call__(self, n):
        self.s = (self.s % C1) * C2 % (1 << 22)
        return self.s % n


def weed(rng, x, y, depth=3, dist=5, ang=15, head=192):
    segs = []

    def F(d, st):
        if d == 0:
            nx, ny = st[0] + dist * SIN[(st[2] + 64) % 256], st[1] + dist * SIN[st[2]]
            segs.append(((st[0] >> 8, st[1] >> 8), (nx >> 8, ny >> 8)))
            return (nx, ny, st[2])
        stack = []
        for c in RULES[rng(3)]:
            if c == 'F':
                st = F(d - 1, st)
            elif c == '[':
                stack.append(st)
            elif c == ']':
                st = stack.pop()
            elif c == '+':
                st = (st[0], st[1], (st[2] - ang) % 256)
            else:
                st = (st[0], st[1], (st[2] + ang) % 256)
        return st
    F(depth, (x * 256, y * 256, head))
    return segs


def weeds(seed=4237, depth=3, origin=(20, 250), dist=5, ang=15, n=9):
    rng = Rng(seed)
    return [weed(rng, origin[0] + 12 * k, origin[1] + (10 if k % 2 else 0), depth, dist, ang) for k in range(n)]


def raster(ws, W=600):
    """Bresenham-ish rasterization of all segments: set of (x, y)."""
    pts = set()
    for w in ws:
        for (x0, y0), (x1, y1) in w:
            n = max(abs(x1 - x0), abs(y1 - y0), 1)
            for i in range(n + 1):
                pts.add((x0 + round((x1 - x0) * i / n), y0 + round((y1 - y0) * i / n)))
    return pts
