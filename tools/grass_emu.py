"""Emulator of drawGrassPatch: RC4 (fastRandom) keyed by a phrase; per tuft x=2*b1, y=b2, type=b3&3 (0 = none).
Sprites are cut from our own render (out/b_grass.png vs out/b_flw.png)."""
import numpy as np
from PIL import Image

ROOT = __import__('os').path.dirname(__import__('os').path.dirname(__import__('os').path.abspath(__file__)))


def codes(s):
    return [ch.encode('cp037')[0] - 64 for ch in s]


def rc4(k, n):
    S = list(range(256)); j = 0; L = len(k)
    for i in range(256):
        j = (j + S[i] + k[i % L]) % 256; S[i], S[j] = S[j], S[i]
    i = j = 0; out = []
    for _ in range(n):
        i = (i + 1) % 256; j = (j + S[i]) % 256; S[i], S[j] = S[j], S[i]; out.append(S[(S[i] + S[j]) % 256])
    return out


OX, OY = 30, 320


def tufts(key, count=70):
    """tuft origins (absolute) and types; sprites from analysis/grass_sprites.npy are relative to the origin"""
    k = codes(key) if isinstance(key, str) else key
    o = rc4(k, 3 * count)
    return [(2 * o[3 * i] + OX, o[3 * i + 1] + OY, o[3 * i + 2] & 3) for i in range(count) if o[3 * i + 2] & 3]


def clean_sprites():
    S = np.load(f'{ROOT}/analysis/grass_sprites.npy', allow_pickle=True)
    return {i + 1: tuple(S[i]) for i in range(4)}


WX0, WX1, WY0, WY1 = 20, 70, 305, 350


def load(base='b_flw', withp='b_grass'):
    n = np.array(Image.open(f'{ROOT}/out/{base}.png').convert('RGB')).astype(int)
    g = np.array(Image.open(f'{ROOT}/out/{withp}.png').convert('RGB')).astype(int)
    return n, g


def sprites(n, g, key="El pasto siempre se ve mas verde del otro lado de la cerca."):
    d = (g != n).any(axis=2)
    spr = {}
    for T in (1, 2, 3):
        best = None
        for x, y, t in tufts(key):
            if t != T or x + WX1 > 600 or y + WY1 > 600:
                continue
            c = d[y + WY0:y + WY1, x + WX0:x + WX1].sum()
            others = sum(1 for (a, b, _) in tufts(key) if (a, b) != (x, y) and abs(a - x) < 30 and abs(b - y) < 25)
            if others == 0 and (best is None or c > best[0]):
                best = (c, x, y)
        _, x, y = best
        m = d[y + WY0:y + WY1, x + WX0:x + WX1]
        ys, xs = np.nonzero(m)
        spr[T] = (ys + WY0, xs + WX0, g[y + WY0:y + WY1, x + WX0:x + WX1][ys, xs])
    return spr


def lawn(n):
    return (n[:, :, 0] <= 3) & (n[:, :, 2] <= 12) & (n[:, :, 1] >= 40)


def render(n, spr, tl):
    im = n.copy(); L = lawn(n)
    for x, y, t in tl:
        ys, xs, col = spr[t]
        Y = ys + y; X = xs + x
        ok = (Y >= 0) & (Y < 600) & (X >= 0) & (X < 600)
        Y, X, c = Y[ok], X[ok], col[ok]
        ok2 = L[Y, X]
        im[Y[ok2], X[ok2]] = c[ok2]
    return im
