"""Helpers: encode DNA prefixes, run executor+builder, score against target."""
import os
import subprocess
import sys

import numpy as np
from PIL import Image

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DNA = os.path.join(ROOT, 'build', 'dna')
BUILD = os.path.join(ROOT, 'build', 'build')
OUT = os.path.join(ROOT, 'out')

# ---------- encoding ----------
BASE = {'I': 'C', 'C': 'F', 'F': 'P', 'P': 'IC'}  # literal base in pattern/template


def nat(n):
    s = ''
    while n:
        s += 'C' if n & 1 else 'I'
        n >>= 1
    return s + 'P'


def lit(s):
    return ''.join(BASE[c] for c in s)


class P:
    """pattern builder"""
    def __init__(self):
        self.s = ''

    def b(self, bases):
        self.s += lit(bases); return self

    def skip(self, n):
        self.s += 'IP' + nat(n); return self

    def search(self, bases):
        self.s += 'IFF' + lit(bases); return self

    def open(self):
        self.s += 'IIP'; return self

    def close(self):
        self.s += 'IIC'; return self

    def end(self):
        return self.s + 'IIC'


class T:
    """template builder"""
    def __init__(self):
        self.s = ''

    def b(self, bases):
        self.s += lit(bases); return self

    def ref(self, n, l=0):
        self.s += 'IP' + nat(l) + nat(n); return self

    def len(self, n):
        self.s += 'IIP' + nat(n); return self

    def end(self):
        return self.s + 'IIC'


def rna_cmd(cmd7):
    """prefix fragment that emits one RNA command"""
    assert len(cmd7) == 7
    return 'III' + cmd7


def set_base(pos, base, old=None):
    """replace original base at absolute position pos with `base`"""
    p = P().open().skip(pos).close()
    if old:
        p.b(old)
    else:
        p.skip(1)
    return p.end() + T().ref(0).b(base).end()


def write_at(pos, bases, oldlen):
    """replace oldlen bases at pos by `bases`"""
    p = P().open().skip(pos).close().skip(oldlen).end()
    return p + T().ref(0).b(bases).end()


# ---------- running ----------
def run(prefix, name='run', layers=False, extra=()):
    os.makedirs(OUT, exist_ok=True)
    pf = os.path.join(OUT, name + '.prefix')
    with open(pf, 'w') as f:
        f.write(prefix)
    rna = os.path.join(OUT, name + '.rna')
    png = os.path.join(OUT, name + '.png')
    r = subprocess.run([DNA, '-d', os.path.join(ROOT, 'data', 'endo.dna'), '-p', pf, '-o', rna, *extra], capture_output=True, text=True)
    info = r.stderr.strip()
    cmd = [BUILD, rna, png]
    if layers:
        cmd += ['--layers', os.path.join(OUT, name + '_layer')]
    r2 = subprocess.run(cmd, capture_output=True, text=True)
    return png, info + ' | ' + r2.stderr.strip()


_target = None


def target():
    global _target
    if _target is None:
        _target = np.array(Image.open(os.path.join(ROOT, 'doc', 'Target-image.png')).convert('RGB')).astype(int)
    return _target


def score(png, tol=24):
    """approximate number of wrong 600x600 pixels: compare 2x2-averaged render with the 300x300 target"""
    img = np.array(Image.open(png).convert('RGB')).astype(int)
    small = img.reshape(300, 2, 300, 2, 3).mean(axis=(1, 3))
    d = np.abs(small - target()).max(axis=2)
    return int((d > tol).sum() * 4)


if __name__ == '__main__':
    pre = sys.argv[1] if len(sys.argv) > 1 else ''
    png, info = run(pre)
    print(info)
    print('approx wrong pixels', score(png))
