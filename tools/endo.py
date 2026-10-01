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


def set_base(pos, base, old=None, off=0):
    """replace original base at absolute position pos with `base`"""
    p = P().open().skip(pos + off).close()
    if old:
        p.b(old)
    else:
        p.skip(1)
    return p.end() + T().ref(0).b(base).end()


def write_at(pos, bases, oldlen, off=0):
    """replace oldlen bases at pos by `bases`"""
    p = P().open().skip(pos + off).close().skip(oldlen).end()
    return p + T().ref(0).b(bases).end()


# ---------- running ----------
def run(prefix, name='run', layers=False, extra=(), timeout=60):
    os.makedirs(OUT, exist_ok=True)
    pf = os.path.join(OUT, name + '.prefix')
    with open(pf, 'w') as f:
        f.write(prefix)
    rna = os.path.join(OUT, name + '.rna')
    png = os.path.join(OUT, name + '.png')
    try:
        r = subprocess.run([DNA, '-d', os.path.join(ROOT, 'data', 'endo.dna'), '-p', pf, '-o', rna, *extra],
                           capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        return png, 'TIMEOUT'
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


G = 13615  # genome start inside endo.dna


def jump_instr(n):
    """pattern !n, empty template: drops n bases from the front (forward jump)"""
    return 'IP' + nat(n) + 'IIC' + 'IIC'


def stub_gene(start, end, off=0):
    """make gene [start,end) (genome coords) return immediately; keeps its length"""
    at = start + 10  # after III+ID
    epi = end - 75
    # instruction length depends on n; iterate to fixed point
    n = 0
    for _ in range(5):
        ins = jump_instr(n)
        n = epi - (at + len(ins))
    ins = jump_instr(n)
    assert at + len(ins) + n == epi
    return write_at(G + at, ins, len(ins), off)


BOOT_FRAME_AT = 13562  # raw offset of the two (0,0) return words in the bootstrap template


def word(n):
    """24-base stack word (23-bit nat, LSB first, terminated by P)"""
    return ''.join('C' if (n >> k) & 1 else 'I' for k in range(23)) + 'P'


def call_after_main(frames, off=0):
    """make main return into a chain of genes: frames = [(addr, len), ...] (genome coords)"""
    old = ('C' * 23 + 'IC') * 2 + 'IIC'
    words = ''.join(word(a) + word(l) for a, l in list(frames) + [(0, 0)])
    raw = lit(words) + 'IIC'  # raw template text to put into the bootstrap
    p = P().open().skip(BOOT_FRAME_AT + off).close().skip(len(old)).end()
    return p + T().ref(0).b(raw).end()


def combine(*patches):
    """patches: callables off -> prefix string; off = length of prefix that follows the instruction"""
    out = ''
    for f in reversed(patches):
        out = f(len(out)) + out
    return out


def no_night(off=0):
    """night overlay gene G+6628979: 44 opaque -> transparent (copies the following 43 transparent cmds)"""
    pos = G + 6628979 + 20
    tr = 'IIIPIPIIPF'
    return P().open().skip(pos + off).close().skip(440).open().skip(430).close().end() + T().ref(0).ref(1).ref(1).b(tr).end()


EXIT_GENE = (2248207, 2248601)  # called by main at the very end; its first instruction truncates the DNA
GENOME_LEN = 7509409


def call_instr(rest, addr, length, ret, retlen):
    """code that calls gene [addr, addr+length) and then returns to genome[ret, ret+retlen);
    `rest` = number of bases following this instruction in the current gene copy (dropped)"""
    p = P().skip(rest).open().skip(addr).open().skip(length).close().skip(GENOME_LEN - addr - length).close().end()
    return p + T().ref(0).ref(1).b(word(ret) + word(retlen)).end()


def call_from_exit(frames, off=0):
    """replace exit gene's first instruction (after its ID) by calls to genes in `frames`; returns to the exit gene's epilogue"""
    s, e = EXIT_GENE
    at = s  # overwrite the gene-ID RNA too (harmless)
    epi = 6630471  # a return instruction in another gene (all epilogues are identical), 235 bases to its end
    RETLEN = 235
    code = ''
    # chain: call f0 -> returns to next call instr ... -> last returns to epi
    # build backwards: each instr lives in genome at known address
    instrs = []
    pos = at
    # we need addresses of each instr: lengths depend on values; iterate
    lens = [0] * len(frames)
    for _ in range(4):
        pos = at
        addrs = []
        for k in range(len(frames)):
            addrs.append(pos)
            pos += lens[k]
        new = []
        for k, (a, l) in enumerate(frames):
            nxt = addrs[k + 1] if k + 1 < len(frames) else epi
            rest = e - (addrs[k] + lens[k])
            ins = call_instr(rest, a, l, nxt, RETLEN if nxt == epi else e - nxt)
            new.append(ins)
        if [len(x) for x in new] == lens:
            break
        lens = [len(x) for x in new]
    code = ''.join(new)
    assert at + len(code) <= e, (len(code), e - at)
    return write_at(G + at, code, len(code), off)
