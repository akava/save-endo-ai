"""Decode the big literal array of the instruction at genome address `at` as 12-bit signed numbers."""
import sys
sys.path.insert(0, __import__('os').path.dirname(__file__))
from disasm import D, G  # noqa


def lits(at):
    s = open('data/endo.dna').read().strip()
    d = D(s, G + at)
    d.pattern()
    t = d.template()
    return max(t.split(), key=len)


def nums(lit, bits=12):
    v = []
    for k in range(0, len(lit) - bits + 1, bits):
        w = lit[k:k + bits]
        if w[-1] != 'P':
            break
        x = sum(1 << j for j, c in enumerate(w[:bits - 1]) if c == 'C')
        v.append(x - (1 << (bits - 1)) if x >= 1 << (bits - 2) else x)
    return v


if __name__ == '__main__':
    print(nums(lits(int(sys.argv[1]))))
