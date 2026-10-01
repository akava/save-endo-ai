"""Linear disassembler: decode (pattern, template) pairs from a DNA string."""
import sys

G = 13615  # genome start in endo.dna


class End(Exception):
    pass


class D:
    def __init__(self, s, pos=0):
        self.s, self.i = s, pos
        self.rna = []

    def peek(self, k=0):
        j = self.i + k
        return self.s[j] if j < len(self.s) else ''

    def nat(self):
        n, bit = 0, 1
        while True:
            c = self.peek()
            if not c:
                raise End
            self.i += 1
            if c == 'P':
                return n
            if c == 'C':
                n |= bit
            bit <<= 1

    def consts(self):
        out = ''
        while True:
            c = self.peek()
            if c == 'C': out += 'I'; self.i += 1
            elif c == 'F': out += 'C'; self.i += 1
            elif c == 'P': out += 'F'; self.i += 1
            elif c == 'I' and self.peek(1) == 'C': out += 'P'; self.i += 2
            else: return out

    def rnacmd(self):
        self.rna.append(self.s[self.i + 3:self.i + 10])
        self.i += 10

    def pattern(self):
        out, lvl = [], 0
        lit = ''
        def fl():
            nonlocal lit
            if lit: out.append(lit); lit = ''
        while True:
            c, c1, c2 = self.peek(), self.peek(1), self.peek(2)
            if c in 'CFP' and c:
                lit += {'C': 'I', 'F': 'C', 'P': 'F'}[c]; self.i += 1
            elif c == 'I' and c1 == 'C':
                lit += 'P'; self.i += 2
            elif c == 'I' and c1 == 'P':
                fl(); self.i += 2; out.append('!%d' % self.nat())
            elif c == 'I' and c1 == 'F':
                fl(); self.i += 3; out.append('?<%s>' % self.consts())
            elif c == 'I' and c1 == 'I' and c2 == 'P':
                fl(); self.i += 3; lvl += 1; out.append('(')
            elif c == 'I' and c1 == 'I' and c2 in 'CF' and c2:
                fl(); self.i += 3
                if lvl == 0:
                    return ''.join(out)
                lvl -= 1; out.append(')')
            elif c == 'I' and c1 == 'I' and c2 == 'I':
                self.rnacmd()
            else:
                raise End

    def template(self):
        out, lit = [], ''
        def fl():
            nonlocal lit
            if lit: out.append(lit); lit = ''
        while True:
            c, c1, c2 = self.peek(), self.peek(1), self.peek(2)
            if c in 'CFP' and c:
                lit += {'C': 'I', 'F': 'C', 'P': 'F'}[c]; self.i += 1
            elif c == 'I' and c1 == 'C':
                lit += 'P'; self.i += 2
            elif c == 'I' and c1 in 'FP' and c1:
                fl(); self.i += 2; l = self.nat(); n = self.nat()
                out.append('\\%d%s' % (n, '^%d' % l if l else ''))
            elif c == 'I' and c1 == 'I' and c2 in 'CF' and c2:
                fl(); self.i += 3; return ' '.join(out)
            elif c == 'I' and c1 == 'I' and c2 == 'P':
                fl(); self.i += 3; out.append('|%d|' % self.nat())
            elif c == 'I' and c1 == 'I' and c2 == 'I':
                self.rnacmd()
            else:
                raise End


def show_lit(t, maxlen=60):
    return t if len(t) <= maxlen else t[:maxlen] + '...(%d)' % len(t)


def disasm(s, start, end, out=sys.stdout):
    d = D(s, start)
    while d.i < end:
        at = d.i
        try:
            p = d.pattern()
            t = d.template()
        except End:
            print('%8d END' % at, file=out)
            return
        rna = (' rna:%d' % len(d.rna)) if d.rna else ''
        d.rna = []
        print('%8d  P: %s\n          T: %s%s' % (at - G, show_lit(p, 200), show_lit(t, 200), rna), file=out)


if __name__ == '__main__':
    s = open(sys.argv[1]).read().strip()
    a = int(sys.argv[2]) + G
    b = int(sys.argv[3]) + G
    disasm(s, a, b)
