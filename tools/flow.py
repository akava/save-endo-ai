"""Compact control-flow view of a gene: calls (named), conditional checks (resolved globals), jumps."""
import json
import re
import subprocess
import sys

GL = 7509409
t = json.load(open('analysis/gene_table.json'))
nm = {o: k for k, (o, l) in t.items()}
import bisect
items = sorted((o, l, k) for k, (o, l) in t.items() if l > 0)
st = [o for o, l, k in items]


def gname(x):
    i = bisect.bisect_right(st, x) - 1
    if i >= 0 and x < items[i][0] + items[i][1]:
        return items[i][2] + ('' if x == items[i][0] else '+%d' % (x - items[i][0]))
    return '?%d' % x


def flow(a, b):
    out = subprocess.run(['python3', 'tools/disasm.py', 'data/endo.dna', str(a), str(b)], capture_output=True, text=True).stdout.split('\n')
    ins = []
    for k, l in enumerate(out):
        m = re.match(r'\s*(\d+)\s+P: (.*)', l)
        if m:
            ins.append([int(m.group(1)), m.group(2), out[k + 1][12:] if k + 1 < len(out) else '', ''])
        elif l.strip().startswith('R:') and ins:
            ins[-1][3] = l.strip()[3:]
    for i, (at, p, tt, r) in enumerate(ins):
        nxt = ins[i + 1][0] if i + 1 < len(ins) else b
        m = re.search(r'\(!(\d+)\(!(\d+)\)!(\d+)\)', p)
        rr = ('   [' + r[:60] + ']') if r else ''
        if m and sum(map(int, m.groups())) == GL:
            print(at, 'CALL', nm.get(int(m.group(1)), m.group(1)), rr)
            continue
        m = re.fullmatch(r'!(\d+)\(!(\d+)([ICFP]+)\)', p)
        if m:
            x = int(m.group(1)) + int(m.group(2)) - (b - nxt)
            lit = m.group(3)
            v = sum(1 << k for k, c in enumerate(lit[:23]) if c == 'C') if len(lit) == 24 else lit
            print(at, 'IF', gname(x), '==', v, '-> skip next', rr)
            continue
        m = re.fullmatch(r'!(\d+)', p)
        if m and tt.strip() == '':
            print(at, 'JUMP ->', nxt + int(m.group(1)), rr)
            continue
        if rr:
            print(at, '...', rr)


if __name__ == '__main__':
    a = int(sys.argv[1]); b = a + t[sys.argv[2]][1] if len(sys.argv) > 2 and sys.argv[2] in t else int(sys.argv[2])
    flow(a, b)
