"""Static map of conditional checks: P '!k(!N lit)' T '\\0' -> genome address checked, named via gene table."""
import bisect
import json
import re

GL = 7509409
genes = dict((s, e) for s, e in json.load(open('analysis/genes.json')) if s)
t = json.load(open('analysis/gene_table.json'))
items = sorted((o, l, k) for k, (o, l) in t.items() if l <= 2000)
starts = [o for o, l, k in items]


def name(x):
    i = bisect.bisect_right(starts, x) - 1
    if i >= 0:
        o, l, k = items[i]
        if x < o + max(l, 1):
            return '%s+%d' % (k, x - o) if x != o else k
    return '?'


cur = None
lines = []
for line in open('out/disasm_all.txt'):
    m = re.match(r'==== GENE (\d+) (\d+)', line)
    if m:
        cur = (int(m.group(1)), int(m.group(2))); lines.append(('G', cur)); continue
    m = re.match(r'\s*(\d+)\s+P: (.*)', line)
    if m:
        lines.append(('P', int(m.group(1)), m.group(2), cur))
out = []
for i, l in enumerate(lines):
    if l[0] != 'P':
        continue
    _, at, pat, (gs, ge) = l
    nxt = lines[i + 1][1] if i + 1 < len(lines) and lines[i + 1][0] == 'P' else None
    m = re.fullmatch(r'!(\d+)\(!(\d+)([ICFP]+)\)', pat)
    if m and nxt:
        k, n, lit = int(m.group(1)), int(m.group(2)), m.group(3)
        x = k + n - (ge - nxt)
        if 0 <= x < GL:
            out.append((gs, at, x, lit, name(x)))
names = {o: k for k, (o, l) in t.items()}
for gs, at, x, lit, nm in out:
    print('%-24s %8d  checks G+%-8d %-28s == %s' % (names.get(gs, gs), at, x, nm, lit))
