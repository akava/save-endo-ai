"""Find instructions whose first skip refers to a given green-zone address (relative addressing inside genes)."""
import bisect
import json
import re
import subprocess
import sys

t = json.load(open('analysis/gene_table.json'))
items = sorted((o, l, k) for k, (o, l) in t.items() if l > 0)
st = [o for o, l, k in items]


def refs(lo, hi, disasm_file):
    lines = open(disasm_file).read().split('\n')
    ins = []
    for i, l in enumerate(lines):
        m = re.match(r'\s*(\d+)\s+P: (.*)', l)
        if m:
            ins.append((int(m.group(1)), m.group(2), lines[i + 1][12:] if i + 1 < len(lines) else ''))
    for k, (at, p, tp) in enumerate(ins):
        nxt = ins[k + 1][0] if k + 1 < len(ins) else at
        j = bisect.bisect_right(st, at) - 1
        if j < 0:
            continue
        o, l, name = items[j]
        if not (o <= at < o + l):
            continue
        m = re.search(r'!(\d+)', p)
        if m:
            x = int(m.group(1)) - (o + l - nxt)
            if lo <= x < hi:
                yield name, at, x, p[:70], tp[:70]


if __name__ == '__main__':
    for r in refs(int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]):
        print(*r)
