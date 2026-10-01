"""Static call graph from out/disasm_all.txt: call = P '!rest(!A(!L)!B)' with A+L+B = genome length."""
import json
import re
import sys
from collections import defaultdict

GL = 7509409
genes = {}
cur = None
calls = defaultdict(list)
for line in open(sys.argv[1] if len(sys.argv) > 1 else 'out/disasm_all.txt'):
    m = re.match(r'==== GENE (\d+) (\d+)', line)
    if m:
        cur = int(m.group(1)); genes[cur] = int(m.group(2)); continue
    m = re.match(r'\s*(\d+)\s+P: (?:!\d+)?\(!(\d+)\(!(\d+)\)!(\d+)\)', line)
    if m:
        at, a, l, b = map(int, m.groups())
        if a + l + b == GL:
            calls[cur].append((at, a, l))
callers = defaultdict(set)
for g, cs in calls.items():
    for at, a, l in cs:
        callers[a].add(g)
json.dump({'calls': {g: cs for g, cs in calls.items()}, 'callers': {a: sorted(c) for a, c in callers.items()}},
          open('analysis/callgraph.json', 'w'))
for g in sorted(genes):
    print(g, genes[g] - g, 'callers', sorted(callers.get(g, [])), 'calls', sorted({a for _, a, _ in calls.get(g, [])}))
