"""Pretty-printer of adaptation trees (help-adaptive-genes / FuunDoc)."""
import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
G = 13615
t = json.load(open(os.path.join(ROOT, 'analysis', 'gene_table.json')))
nm = {o: k for k, (o, l) in t.items()}
D = open(os.path.join(ROOT, 'data', 'endo.dna')).read().strip()
AAT, INTBOX = 7327388, 7325992


def w(pos):
    return sum(1 << k for k, c in enumerate(D[G + pos:G + pos + 23]) if c == 'C')


def name(x):
    return nm.get(x, '@%d' % x).replace('_adaptation', '')


def gene(o):
    """green-zone adaptation gene at o: [size][payload][rest]; returns expr of its payload"""
    size = w(o)
    return branch(o + 24, size)


def branch(p, size, depth=0):
    if size <= 0 or depth > 60:
        return '?'
    sel = w(p)
    if size == 24:            # leaf: just a code pointer / tree pointer
        x = w(p)
        if x in nm and nm[x].endswith('_adaptation'):
            return name(x)
        if x not in nm:     # unnamed: code pointer
            return '@%d' % x
        return name(x)
    extra = w(p + 24)
    q = p + 48 if sel in (AAT, INTBOX) else p + 24
    args = []
    while q < p + size:
        s = w(q)
        if q + 24 + s > p + size or s % 24:
            args.append('<junk>')
            break
        args.append(branch(q + 24, s, depth + 1))
        q += 24 + s
    if sel == INTBOX:
        return str(extra)
    if sel == AAT:
        head = gene_expr(extra)
    elif sel in nm:
        return '%s(%s)' % (name(sel), ', '.join(args))
    if sel == AAT:
        return '%s(%s)' % (head, ', '.join(args)) if args else head
    return '%s[%s](%s)' % (name(sel), extra, ', '.join(args))


def gene_or_code(x):
    # heuristic: an adaptation gene starts with a small size word
    s = w(x)
    if 0 < s < 100000 and s % 24 == 0:
        return gene(x)
    return '@%d' % x


def gene_expr(x):
    if x in nm and nm[x].endswith('_adaptation') and t[nm[x]][1] == 72:
        return name(x)
    if x in nm:
        return name(x)
    return gene_or_code(x)


if __name__ == '__main__':
    k = sys.argv[1]
    print(gene(t[k][0] if k in t else int(k)))
