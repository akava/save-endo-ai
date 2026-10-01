"""Parse the printGeneTable data (inside gene G+2640835): name -> (offset, size) relative to the green zone start."""
import json
import re

G = 13615


def dec(seg):
    o, i = [], 0
    while i < len(seg):
        c = seg[i]
        if c == 'C': o.append('I'); i += 1
        elif c == 'F': o.append('C'); i += 1
        elif c == 'P': o.append('F'); i += 1
        elif seg[i:i + 2] == 'IC': o.append('P'); i += 2
        else: o.append('|'); i += 1
    return ''.join(o)


def w(s):
    return sum(1 << k for k, c in enumerate(s[:23]) if c == 'C')


def name(chars):
    return bytes((int(g[:8][::-1].replace('I', '0').replace('C', '1'), 2) + 64) % 256 for g in chars).decode('cp037')


if __name__ == '__main__':
    dna = open('data/endo.dna').read()
    d = dec(dna[G + 2640835:G + 2801284])
    out = {}
    # record: word(offset) word(size) <code> name-chars 255
    rec = re.compile(r'([IC]{23}P)([IC]{23}P)([^IC]*(?:[IC]{1,8}[^IC][^IC]*)*?)((?:[IC]{8}P)+?)CCCCCCCCP')
    for m in rec.finditer(d):
        nm = name(re.findall(r'[IC]{8}P', m.group(4)))
        out[nm] = (w(m.group(1)), w(m.group(2)))
    for k, v in out.items():
        print('%-40s %9d %8d' % (k, v[0], v[1]))
    json.dump(out, open('analysis/gene_table.json', 'w'), indent=0)
