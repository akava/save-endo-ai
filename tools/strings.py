"""Extract string literals passed to drawString: template literals made of 9-base groups (8 bits + P)."""
import re
import sys


def raw_strings(dna, a, b):
    """yield (offset, [codes]) for runs of >=3 groups [IC]{8}P in the *encoded* template (literal bases)"""
    sys.path.insert(0, 'tools')
    seg = dna[a:b]
    # template literal bases: C->I, F->C, P->F, IC->P ; decode linear text
    out = []
    i = 0
    dec = []
    pos = []
    while i < len(seg):
        c = seg[i]
        if c == 'C': dec.append('I'); pos.append(i); i += 1
        elif c == 'F': dec.append('C'); pos.append(i); i += 1
        elif c == 'P': dec.append('F'); pos.append(i); i += 1
        elif seg[i:i + 2] == 'IC': dec.append('P'); pos.append(i); i += 2
        else: dec.append('|'); pos.append(i); i += 1
    d = ''.join(dec)
    for m in re.finditer(r'(?:[IC]{8}P){2,}', d):
        codes = [int(g[:8][::-1].replace('I', '0').replace('C', '1'), 2) for g in re.findall(r'[IC]{8}P', m.group())]
        yield a + pos[m.start()], codes


if __name__ == '__main__':
    dna = open('data/endo.dna').read()
    G = 13615
    a, b = int(sys.argv[1]) + G, int(sys.argv[2]) + G
    for off, codes in raw_strings(dna, a, b):
        print(off - G, len(codes), codes)


def text(codes):
    return bytes((c + 64) % 256 for c in codes if c != 255).decode('cp037')
