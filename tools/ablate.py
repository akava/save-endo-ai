"""Ablation montage: for each child of a gene (by address) in an RNA stream, drop its drawing commands and render.

usage: ablate.py RNA GENE_ADDR OUT.png [mindraw]
"""
import os
import subprocess
import sys

from PIL import Image, ImageDraw

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rnatree import KNOWN, load, tree, ROOT  # noqa: E402

LAYER = {'PCCPFFP', 'PFFPCCP', 'PFFICCF'}
BUILD = os.path.join(ROOT, 'build', 'build')


def find(n, addr, occurrence=0):
    hits = []

    def rec(x):
        if x.addr == addr:
            hits.append(x)
        for k in x.kids:
            rec(k)
    rec(n)
    return hits[occurrence] if len(hits) > occurrence else None


def render(cmds, path):
    tmp = path + '.rna'
    open(tmp, 'w').write(''.join(cmds))
    subprocess.run([BUILD, tmp, path], capture_output=True)
    os.remove(tmp)


def drop(cmds, a, b):
    return [c for i, c in enumerate(cmds) if not (a <= i < b and c in KNOWN and c not in LAYER)]


def isolate(cmds, a, b):
    return [c for i, c in enumerate(cmds) if (a <= i < b) or c in LAYER or c not in KNOWN]


def montage(cmds, items, out, W=180, cols=6):
    rows = (len(items) + cols - 1) // cols
    m = Image.new('RGB', (W * cols, (W + 14) * rows), 'white')
    dr = ImageDraw.Draw(m)
    tmp = os.path.join(ROOT, 'out', '_abl.png')
    for k, (label, cs) in enumerate(items):
        render(cs, tmp)
        im = Image.open(tmp).resize((W, W))
        x, y = (k % cols) * W, (k // cols) * (W + 14)
        m.paste(im, (x, y + 14))
        dr.text((x + 2, y + 1), label, fill='black')
    m.save(out)


if __name__ == '__main__':
    rna, addr, out = sys.argv[1], int(sys.argv[2]), sys.argv[3]
    mindraw = int(sys.argv[4]) if len(sys.argv) > 4 else 1
    mode = sys.argv[5] if len(sys.argv) > 5 else 'drop'
    op = isolate if mode == 'iso' else drop
    cmds = load(rna)
    t = tree(cmds)
    n = find(t, addr)
    kids = [k for k in n.kids if k.draw >= mindraw]
    items = [('all', cmds)] + [('%d %s d=%d' % (i, k.addr, k.draw), op(cmds, k.start, k.end)) for i, k in enumerate(kids)]
    montage(cmds, items, out)
    for i, k in enumerate(kids):
        print(i, k.addr, k.gid, k.start, k.end, k.draw)
