"""Reconstruct gene call tree from an RNA stream (gene-ID markers on entry, CFPICFP on return)."""
import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
KNOWN = set('PIPIIIC PIPIIIP PIPIICC PIPIICF PIPIICP PIPIIFC PIPIIFF PIPIIPC PIPIIPF PIPIIPP PIIPICP '
            'PIIIIIP PCCCCCP PFFFFFP PCCIFFP PFFICCP PIIPIIP PCCPFFP PFFPCCP PFFICCF'.split())
RET = 'CFPICFP'
IDS = json.load(open(os.path.join(ROOT, 'analysis', 'gene_ids.json')))


class Node:
    def __init__(self, gid, start, parent):
        self.gid, self.start, self.end, self.parent = gid, start, None, parent
        self.kids = []
        self.draw = 0  # known commands inside (inclusive of kids)

    @property
    def addr(self):
        return IDS[self.gid][0] if self.gid in IDS else None

    def name(self):
        a = self.addr
        return '%s@%s' % (self.gid, a if a is not None else '?')


def load(path):
    r = open(path).read()
    return [r[i:i + 7] for i in range(0, len(r) - 6, 7)]


def tree(cmds):
    root = Node('ROOT', 0, None)
    cur = root
    for i, c in enumerate(cmds):
        if c in KNOWN:
            n = cur
            while n:
                n.draw += 1
                n = n.parent
        elif c == RET:
            if cur.parent:
                cur.end = i + 1
                cur = cur.parent
        else:
            n = Node(c, i, cur)
            cur.kids.append(n)
            cur = n
    while cur:
        cur.end = len(cmds)
        cur = cur.parent
    root.end = len(cmds)
    return root


def show(n, depth, maxdepth, mindraw, out=sys.stdout):
    if depth > 0:
        print('%s%s [%d..%d) draw=%d kids=%d' % ('  ' * (depth - 1), n.name(), n.start, n.end, n.draw, len(n.kids)), file=out)
    if depth < maxdepth:
        # collapse runs of identical gene with no draw
        for k in n.kids:
            if k.draw >= mindraw:
                show(k, depth + 1, maxdepth, mindraw, out)


def remove(cmds, ranges):
    """return cmds with draw commands in given [a,b) ranges removed"""
    drop = bytearray(len(cmds))
    for a, b in ranges:
        for i in range(a, b):
            drop[i] = 1
    return [c for i, c in enumerate(cmds) if not drop[i]]


if __name__ == '__main__':
    path = sys.argv[1]
    maxdepth = int(sys.argv[2]) if len(sys.argv) > 2 else 3
    mindraw = int(sys.argv[3]) if len(sys.argv) > 3 else 1
    t = tree(load(path))
    show(t, 0, maxdepth, mindraw)
