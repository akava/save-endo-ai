"""Build adaptation-tree words."""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from endo import word  # noqa

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GT = json.load(open(os.path.join(ROOT, 'analysis', 'gene_table.json')))
AAT, INTBOX = 7327388, 7325992


def a(name, *args):
    """activateAdaptationTree node: [AAT, gene, (size, sub)...] as list of ints"""
    out = [AAT, GT[name + "_adaptation"][0] if name + "_adaptation" in GT else GT[name][0]]
    for s in args:
        out += [24 * len(s)] + s
    return out


def direct(name, *args):
    out = [GT[name][0]]
    for s in args:
        out += [24 * len(s)] + s
    return out


def intbox(v):
    return [INTBOX, v]


def bases(ws):
    return ''.join(word(x) for x in ws)


def bioMul_fixed():
    tree = a('caseVar1', a('bioZero'), a('apply2', a('bioAdd'), a('var2'), a('apply2', a('bioMul'), a('var1'), a('var2'))))
    filler_words = (1128 - 24 * 3 - 24 * len(tree)) // 24
    assert filler_words % 2 == 0
    fill = []
    for _ in range(filler_words // 2):
        fill += intbox(0)
    return direct('k', tree, fill)
