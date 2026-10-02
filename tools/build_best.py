"""Build the best prefix from named patches (reproducible). usage: build_best.py [OUT] [-patch ...] [+name=...]"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from endo import *  # noqa

dna = open(os.path.join(ROOT, 'data', 'endo.dna')).read()
M23 = (1 << 23) - 1


def lit_word_patch(instr_at, old, new):
    """rewrite a 24-bit word literal (template literal) found at/after genome address instr_at"""
    p = dna.index(lit(word(old & M23)), G + instr_at)
    L = len(lit(word(old & M23)))
    return lambda o: write_at(p, lit(word(new & M23)), L, o)


def natfix(n, L):
    s = nat(n)[:-1]
    assert len(s) <= L - 1
    return s + 'I' * (L - 1 - len(s)) + 'P'


def readnat(s, i):
    j = s.index('P', i)
    return sum(1 << k for k, c in enumerate(s[i:j]) if c == 'C'), j + 1


def retarget_call(a, A2, L2):
    """change the callee of the call instruction at genome address a (same instruction length)"""
    GL = 7509409
    raw = dna[G + a:G + a + 200]
    i = 2; _, i = readnat(raw, i); i += 5; aS = i; _, i = readnat(raw, i); aE = i
    i += 5; lS = i; _, i = readnat(raw, i); lE = i; i += 5; bS = i; _, i = readnat(raw, i); bE = i
    new = natfix(A2, aE - aS) + raw[aE:lS] + natfix(L2, lE - lS) + raw[lE:bS] + natfix(GL - A2 - L2, bE - bS)
    old = raw[aS:bE]
    x = next(k for k in range(len(old)) if old[k] != new[k]); y = max(k for k in range(len(old)) if old[k] != new[k]) + 1
    return lambda o: write_at(G + a + aS + x, new[x:y], y - x, o)


def fix_from(start, good):
    return fix_bases(G + start, good, dna[G + start:G + start + len(good)])


def tab(vals):
    return lit(''.join(word(v) for v in vals))


# ---------------- patches ----------------
def P_day():          # night-or-day = F (the official "turn to the sun" prefix does the same)
    return [lambda o: set_base(G + 1295, 'F', off=o)]


def P_hills():        # hillsEnabled + neutralize surfaceTransform's clear->PIPIPIF search literal
    s = G + 7159733
    def q(x): return ''.join({'I': 'C', 'C': 'F', 'F': 'P', 'P': 'IC'}[c] for c in x)
    c2 = q(q('IIIPIIPICP')); d2 = q(q('IIIPIIPICF'))
    return [lambda o: set_base(G + 210026, 'P', off=o), lambda o: write_at(s + 119, d2, len(c2), o)]


def P_bio():          # BMU cow branch: enableBioMorph + weather==2 -> ==0
    return [lambda o: set_base(G + 211299, 'P', off=o), lambda o: set_base(G + 897065, 'C', off=o)]


def P_caravan():      # vmuMode=31, registration code, position
    k = key128('Out_of_Band_II')[:15 * 9]
    return [lambda o: write_at(G + 210027, word(31), 24, o), lambda o: write_at(G + 210051, k, len(k), o),
            lit_word_patch(5048483, 390, 265), lit_word_patch(5048559, 230, 215)]


def P_clouds():       # cloud gene repaired from reversed duolc (minimal set), clouds on, no cloudy side effect
    t = json.load(open(os.path.join(ROOT, 'analysis', 'gene_table.json')))
    co, cl = t['cloud']; do, dl = t['duolc']
    good = dna[G + do:G + do + dl][::-1]
    pats = fix_bases(G + co, good, dna[G + co:G + co + cl])
    pats = [p for i, p in enumerate(pats) if i not in (2, 3, 6, 20, 30, 36)]
    p1 = dna.index('CFCCCCCCCCCCCCCCCCCCCCCIC', G + 7314955)
    return [lambda o: set_base(p1 + 1, 'C', off=o), lambda o: set_base(G + 6356797 + 5, 'I', off=o)] + pats


def P_cloudpos():
    return [lit_word_patch(6066484, 46, 28),
            lit_word_patch(6067326, 180, 176), lit_word_patch(6067402, 46, 54), lit_word_patch(6067725, 15, 11),
            lit_word_patch(6068330, 46, 32), lit_word_patch(6068653, 15, 20)]


def P_box():          # cargobox decompressor table: magenta->yellow, green->blue
    s = G + 2223573
    return [lambda o: write_at(s + 876 + 8, 'CF', 2, o), lambda o: set_base(s + 1004 + 9, 'P', off=o)]


def P_pears():
    return [retarget_call(a, 5800492, 14123) for a in (5044497, 5045287, 5047240)]


def P_cow():          # real cow: transparent shadow, CLIP->COMPOSE, no endocow hybrid
    return [lambda o: set_base(G + 897135 + 29, 'F', off=o), lambda o: write_at(G + 914550 + 3, 'PFFPCCP', 7, o),
            lambda o: kill_instr(G + 914970, 185, o)]


def P_nolambda():
    return [lambda o: kill_instr(G + 5050400, 185, o)]


def P_ducks():        # ducks branch on, motherDuck wing fix, chick position
    wing = G + 2921222 + len(lit('ICCIIIIIIIIP' + 'CIIIICCIIIIP' + 'CCIIICIIIIIP')) + 1
    return [lambda o: kill_instr(G + 5043225, 33, o), lambda o: set_base(wing, 'C', off=o),
            lit_word_patch(5045845, 525, 170), lit_word_patch(5045921, 260, 410)]


def P_whale():        # smiling whale: skip the eye-cross branch, keep the RNA before the check
    return [lambda o: kill_instr(G + 2516301, 72, o)]


def P_whale_pos():
    return [lit_word_patch(5067562, 391, 410), lit_word_patch(5067638, 176, 206)]


def P_balloon():
    return [lit_word_patch(5070089, 255, 155), lit_word_patch(5070165, 305, 324)]


def P_blades():
    return [lambda o: write_at(G + 823776, word(5), 24, o)]


def P_text():         # "Endo has morphed!" with the non-cloudy styling
    S = G + 5071568
    old = enc_str('Endo hat gemorpht'); new = enc_str('Endo has morphed|')
    a = next(i for i in range(len(old)) if old[i] != new[i]); b = max(i for i in range(len(old)) if old[i] != new[i]) + 1
    w1 = lit('ICCCCICCCCCICCCIICCIICIP' + 'ICCCIIIIIIICIIIIIIIIIIIP'); w2 = lit('IIICCIIIIICCIICCCIIIIICP' + 'CICICCIIICICCIIIIIIIIIIP')
    p = dna.index(w1, G + 5070997)
    bt = tab([7, 0, 0, 0, 1, 0, 0, 0, 0, 0]); nt = tab([29, 0, 28, 0, 0, 0, 0, 46, 0, 0])
    q = dna.index(bt, G + 7292625)
    return [lambda o: kill_instr(G + 5070888, 33, o), lambda o: write_at(S + a, new[a:b], b - a, o),
            lambda o: write_at(p, w2, len(w2), o), lambda o: write_at(q, nt, len(bt), o)]


def ADAPTER_ecc():    # correctErrors(cow-spot-middle) via the Adapter, runs before the patches
    return ''.join(push_arg(word(a)) for a in [890971, 893863, 2868]) + adapter_call(5995507, 59614)


import json  # noqa: E402
ORDER = ['day', 'hills', 'bio', 'caravan', 'clouds', 'cloudpos', 'box', 'pears', 'cow', 'nolambda', 'ducks', 'whale',
         'whale_pos', 'balloon', 'blades', 'text']


def build(names=ORDER, adapters=('ecc',)):
    pats = []
    for n in names:
        pats += globals()['P_' + n]()
    head = ''.join(globals()['ADAPTER_' + a]() for a in adapters)
    return head + combine(*pats)


if __name__ == '__main__':
    out = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, 'out', 'best.prefix')
    names = [n for n in ORDER if '-' + n not in sys.argv[2:]]
    pre = build(names)
    open(out, 'w').write(pre)
    print(len(pre))
