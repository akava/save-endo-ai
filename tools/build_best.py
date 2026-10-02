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


# tuned numbers (positions etc.), overridable for re-tuning
PARAMS = dict(caravan=(267, 210), chick=(170, 410), whale=(410, 200), balloon=(159, 324), blades=5,
              cloud1=(20, 24, 15), cloud2=(176, 54, 11), cloud3=(340, 30, 20),
              h2=(224, 209), h3=(350, 257), h3s=(104, 60, 8),
              h1y=(238,), h1p=(-30, 3348), h2p=(-21, 1848), h1s=(488, 7, 5))

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
            lit_word_patch(5048483, 390, PARAMS['caravan'][0]), lit_word_patch(5048559, 230, PARAMS['caravan'][1])]


def P_clouds():       # cloud gene repaired from reversed duolc (minimal set), clouds on, no cloudy side effect
    t = json.load(open(os.path.join(ROOT, 'analysis', 'gene_table.json')))
    co, cl = t['cloud']; do, dl = t['duolc']
    good = dna[G + do:G + do + dl][::-1]
    pats = fix_bases(G + co, good, dna[G + co:G + co + cl])
    pats = [p for i, p in enumerate(pats) if i not in (2, 3, 6, 20, 30, 36)]
    p1 = dna.index('CFCCCCCCCCCCCCCCCCCCCCCIC', G + 7314955)
    return [lambda o: set_base(p1 + 1, 'C', off=o), lambda o: set_base(G + 6356797 + 5, 'I', off=o)] + pats


def P_cloudpos():
    out = []
    for (xa, ya, sa), (x0, y0, s0), key in (((6066408, 6066484, 6066807), (20, 46, 15), 'cloud1'),
                                           ((6067326, 6067402, 6067725), (180, 46, 15), 'cloud2'),
                                           ((6068254, 6068330, 6068653), (340, 46, 15), 'cloud3')):
        x, y, sz = PARAMS[key]
        for a, old, new in ((xa, x0, x), (ya, y0, y), (sa, s0, sz)):
            if old != new:
                out.append(lit_word_patch(a, old, new))
    return out


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
            lit_word_patch(5045845, 525, PARAMS['chick'][0]), lit_word_patch(5045921, 260, PARAMS['chick'][1])]


def P_whale():        # smiling whale: skip the eye-cross branch, keep the RNA before the check
    return [lambda o: kill_instr(G + 2516301, 72, o)]


def P_whale_pos():
    return [lit_word_patch(5067562, 391, PARAMS['whale'][0]), lit_word_patch(5067638, 176, PARAMS['whale'][1])]


def sun_tail_fixed(shift=None):
    """sun's mutated tail rebuilt from lightningBolt's identical tail (830 bases), keeping sun's own words"""
    S = G + 2126041; L = G + 4305372; n = 830
    a = dna[S:S + n]; b = dna[L:L + n]
    fix = list(b)
    for i in range(100, 185): fix[i] = a[i]
    for base in (248, 324):
        for i in range(base + 40, base + 76): fix[i] = a[i]
    if shift:  # moveTo(x, y) before the fill, shifted when the origin is not set by the caller
        for base, d in ((248, shift[0]), (324, shift[1])):
            blk = ''.join(fix[base:base + 76]); import re as _re
            m = list(_re.finditer(r'(?:[CF]{23}IC)', blk))[-1]
            raw = blk[m.start():m.end()]
            dec = ''.join({'C': 'I', 'F': 'C'}[c] for c in raw[:23])
            v = sum(1 << k for k, c in enumerate(dec) if c == 'C') + d
            newraw = lit(word(v))
            fix[base + m.start():base + m.end()] = list(newraw)
    lbw = lit(word(4305957)) + lit(word(245)); sw = lit(word(2126626)) + lit(word(245))
    blk = ''.join(fix[400:585]); k = blk.find(lbw); blk = blk[:k] + sw + blk[k + len(lbw):]
    fix[400:585] = list(blk)
    return ''.join(fix)


def P_sun():          # repair sun; paint it with colorSoftYellow (called instead of setOrigin in sky-day-bodies)
    pats = []
    p3 = dna.index('FFCCCCCCCCCCCCCCCCCCCCCIC', G + 7316761)
    pats.append(lambda o: write_at(p3, 'CC', 2, o))                       # sun block: weather==3 -> ==0
    good = sun_tail_fixed((480, 20)); S = 2126041
    pats += fix_bases(G + S, good, dna[G + S:G + S + 830])
    pats.append(lambda o: set_base(G + 2125111 + 19, 'I', off=o))         # sun's own 'clear' -> unknown RNA
    pats.append(lambda o: set_base(G + 2125111 + 29, 'I', off=o))         # sun's own 'yellow' -> unknown RNA
    def n12(v): return ''.join('C' if (v >> k) & 1 else 'I' for k in range(11)) + 'P'
    q = dna.index(lit(n12(25) + n12(100) + n12(50)), G + 2125262)
    q += len(lit(n12(25)))
    new = lit(n12(580) + n12(70)); old = lit(n12(100) + n12(50))
    if len(new) == len(old):
        pats.append(lambda o: write_at(q, new, len(old), o))
    else:
        raise Exception('polygon start length %d %d' % (len(old), len(new)))
    pats.append(lambda o: kill_instr(G + 7316909, 215, o))                # no setOrigin args frame
    pats.append(retarget_call(7317124, 6069368, 690))                     # setOrigin -> colorSoftYellow
    pats.append(lambda o: kill_instr(G + 7317494, 185, o))                # no resetOrigin
    return pats


def P_balloon():
    return [lit_word_patch(5070089, 255, PARAMS['balloon'][0]), lit_word_patch(5070165, 305, PARAMS['balloon'][1])]


def P_blades():
    return [lambda o: write_at(G + 823776, word(PARAMS['blades']), 24, o)]


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


def P_nopatch():      # no random grass patch (target's tufts differ in layout and colour; key unknown).
    # Skip the argument pushes too, otherwise the leftover blue-zone args break scenario's return into main.
    return [lambda o: kill_instr(G + 5035669, 5036709 - 5035669, o)]


TAIL_ALPHA = (4, 3, 2, 1, 7, 3)   # removed white, cyan, blue, black; added opaque, transparent


def P_tailalpha():   # target tail is translucent (alpha 178 = 7:3): skip bmu's checkIntegrity(cow-tail) jump and
    # replace colour commands of the tail bucket (w50 c20 b22 k15) by alpha ones, keeping the premultiplied colour
    rw, rc, rb, rk, O, T = TAIL_ALPHA
    base = G + 4892541 + 20
    al = ['IIIPIPIIPP'] * O + ['IIIPIPIIPF'] * T
    w1 = ''.join(al[:rw + rc]); w2 = ''.join(al[rw + rc:])
    pats = [lambda o: kill_instr(G + 899780, 33, o)]
    if w1:
        pats.append(lambda o: write_at(base + (50 - rw) * 10, w1, len(w1), o))
    if w2:
        pats.append(lambda o: write_at(base + (92 - rb) * 10, w2, len(w2), o))
    return pats


HILL_LITS = dict(h2=((7164867, 200), (7164943, 209)), h3=((7168555, 350), (7168631, 257)),
                 h3s=((7168957, 104), (7169033, 65), (7169109, 8)), h1y=((7161155, 242),),
                 h1p=((7161480, -28), (7161556, 3348)), h2p=((7165268, -21), (7165344, 1848)),
                 h1s=((7162007, 408), (7162083, 7), (7162159, 4)))


def P_hill2():        # surfaceTransform hill parameters (moveTo / functionSine literals), fitted to the target ridges
    out = []
    for k, lits in HILL_LITS.items():
        for (at, old), new in zip(lits, PARAMS[k]):
            if new != old:
                out.append(lit_word_patch(at, old, new))
    return out


def P_seed():        # randomInt seed (G+818624) = 8128, the 'fourth perfect number' of help-beautiful-numbers: weeds as in target
    return [lambda o: write_at(G + 818624, word(8128), 24, o)]


def jmp_at(pos, target, L):
    """jump instruction of exactly L bases at genome address pos to genome address target"""
    return 'IP' + natfix(target - (pos + L), L - 8) + 'IICIIC'


SPIRO_SUN = (530, 70)


def P_spiro():        # help page 180878 ('Synthesis of complex structures (2)': 'give it a try yourself'):
    # after scenario, main jumps into the page's three yellow spirographs, drawn small (arg0 4->1) at the sun's centre
    pats = [lambda o: write_at(G + 6536833, jmp_at(6536833, 6558851, 33), 33, o),
            lambda o: write_at(G + 6564817, jmp_at(6564817, 6570038, 32), 32, o)]
    at = G + 6558851
    for sp in [(4, 11, 9, 939, 0, 4), (4, 12, 8, 768, 0, 1), (4, 12, 8, 768, 128, 1)]:
        p1 = dna.index(lit(word(271)), at); p2 = dna.index(lit(word(293)), p1); p3 = dna.index(lit(word(4)), p2 + 10)
        for p, old, new in ((p1, 271, SPIRO_SUN[0]), (p2, 293, SPIRO_SUN[1]), (p3, 4, 1)):
            pats.append(lambda o, p=p, old=old, new=new: write_at(p, lit(word(new)), len(lit(word(old))), o))
        at = p3 + 10
    return pats


FISH_DX = 18


def P_fish():         # river fish (goldenFish adaptation tree): left fish coloured like the right ones (yellow:red 1:1,
    # one 'black' -> 'red' in both buckets of goldfishLeft), L/R drawings swapped, mkEmp(0,24) -> mkEmp(FISH_DX,24),
    # threeFish(R, emptyBox) -> threeFish(R, L): our threeFish draws its 2nd argument twice, then the first
    oL = GENES['mkGoldfishL_adaptation'][0]; oR = GENES['mkGoldfishR_adaptation'][0]
    wL = dna[G + oL + 24:G + oL + 48]; wR = dna[G + oR + 24:G + oR + 48]
    o = GENES['goldenFish_adaptation'][0]
    seq = ''.join(word(x) for x in [7327388, GENES['mkEmp_adaptation'][0], 48, 7325992, 0, 48, 7325992, 24])
    zero = dna.index(seq, G + o) + 24 * 4
    box = dna.index(''.join(word(x) for x in [7327388, GENES['emptyBox_adaptation'][0]]), G + o) + 24
    return [lambda off: write_at(box, word(GENES['mkGoldfishL_adaptation'][0]), 24, off),   # threeFish(R, L)
            lambda off: set_base(G + 2626409 + 89, 'P', None, off), lambda off: set_base(G + 2627598 + 89, 'P', None, off),
            lambda off: write_at(G + oL + 24, wR, 24, off), lambda off: write_at(G + oR + 24, wL, 24, off),
            lambda off: write_at(zero, word(FISH_DX), 24, off)]


def ADAPTER_ecc():    # correctErrors(cow-spot-middle) via the Adapter, runs before the patches
    return ''.join(push_arg(word(a)) for a in [890971, 893863, 2868]) + adapter_call(5995507, 59614)


def ADAPTER_tail():   # decrypt the cow tail (RC4 key '9546') in place, before everything else
    return crypt_call('9546', 4892541, 5212)


import json  # noqa: E402
GENES = json.load(open(os.path.join(ROOT, 'analysis', 'gene_table.json')))
ORDER = ['day', 'hills', 'bio', 'caravan', 'clouds', 'cloudpos', 'box', 'pears', 'cow', 'nolambda', 'ducks', 'whale',
         'whale_pos', 'balloon', 'blades', 'text', 'sun', 'nopatch', 'tailalpha', 'hill2', 'seed', 'spiro', 'fish']


def build(names=ORDER, adapters=('tail', 'ecc')):
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
