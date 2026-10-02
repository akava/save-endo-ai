"""Build the best prefix from named patches (reproducible). usage: build_best.py [OUT] [-patch ...] [+name=...]"""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from endo import *  # noqa

dna = open(os.path.join(ROOT, 'data', 'endo.dna')).read()
M23 = (1 << 23) - 1


def lit_word_patch(instr_at, old, new):
    """rewrite a 24-bit word literal (template literal) found at/after genome address instr_at"""
    a, b = lit(word(old & M23)), lit(word(new & M23))
    p = dna.index(a, G + instr_at)
    if len(a) != len(b):
        return lambda o: write_at(p, b, len(a), o)
    diff = [k for k in range(len(a)) if a[k] != b[k]]
    if not diff:
        return lambda o: ''
    i, j = diff[0], diff[-1] + 1                       # rewrite only the changed span
    return lambda o: write_at(p + i, b[i:j], j - i, o)


def wdiff(pos, new, off):
    """write `new` at absolute position pos (original DNA there), rewriting only the changed span"""
    old = dna[pos:pos + len(new)]
    diff = [k for k in range(len(new)) if old[k] != new[k]]
    if not diff:
        return ''
    i, j = diff[0], diff[-1] + 1
    return write_at(pos + i, new[i:j], j - i, off)


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
PARAMS = dict(caravan=(267, 210), chick=(171, 410), whale=(410, 200), balloon=(198, 324), blades=5,
              cloud1=(20, 25, 15), cloud2=(180, 55, 10), cloud3=(340, 30, 20),
              h2=(200, 235), h3=(350, 257), h3s=(104, 60, 8),
              h1y=(218,), h1p=(-21, 1848), h2p=(-28, 3348), h1s=(408, 7, 4), h2s=(328, 13, 3), h1x=(0,))
# hills: the parabolas of hills 1 and 2 are swapped back (shoutOut: "we sabotaged the Fuun DNA by swapping some
# parabolas"); with tools/hillmodel.py the target ridges then fit exactly with the original sines and x positions

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
    return [lambda o: wdiff(G + 210027, word(31), o), lambda o: wdiff(G + 210051, k, o),
            lit_word_patch(5048483, 390, PARAMS['caravan'][0]), lit_word_patch(5048559, 230, PARAMS['caravan'][1])]


def cheapest_spans(good, cur, need):
    """group the positions `need` (where cur must become good) into writes of minimal total prefix length (DP)"""
    def cost(a, b):
        return len(set_base(a, good[a], off=5000)) if a == b else len(write_at(a, good[a:b + 1], b - a + 1, 5000))
    n = len(need); best = [0] + [10 ** 9] * n; prev = [0] * (n + 1)
    for j in range(1, n + 1):
        for i in range(j):
            c = best[i] + cost(need[i], need[j - 1])
            if c < best[j]:
                best[j], prev[j] = c, i
    segs = []; j = n
    while j > 0:
        i = prev[j]; segs.append((need[i], need[j - 1])); j = i
    return segs[::-1]


def _cloud_need():
    t = json.load(open(os.path.join(ROOT, 'analysis', 'gene_table.json')))
    co, cl = t['cloud']; do, dl = t['duolc']
    good = dna[G + do:G + do + dl][::-1]; cur = dna[G + co:G + co + cl]
    groups = []
    for i in (i for i in range(cl) if good[i] != cur[i]):
        if groups and i - groups[-1][1] <= 12:
            groups[-1][1] = i
        else:
            groups.append([i, i])
    keep = [g for k, g in enumerate(groups) if k not in (2, 3, 6, 20, 30, 36)]   # minimal set found by ablation
    return [i for a, b in keep for i in range(a, b + 1) if good[i] != cur[i]]


def P_clouds():       # cloud gene repaired from reversed duolc (minimal set), clouds on, no cloudy side effect
    t = json.load(open(os.path.join(ROOT, 'analysis', 'gene_table.json')))
    co, cl = t['cloud']; do, dl = t['duolc']
    good = dna[G + do:G + do + dl][::-1]
    cur = dna[G + co:G + co + cl]
    pats = [lambda o, a=a, b=b: (set_base(G + co + a, good[a], off=o) if a == b else
                                 write_at(G + co + a, good[a:b + 1], b - a + 1, o)) for a, b in cheapest_spans(good, cur, CLOUD_FIX)]
    p1 = dna.index('CFCCCCCCCCCCCCCCCCCCCCCIC', G + 7314955)
    return [lambda o: set_base(p1 + 1, 'C', off=o), lambda o: set_base(G + 6356797 + 5, 'I', off=o)] + pats


CLOUD_FIX = _cloud_need()


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


def P_ducks():        # ducks branch on, motherDuck wing fix; the chick is drawn over the windmill in the target:
    # its call is killed (the argument pushes stay on the blue-zone stack, all calls in between are balanced) and the
    # unused lambda-id call slot after the windmill calls `chick` instead, with the chick's origin
    wing = G + 2921222 + len(lit('ICCIIIIIIIIP' + 'CIIIICCIIIIP' + 'CCIIICIIIIIP')) + 1
    return [lambda o: kill_instr(G + 5043225, 33, o), lambda o: set_base(wing, 'C', off=o),
            lambda o: kill_instr(G + 5046450, 185, o),
            lit_word_patch(5050063, 45, PARAMS['chick'][0]), lit_word_patch(5050139, 275, PARAMS['chick'][1]),
            retarget_call(5050400, 5512462, 12361)]


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
        pats.append(lambda o: wdiff(q, new, o))
    else:
        raise Exception('polygon start length %d %d' % (len(old), len(new)))
    pats.append(lambda o: kill_instr(G + 7316909, 215, o))                # no setOrigin args frame
    pats.append(retarget_call(7317124, 6069368, 690))                     # setOrigin -> colorSoftYellow
    pats.append(lambda o: kill_instr(G + 7317494, 185, o))                # no resetOrigin
    return pats


def P_balloon():
    return [lit_word_patch(5070089, 255, PARAMS['balloon'][0]), lit_word_patch(5070165, 305, PARAMS['balloon'][1])]


def P_blades():
    return [lambda o: wdiff(G + 823776, word(PARAMS['blades']), o)]


BALLOON_TABLE = [35, 0, 0, 0, 5, 0, 0, 0, 0, 0]   # 35 black + 5 blue = (0, 0, 31): the navy mu of the target
TEXT_TABLE = [35, 0, 34, 0, 0, 0, 0, 55, 0, 0]     # (113, 183, 113): the bottom text; same black count as the balloon


def P_text():         # "Endo has morphed!": German branch with the string rewritten, coloured by useColorTable.
    # The German branch writes rotateColorVar and charColorCallback (germanColors). Both writes are dropped:
    # `balloon` has just set charColorCallback = useColorTable. Their 211 bases (with the jump to the second
    # branch) now hold one write into colorTable. It changes only entries 2..7 of the balloon's table, so the
    # two tables are chosen with the same black count; the balloon keeps its own colours for the mu.
    S = G + 5071568
    old = enc_str('Endo hat gemorpht'); new = enc_str('Endo has morphed|')
    a = next(i for i in range(len(old)) if old[i] != new[i]); b = max(i for i in range(len(old)) if old[i] != new[i]) + 1
    ob = ''.join(word(v) for v in BALLOON_TABLE); nb = ''.join(word(v) for v in TEXT_TABLE)
    diff = [k for k in range(240) if ob[k] != nb[k]]; i, j = diff[0], diff[-1] + 1
    R0, R1 = 5070888, 5071099
    N = 826012 - 240 + i - (R1 - 5070997)            # 5070921's write (ends 5070997) reaches rotateColorVar = colorTable+240
    tw = write_at(N, nb[i:j], j - i)
    X = R1 - len(tw)
    reg = (jmp_at(R0, X, 33) + dna[G + R0 + 33:G + X] if X - R0 >= 33 else jmp_at(R0, X, X - R0)) + tw
    assert len(reg) == R1 - R0
    q = dna.index(tab([7, 0, 0, 0, 1, 0, 0, 0, 0, 0]), G + 7292625)
    return [lambda o: write_at(S + a, new[a:b], b - a, o), lambda o: wdiff(G + R0, reg, o),
            lambda o: wdiff(q, tab(BALLOON_TABLE), o)]


def P_bubble():       # balloon outline of the target: polygon traced from the target (analysis/balloon_target_polygon.json,
    # tools/polyfit.py), the gradient sized to its bounding box, fill point and text inside the new circle.
    v = [tuple(p) for p in json.load(open(os.path.join(ROOT, 'analysis', 'balloon_target_polygon.json')))]
    nums = [86, v[0][0], v[0][1]]                   # the count must stay 86 (literal length): pad with zero steps
    for p0, p1 in zip(v, v[1:] + v[:1]):
        nums += [p1[0] - p0[0], p1[1] - p0[1]]
    nums += [0, 0] * (86 - len(v))
    new = lit(''.join(''.join('C' if (x & 0x7ff) >> k & 1 else 'I' for k in range(11)) + 'P' for x in nums))
    p = G + 7289758                                  # the polygon literal of `balloon` (drawPolyline argument)
    assert len(new) <= 2275
    W, H = BUBBLE['wh']; FX, FY = BUBBLE['fill']; TX, TY = BUBBLE['text']
    return [lambda o: wdiff(p, new, o), lit_word_patch(7288815, 145, W), lit_word_patch(7288891, 116, H),
            lit_word_patch(7292288, 106, FX), lit_word_patch(7292364, 37, FY),
            lit_word_patch(7293325, 98, TX), lit_word_patch(7293401, 20, TY)]


BUBBLE = dict(wh=(101, 169), fill=(64, 19), text=(56, 2))


def P_mu():           # balloon character L (lambda) -> M (mu); the glyph itself is decrypted by ADAPTER_mu
    return [lambda o: wdiff(G + 7292625 + 927, 'CCFCFCCFIC', o)]


def P_nopatch():      # no random grass patch (target's tufts differ in layout and colour; key unknown).
    # Skip the argument pushes too, otherwise the leftover blue-zone args break scenario's return into main.
    return [lambda o: kill_instr(G + 5035669, 5036709 - 5035669, o)]


def P_grass():         # grass patch of the target without its key: initFastRandom jumps straight to its `ret` (no
    # key schedule), and the RC4 state frs is written by the prefix. The state is built by tools/rc4craft.py so that
    # the stream draws the target tufts (analysis/target_patch_tufts.json); count = number of byte triples used.
    st = json.load(open(os.path.join(ROOT, 'analysis', 'grass_rc4_state.json')))
    frs = ''.join(''.join('C' if (v >> k) & 1 else 'I' for k in range(8)) + 'P' for v in st['S0'])
    return [lambda o: write_at(G + 4526406, jmp_at(4526406, 4532360, 98), 98, o),
            lambda o: write_at(G + 818666, frs, 2304, o), lit_word_patch(5035767, 70, st['count'])]


def P_biomorph():      # InitialBioMorph.hs: "enableBioMorph = False", "bioMul x y = Zero -- this is broken and must be
    # fixed!!". With enableBioMorph_adaptation -> true and bioMul repaired, `biomorph` fills bioMorphPerturb (59 words),
    # which drawGrassPatch adds to the 59 characters of its key: the grass patch of the target, no key search needed.
    from adapt_build import bases, bioMul_fixed
    fix = bases(bioMul_fixed())
    return [lambda o: wdiff(G + 7442771 + 24, fix, o), lambda o: wdiff(G + 7455983 + 24, word(7340337), o)]


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
                 h1s=((7162007, 408), (7162083, 7), (7162159, 4)), h2s=((7165795, 328), (7165871, 13), (7165947, 3)), h1x=((7161079, 0),))


def P_hill2():        # surfaceTransform hill parameters (moveTo / functionSine literals), fitted to the target ridges
    out = []
    for k, lits in HILL_LITS.items():
        for (at, old), new in zip(lits, PARAMS[k]):
            if new != old:
                out.append(lit_word_patch(at, old, new))
    return out


def P_seed():        # randomInt seed (G+818624) = 8128, the 'fourth perfect number' of help-beautiful-numbers: weeds as in target
    return [lambda o: wdiff(G + 818624, word(8128), o)]


def jmp_at(pos, target, L):
    """jump instruction of exactly L bases at genome address pos to genome address target"""
    return 'IP' + natfix(target - (pos + L), L - 8) + 'IICIIC'


SPIRO_SUN = (530, 70)


def P_spiro():        # help page 180878 ('Synthesis of complex structures (2)': 'give it a try yourself'):
    # after scenario, main jumps into the page's three yellow spirographs, drawn small (arg0 4->1) at the sun's centre
    # the target's spirographs are drawn before anticompressant (its faint overlay tints them): main's anticompressant
    # call jumps to the spirographs; after them a spare setOrigin call slot calls anticompressant, then main's end
    pats = [lambda o: wdiff(G + 6536648, jmp_at(6536648, 6558851, 185), o),
            lambda o: wdiff(G + 6564817, jmp_at(6564817, 6564969, 32), o),
            retarget_call(6564969, 5603524, 11074),
            lambda o: wdiff(G + 6565154, jmp_at(6565154, 6570038, 65), o)]
    at = G + 6558851
    for sp in [(4, 11, 9, 939, 0, 4), (4, 12, 8, 768, 0, 1), (4, 12, 8, 768, 128, 1)]:
        p1 = dna.index(lit(word(271)), at); p2 = dna.index(lit(word(293)), p1); p3 = dna.index(lit(word(4)), p2 + 10)
        for p, old, new in ((p1, 271, SPIRO_SUN[0]), (p2, 293, SPIRO_SUN[1]), (p3, 4, 1)):
            pats.append(lambda o, p=p, old=old, new=new: wdiff(p, lit(word(new)), o) if len(lit(word(new))) == len(lit(word(old))) else write_at(p, lit(word(new)), len(lit(word(old))), o))
        at = p3 + 10
    return pats


FISH_DX = 18


def P_fish():         # river fish (goldenFish adaptation tree): L/R drawings swapped (this also explains the colours:
    # goldfishLeft is yellow9 red8 black, goldfishRight yellow:red 1:1), mkEmp(0,24) -> mkEmp(FISH_DX,24),
    # threeFish(R, emptyBox) -> threeFish(R, L): our threeFish draws its 2nd argument twice, then the first
    oL = GENES['mkGoldfishL_adaptation'][0]; oR = GENES['mkGoldfishR_adaptation'][0]
    wL = dna[G + oL + 24:G + oL + 48]; wR = dna[G + oR + 24:G + oR + 48]
    o = GENES['goldenFish_adaptation'][0]
    seq = ''.join(word(x) for x in [7327388, GENES['mkEmp_adaptation'][0], 48, 7325992, 0, 48, 7325992, 24])
    zero = dna.index(seq, G + o) + 24 * 4
    box = dna.index(''.join(word(x) for x in [7327388, GENES['emptyBox_adaptation'][0]]), G + o) + 24
    return [lambda off: wdiff(box, word(GENES['mkGoldfishL_adaptation'][0]), off),   # threeFish(R, L)
            lambda off: wdiff(G + oL + 24, wR, off), lambda off: wdiff(G + oR + 24, wL, off),
            lambda off: wdiff(zero, word(FISH_DX), off)]


def P_flowers():      # flowerbed: yellow and lilac flowers are swapped vs target: swap positions 1<->2 and 3<->4
    return [lit_word_patch(4566269, 0, 34), lit_word_patch(4567420, 34, 0),
            lit_word_patch(4568581, 17, 58), lit_word_patch(4568657, 24, 12),
            lit_word_patch(4569742, 58, 17), lit_word_patch(4569818, 12, 24)]


RNA_CW, RNA_CCW, RNA_MOVE = 'IIIPFFFFFP', 'IIIPCCCCCP', 'IIIPIIIIIP'


def rna_moves(dx, dy):
    """RNA turtle moves by (dx, dy) starting and ending with heading east"""
    out = ''
    if dx > 0:
        out += RNA_MOVE * dx
    elif dx < 0:
        out += RNA_CW * 2 + RNA_MOVE * -dx + RNA_CW * 2
    if dy > 0:
        out += RNA_CW + RNA_MOVE * dy + RNA_CCW
    elif dy < 0:
        out += RNA_CCW + RNA_MOVE * -dy + RNA_CW
    return out


CUP = dict(s=(26, -20), w=(-13, 10))
FOUNTAIN = (52, -64)    # moveTo of the fountain picture, relative to the cup's origin

PGT_END = 2640835 + 160449           # end of printGeneTable
FOUNT_B0, FOUNT_B1 = 2734743, 2787365  # jump over the hidden picture / its final jump onwards into printGeneTable


def fountain_block():   # the fountain: a compressed RNA picture hidden as dead code in printGeneTable (a jump skips it).
    # It is called as a gene from its skip instruction (now a no-op) to the gene end; its final jump is replaced by a
    # `ret`. It is stored upside down: swapping the dictionary entries ccw <-> cw mirrors every turn of the picture.
    ret = dna[G + 7293764:G + 7293934][47:]                               # tail of a plain 'ret' (balloon's)
    endi = 'IP' + natfix(PGT_END - (FOUNT_B1 + 2 + 25 + len(ret)), 25) + ret
    dic = dna[G + 2735354:G + 2735354 + 320]
    assert dic[246:256] == 'IIIPCCCCCP' and dic[310:320] == 'IIIPFFFFFP'
    sw = dic[:246] + dic[310:320] + dic[256:310] + dic[246:256]
    return [lambda o: write_at(G + FOUNT_B0, jmp_at(0, 33, 33), 33, o), lambda o: write_at(G + FOUNT_B1, endi, len(endi), o),
            lambda o: wdiff(G + 2735354, sw, o),
            lit_word_patch(2734839, 550, FOUNTAIN[0]), lit_word_patch(2734915, 595, FOUNTAIN[1])]


def P_cup():          # whale in a cup of water: instead of `crater`, call ufo's rain branch (water clipped into the
    # ufo dome + translucent dome) with the RNA turtle turned by 180 degrees for the dome (it becomes a cup).
    # The turn makes the real turtle position differ from the one the DNA believes: compensated by RNA moves.
    # At the end the fountain picture is called, which returns into ufo's 'compose ret'.
    UFO, UEND = 6630730, 6630730 + 11528
    sx, sy = CUP['s']; wx, wy = CUP['w']
    start = 'IP' + 'I' * (128 - 9) + 'P' + 'IICIIC'                       # neutralise the weather checks
    head = 'IIIPIIPIIP' + 'IIIPFFPCCP' + RNA_CW * 2 + rna_moves(-14 - 2 * wx, 34 - 2 * wy)  # fill compose, turn back
    endi = head
    for _ in range(6):                                                    # call length depends on its own length
        new = head + call_instr(UEND - (6637926 + len(endi)), FOUNT_B0, PGT_END - FOUNT_B0, 6642013, UEND - 6642013)
        done = len(new) == len(endi); endi = new
        if done:
            break
    a = G + 6633668                                                       # 'addbmp clear white' before the dome mask
    return [lambda o: wdiff(G + 6632760, start, o), lambda o: write_at(G + 6637926, endi, len(endi), o),
            lambda o: write_at(a + 10, RNA_CW * 2, 20, o),                # clear white -> cw cw (mask needs only alpha)
            retarget_call(5068846, 6632760, UEND - 6632760),
            lit_word_patch(6632961, 55, 55 + wx), lit_word_patch(6633037, 42, -9 + wy),
            lit_word_patch(5068509, 392, 352 + sx), lit_word_patch(5068585, 230, 284 + sy)] + fountain_block()


def ADAPTER_ecc():    # correctErrors(cow-spot-middle) via the Adapter, runs before the patches
    return ''.join(push_arg(word(a)) for a in [890971, 893863, 2868]) + adapter_call(5995507, 59614)


def ADAPTER_mu():     # decrypt charInfo_Tempus-Bold-Huge_M (RC4 key 'no1@Ax3' from the sticky note on the monitor)
    return crypt_call('no1@Ax3', 502139, 3665)


def ADAPTER_tail():   # decrypt the cow tail (RC4 key '9546') in place, before everything else
    return crypt_call('9546', 4892541, 5212)


import json  # noqa: E402
GENES = json.load(open(os.path.join(ROOT, 'analysis', 'gene_table.json')))
ORDER = ['day', 'hills', 'bio', 'caravan', 'clouds', 'cloudpos', 'box', 'pears', 'cow', 'ducks', 'whale',
         'whale_pos', 'balloon', 'blades', 'text', 'sun', 'biomorph', 'tailalpha', 'hill2', 'seed', 'spiro', 'fish', 'flowers', 'cup', 'bubble', 'mu']


def build(names=ORDER, adapters=('tail', 'ecc', 'mu')):
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
