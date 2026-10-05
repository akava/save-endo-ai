"""Build the best prefix from named patches (reproducible). usage: build_best.py [OUT] [-patch ...] [+name=...]"""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from endo import *  # noqa
_WRITE_AT, _SET_BASE = write_at, set_base   # real writers (build_merged temporarily records the patches' calls)

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


HILLS_HINT = 'IIPIFFCPICPFFIICPIICIPPPICIIC'   # yellow letters of contest-1998..2002, read from our renders:
# (?IFPFCC)F -> \0 P, i.e. hillsEnabled (G+210026, right after the first IFPFCC) := P


def P_hills():        # hillsEnabled (the hint prefix itself) + neutralize surfaceTransform's clear->PIPIPIF search literal
    s = G + 7159733
    def q(x): return ''.join({'I': 'C', 'C': 'F', 'F': 'P', 'P': 'IC'}[c] for c in x)
    c2 = q(q('IIIPIIPICP')); d2 = q(q('IIIPIIPICF'))
    return [lambda o: HILLS_HINT, lambda o: write_at(s + 119, d2, len(c2), o)]


def P_bio():          # BMU cow branch: enableBioMorph + weather==2 -> ==0
    return [lambda o: set_base(G + 211299, 'P', off=o), lambda o: set_base(G + 897065, 'C', off=o)]


def P_caravan():      # vmuMode=31, registration code, position
    k = key128('Out_of_Band_II')[:15 * 9]
    return [lambda o: wdiff(G + 210027, word(31), o), lambda o: wdiff(G + 210051, k, o),
            lit_word_patch(5048483, 390, PARAMS['caravan'][0]), lit_word_patch(5048559, 230, PARAMS['caravan'][1])]


def cheapest_spans(good, cur, need):
    """group the positions `need` (where cur must become good) into writes of minimal total prefix length (DP)"""
    def cost(a, b):
        return len(_SET_BASE(a, good[a], off=5000)) if a == b else len(_WRITE_AT(a, good[a:b + 1], b - a + 1, 5000))
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
    return [i for a, b in keep for i in range(a, b + 1) if good[i] != cur[i] and i not in (1111, 1693, 5440)]  # +ablation


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


def sun_restored():
    """the intended repair: sun = sunflower XOR flower (I=0, C=1, F=2, P=3), the name 'sun+flower' is the hint"""
    V = {'I': 0, 'C': 1, 'F': 2, 'P': 3}
    sf = dna[G + GENES['sunflower'][0]:][:GENES['sun'][1]]; fl = dna[G + GENES['flower'][0]:][:GENES['sun'][1]]
    return ''.join('ICFP'[V[a] ^ V[b]] for a, b in zip(sf, fl))


def sun_tail_fixed(shift=None):
    """sun's tail (830 bases from G+2126041) as restored from sunflower XOR flower; the moveTo before the fill is
    optionally shifted (the sun is painted via colorSoftYellow, which does not set the origin)"""
    S = G + 2126041; n = 830
    fix = list(sun_restored()[S - G - GENES['sun'][0]:][:n])
    if shift:  # moveTo(x, y) before the fill, shifted when the origin is not set by the caller
        for base, d in ((248, shift[0]), (324, shift[1])):
            blk = ''.join(fix[base:base + 76]); import re as _re
            m = list(_re.finditer(r'(?:[CF]{23}IC)', blk))[-1]
            raw = blk[m.start():m.end()]
            dec = ''.join({'C': 'I', 'F': 'C'}[c] for c in raw[:23])
            v = sum(1 << k for k, c in enumerate(dec) if c == 'C') + d
            newraw = lit(word(v))
            fix[base + m.start():base + m.end()] = list(newraw)
    return ''.join(fix)


def P_sun():          # repair sun; paint it with colorSoftYellow (called instead of setOrigin in sky-day-bodies)
    pats = []
    p3 = dna.index('FFCCCCCCCCCCCCCCCCCCCCCIC', G + 7316761)
    pats.append(lambda o: write_at(p3, 'CC', 2, o))                       # sun block: weather==3 -> ==0
    good = sun_tail_fixed((480, 20)); S = 2126041
    cur = dna[G + S:G + S + 830]
    pats += [lambda o, a=a, b=b: (set_base(G + S + a, good[a], off=o) if a == b else
                                  write_at(G + S + a, good[a:b + 1], b - a + 1, o))
             for a, b in cheapest_spans(good, cur, [i for i in range(830) if good[i] != cur[i]])]
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
    tw = _WRITE_AT(N, nb[i:j], j - i)
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
    # order by return addresses (a call pushes 'continue at A, rest length L' of the caller): scenario returns into
    # the spirograph block instead of main's anticompressant call; the spare setOrigin slot after the third spirograph
    # calls anticompressant and returns to the spirograph page's final 'compose' and jump to main's end
    END = GENES['main'][0] + GENES['main'][1]
    pats = [lit_word_patch(6536463, 6536648, 6558851), lit_word_patch(6536463, END - 6536648, END - 6558851),
            lambda o: wdiff(G + 6564817, jmp_at(6564817, 6564969, 32), o),
            retarget_call(6564969, 5603524, 11074),
            lit_word_patch(6564969, 6565154, 6570038), lit_word_patch(6564969, END - 6565154, END - 6570038)]
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
CUP_FIX = 1                # 'moves', or +1/-1: fix curX/curY (y sign) instead of RNA moves
CUP_BELIEVED = (84, 41)    # curX/curY when the cup's tail runs (dump with build/dna --stop-rna)
FOUNTAIN = (199, -375)    # moveTo of the fountain picture, relative to the origin of that tuft (231, 575)

PGT_END = 2640835 + 160449           # end of printGeneTable
FOUNT_B0, FOUNT_B1 = 2734743, 2787365  # jump over the hidden picture / its final jump onwards into printGeneTable


def fountain_block():   # the fountain: a compressed RNA picture hidden as dead code in printGeneTable (a jump skips it).
    # It is called as a gene (scenario's last call of `grass1`, right before the whale, now leads here) from its skip
    # instruction (now a no-op) to the gene end; its final jump is replaced by a call of `grass1` (no arguments) and a
    # `ret`, so the whale is drawn over the fountain's stem, as in the target.
    # and a `ret`. It is stored upside down: swapping the dictionary entries ccw <-> cw mirrors every turn.
    ret = dna[G + 7293764:G + 7293934][47:]                               # tail of a plain 'ret' (balloon's)
    WH, WL = 5451660, 3186                                                # grass1 (the last fixed tuft before the whale)
    ci = ''
    for _ in range(6):                                                    # the call's length depends on itself
        r = FOUNT_B1 + len(ci)
        new = call_instr(PGT_END - r, WH, WL, r, PGT_END - r)
        done = len(new) == len(ci); ci = new
        if done:
            break
    r = FOUNT_B1 + len(ci)
    endi = ci + 'IP' + natfix(PGT_END - (r + 2 + 25 + len(ret)), 25) + ret
    dic = dna[G + 2735354:G + 2735354 + 320]
    assert dic[246:256] == 'IIIPCCCCCP' and dic[310:320] == 'IIIPFFFFFP'
    sw = dic[:246] + dic[310:320] + dic[256:310] + dic[246:256]
    return [lambda o: write_at(G + FOUNT_B0, jmp_at(0, 33, 33), 33, o), lambda o: write_at(G + FOUNT_B1, endi, len(endi), o),
            lambda o: wdiff(G + 2735354, sw, o), retarget_call(5067109, FOUNT_B0, PGT_END - FOUNT_B0),
            lit_word_patch(2734839, 550, FOUNTAIN[0]), lit_word_patch(2734915, 595, FOUNTAIN[1])]


def P_cup():          # whale in a cup of water: instead of `crater`, call ufo's rain branch (water clipped into the
    # ufo dome + translucent dome) with the RNA turtle turned by 180 degrees for the dome (it becomes a cup).
    # The turn makes the real turtle position differ from the one the DNA believes: compensated by RNA moves.
    UFO, UEND = 6630730, 6630730 + 11528
    sx, sy = CUP['s']; wx, wy = CUP['w']
    start = jmp_at(6632760, 6632760 + 128, 33)                            # jump over the weather checks
    if CUP_FIX == 'moves':
        head = 'IIIPIIPIIP' + 'IIIPFFPCCP' + RNA_CW * 2 + rna_moves(-14 - 2 * wx, 34 - 2 * wy)  # fill compose, turn back
    else:   # tell the DNA where the turtle really is instead of moving it: curX/curY := believed - (dx, dy)
        dx, dy = -14 - 2 * wx, 34 - 2 * wy
        head = 'IIIPIIPIIP' + 'IIIPFFPCCP' + RNA_CW * 2
        for gaddr, v in ((804464, CUP_BELIEVED[0] - dx), (804488, CUP_BELIEVED[1] - dy * CUP_FIX)):
            ins = ''
            for _ in range(6):
                nxt = 6637926 + len(head) + len(ins)
                new = P().open().skip(6642258 - nxt + gaddr).close().skip(24).end() + T().ref(0).b(word(v & M23)).end()
                done = len(new) == len(ins); ins = new
                if done:
                    break
            head += ins
    L = len(head) + 2 + 24 + 6
    endi = head + 'IP' + natfix(6642013 - (6637926 + L), 24) + 'IICIIC'  # ...and jump to ufo's 'compose ret'
    a = G + 6633668                                                       # 'addbmp clear white' before the dome mask
    return [lambda o: wdiff(G + 6632760, start, o), lambda o: write_at(G + 6637926, endi, len(endi), o),
            lambda o: write_at(a + 10, RNA_CW * 2, 20, o),                # clear white -> cw cw (mask needs only alpha)
            retarget_call(5068846, 6632760, UEND - 6632760),
            lit_word_patch(6632961, 55, 55 + wx), lit_word_patch(6633037, 42, -9 + wy),
            lit_word_patch(5068509, 392, 352 + sx), lit_word_patch(5068585, 230, 284 + sy)] + fountain_block()


def push_key(key, off):
    """push an int9[128] key argument: its characters as a literal, the rest (terminators) copied from
    giveMeAPresent (G+93: 128 x 255) instead of 1100+ literal bases"""
    k = key128(key)[:9 * len(key)]
    rest = 1152 - len(k)
    p = P().open().skip(off + G + 93 + len(k)).close().open().skip(rest).close().open().search(BLUE_MARK).close().end()
    return p + T().ref(0).ref(1).ref(2).b(k).ref(1).end()


def crypt_call2(key, offset, size, off=0):
    """crypt_call with the cheap key push (order key, offset, size as in crypt_call); off = prefix length after it"""
    tail = push_arg(word(offset)) + push_arg(word(size)) + adapter_call(5086510, 20482)
    return push_key(key, len(tail) + off) + tail


ECC_FIX = {2069: 'P', 2132: 'C', 2134: 'I', 2141: 'C', 2182: 'C', 2283: 'I', 2306: 'C', 2364: 'P', 2379: 'C', 2450: 'I', 2535: 'C', 2610: 'C', 2615: 'C', 2715: 'I', 2730: 'C', 2827: 'I', 2861: 'I'}   # cow-spot-middle offset -> base, as restored by the program's own correctErrors(cow-spot-middle,
# cow-spot-middle-ecc, 2868) (Hamming codes, help page 84); 17 point mutations


def P_ecc():          # the 17 mutated bases of cow-spot-middle put back (instead of a new correctErrors call)
    return [lambda o, k=k, c=c: set_base(G + 890971 + k, c, off=o) for k, c in sorted(ECC_FIX.items())]


def ADAPTER_ecc(off=0):    # correctErrors(cow-spot-middle) via the Adapter, runs before the patches
    return ''.join(push_arg(word(a)) for a in [890971, 893863, 2868]) + adapter_call(5995507, 59614)


def ADAPTER_mu(off=0):     # decrypt charInfo_Tempus-Bold-Huge_M (RC4 key 'no1@Ax3' from the sticky note on the monitor)
    return crypt_call2('no1@Ax3', 502139, 3665, off)


def P_tailmain():     # cow tail decrypted by main itself: its vmu-code branch runs crypt(giveMeAPresent, offset, size)
    # when giveMeAPresent is not empty. Key '9546' (from E.T.) into giveMeAPresent, the branch's offset/size literals
    # -> cow-tail, and the branch's 'CALL vmu-code' (draws the registration page) and its 'ret' are disabled, so main
    # goes on to the scene with the tail already decrypted.
    k = key128('9546')[:9 * 4]
    return [lambda o: wdiff(G + 93, k, o),
            lit_word_patch(6530931, 2176355, 4892541), lit_word_patch(6530931, 33012, 5212),
            lambda o: kill_instr(G + 6531739, 185, o), lambda o: kill_instr(G + 6531981, 160, o)]


TAIL_KEEP_FRAME = (lambda o: kill_instr(G + 6531924, 57, o))   # the branch's cleanup deletes main's whole frame


def P_mumain():       # mu glyph decrypted by main's second crypt branch (help-beautiful-numbers, entered when
    # goodVibrations != 0). Its key is built from goodVibrations, so instead the branch's key push reads a second key
    # stored in giveMeAPresent right after the tail's: '9546' 255 'no1@Ax3' 255 (crypt reads up to the terminator).
    # Offset/size literals -> charInfo_Tempus-Bold-Huge_M; drawString(key), makeDarkness, the page call and 'ret' off.
    # Both branches end main: their cleanup (!1233 at the blue zone start) deletes main's whole frame, which the
    # second branch's entry check and the scene path still need, so both cleanups are disabled too.
    def nat_fixed(v, n):
        b = ''
        while v:
            b += 'C' if v & 1 else 'I'; v >>= 1
        return b + 'I' * (n - 1 - len(b)) + 'P'
    end = GENES['main'][0] + GENES['main'][1]
    k = 93 + 45 + (end - 6535093)            # skip of the key push at +6534930 (relative to main's end)
    return [lambda o: wdiff(G + 93 + 45, key128('no1@Ax3')[:9 * 7], o),
            lambda o: wdiff(G + 6534930 + 5, nat_fixed(k, 25), o),
            lit_word_patch(6534760, 941328, 502139), lit_word_patch(6534760, 32301, 3665),
            lambda o: wdiff(G + 1281, key128(']')[:9], o),
            lambda o: kill_instr(G + 6534390, 185, o), lambda o: kill_instr(G + 6534575, 185, o),
            lambda o: kill_instr(G + 6535568, 185, o), lambda o: kill_instr(G + 6535753, 57, o),
            lambda o: kill_instr(G + 6535810, 160, o), TAIL_KEEP_FRAME]


def ADAPTER_tail(off=0):   # decrypt the cow tail in place, before everything else; key '9546' is written in the
    # least significant bit of E.T.'s portrait on help page 112 (steganography, help page 3)
    return crypt_call2('9546', 4892541, 5212, off)


import json  # noqa: E402
GENES = json.load(open(os.path.join(ROOT, 'analysis', 'gene_table.json')))
ORDER = ['day', 'hills', 'bio', 'caravan', 'clouds', 'cloudpos', 'box', 'pears', 'cow', 'ducks', 'whale',
         'whale_pos', 'balloon', 'blades', 'text', 'sun', 'biomorph', 'tailalpha', 'hill2', 'seed', 'spiro', 'fish', 'flowers', 'cup', 'bubble', 'mu', 'tailmain', 'ecc', 'mumain']


def build(names=ORDER, adapters=('tail', 'ecc', 'mu')):
    pats = []
    for n in names:
        pats += globals()['P_' + n]()
    out = combine(*pats)
    for a in reversed(adapters):
        out = globals()['ADAPTER_' + a](len(out)) + out
    return out




NOMERGE = False


def build_merged(names=ORDER, adapters=(), plain=()):
    """build() with all same-length replacements merged globally: the patches are recorded, overlaid on the original
    DNA, and the changed bases are written as spans chosen by DP over the real write cost."""
    import endo as E
    rec = []
    real_w, real_s = E.write_at, E.set_base

    def w_rec(pos, bases, oldlen, off=0):
        rec.append((pos, bases, oldlen)); return ''

    def s_rec(pos, base, old=None, off=0):
        rec.append((pos, base, 1)); return ''
    g = globals()
    others = []
    E.write_at, E.set_base, g['write_at'], g['set_base'] = w_rec, s_rec, w_rec, s_rec
    owner = []
    raw = {}
    try:
        for n in names:
            if n in plain:
                continue
            for p in g['P_' + n]():
                k = len(rec)
                out = p(0)
                owner += [n] * (len(rec) - k)
                if out:
                    others.append(p)
                    raw[n] = raw.get(n, 0) + len(out)
                if out and len(rec) > k:
                    raise Exception('mixed patch ' + n)
    finally:
        E.write_at, E.set_base, g['write_at'], g['set_base'] = real_w, real_s, real_w, real_s
    final = {}
    keep = []
    per = {}                                   # intervention stats: bases of Endo's DNA each patch changes
    for (pos, bases, oldlen), n in zip(rec, owner):
        if len(bases) != oldlen:
            per.setdefault(n, set()).update(range(pos, pos + oldlen))
        else:
            per.setdefault(n, set()).update(pos + i for i, c in enumerate(bases) if dna[pos + i] != c)
    last = {}
    for (pos, bases, oldlen), n in zip(rec, owner):
        for i in range(oldlen):
            last[pos + i] = n
    for pos, bases, oldlen in rec:
        if len(bases) != oldlen:
            keep.append((pos, bases, oldlen)); continue
        for i, c in enumerate(bases):
            final[pos + i] = c
    TA, TL = G + 4892541, 5212
    if 'tail' not in adapters:                 # tail decrypted later by main: store our plaintext edits encrypted
        from fcrypt import crypt as rc4
        tplain = list(rc4('9546', dna[TA:TA + TL]))
        hit = False
        for p in list(final):
            if TA <= p < TA + TL:
                tplain[p - TA] = final[p]; hit = True
        if hit:
            enc = rc4('9546', ''.join(tplain))
            for i, c in enumerate(enc):
                if dna[TA + i] != c:
                    final[TA + i] = c
                elif TA + i in final:
                    del final[TA + i]
    prot = ([(TA, TA + TL)] if 'tail' in adapters else []) + ([(G + 502139, G + 502139 + 3665)] if 'mu' in adapters else [])
    chg = {}
    for p_, c in final.items():
        if dna[p_] != c:
            chg[last.get(p_, '?')] = chg.get(last.get(p_, '?'), 0) + 1
    global LAST_STATS
    LAST_STATS = dict(changed=chg, raw=raw, adapters={a: len(g['ADAPTER_' + a](0)) for a in adapters})
    ch = sorted(p for p, c in final.items() if dna[p] != c or any(a <= p < b for a, b in prot))
    runs = []
    for p in ch:
        if runs and p == runs[-1][1] + 1:
            runs[-1][1] = p
        else:
            runs.append([p, p])

    def seg(a, b):
        return ''.join(final.get(i, dna[i]) for i in range(a, b + 1))
    # cost of writing [a, b] = fixed part (depends on a and the length) + quoted length of the bases
    R = len(runs)
    pc = {}                                   # quoted extra per run interval via cumulative counts over run bounds
    lo, hi = (runs[0][0], runs[-1][1]) if runs else (0, 0)
    content = seg(lo, hi) if runs else ''
    cum = [0]
    for c in content:
        cum.append(cum[-1] + (2 if c == 'P' else 1))

    def cost(a, b):
        return len(real_w(a, '', b - a + 1, 20000)) - 3 + (cum[b - lo + 1] - cum[a - lo]) + 3
    barrier = [NOMERGE or any(a <= runs[k + 1][0] - 1 and runs[k][1] + 1 < b for a, b in prot) for k in range(R - 1)]
    best = [0] + [10 ** 9] * R; prev = [0] * (R + 1)
    for j in range(1, R + 1):
        for i in range(j - 1, -1, -1):
            if runs[j - 1][1] - runs[i][0] > 3000 or (i < j - 1 and barrier[i]):
                break
            c = best[i] + cost(runs[i][0], runs[j - 1][1])
            if c < best[j]:
                best[j], prev[j] = c, i
    spans = []; j = R
    while j > 0:
        i = prev[j]; spans.append((runs[i][0], runs[j - 1][1])); j = i
    pats = [lambda o, a=a, b=b: real_w(a, seg(a, b), b - a + 1, o) for a, b in spans[::-1]]
    pats += [lambda o, pos=pos, bases=bases, oldlen=oldlen: real_w(pos, bases, oldlen, o) for pos, bases, oldlen in keep]
    pats += others
    for n in plain:
        pats += g['P_' + n]()
    out = combine(*pats)
    for a in reversed(adapters):
        out = g['ADAPTER_' + a](len(out)) + out
    return out


if __name__ == '__main__':
    out = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, 'out', 'best.prefix')
    names = [n for n in ORDER if '-' + n not in sys.argv[2:]]
    pre = build_merged(names)  # build(names) writes each patch separately
    open(out, 'w').write(pre)
    print(len(pre))
