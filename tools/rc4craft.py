"""Build an RC4 initial state (permutation, i = j = 0) whose output stream draws given grass tufts.
Each tuft is a triple (b1, b2, type): bytes b1, b2 and a byte with (b & 3) == type; tufts may come in any order.
Depth-first search with randomized candidate order; the free entries of the state are assigned lazily."""
import random
import sys

sys.setrecursionlimit(100000)


def prga(S0, n):
    S = list(S0); i = j = 0; out = []
    for _ in range(n):
        i = (i + 1) % 256; j = (j + S[i]) % 256; S[i], S[j] = S[j], S[i]; out.append(S[(S[i] + S[j]) % 256])
    return out


def craft(tufts, seed=0, width=6, budget=2000000, max_dummies=40):
    rnd = random.Random(seed)
    cur = [None] * 256; orig = list(range(256)); S0 = [None] * 256; used = [False] * 256
    order = []
    nodes = [0]
    dummies = [0]
    maxstep = [0]; maxch = [0]; best = [None]
    steps = [0]

    def assign(p, v):
        cur[p] = v; S0[orig[p]] = v; used[v] = True

    def unassign(p, v):
        cur[p] = None; S0[orig[p]] = None; used[v] = False

    ALL = set(range(256))

    def allowed_for(step, chosen):
        k = step % 3
        if k == 0:
            return ALL if dummies[0] < max_dummies else {t[0] for idx, t in enumerate(tufts) if idx not in chosen}
        if order[-1] is None:
            return ALL if k == 1 else {v for v in range(256) if v & 3 == 0}
        t = tufts[order[-1]]
        return {t[1]} if k == 1 else {v for v in range(256) if v & 3 == t[2]}

    def rec(step, i, j, chosen):
        nodes[0] += 1
        if step > maxstep[0]: maxstep[0] = step
        if step % 3 == 0 and len(chosen) > maxch[0]:
            maxch[0] = len(chosen); best[0] = (list(S0), list(order), step)
        if nodes[0] > budget:
            return False
        if len(chosen) == len(tufts) and step % 3 == 0:
            steps[0] = step
            return True
        allowed = allowed_for(step, chosen)
        pref = {t[0] for idx, t in enumerate(tufts) if idx not in chosen} if step % 3 == 0 else allowed
        i2 = (i + 1) % 256
        if cur[i2] is not None:
            cands_i = [cur[i2]]
        else:
            fv = [v for v in range(256) if not used[v]]; rnd.shuffle(fv)
            freej = [v for v in fv if cur[(j + v) % 256] is None and (j + v) % 256 != i2]
            cands_i = (freej + [v for v in fv if v not in freej])[:width]
        for vi in cands_i:
            newi = cur[i2] is None
            if newi:
                assign(i2, vi)
            j2 = (j + cur[i2]) % 256
            if cur[j2] is not None:
                cands_j = [cur[j2]]
            else:
                fv = [v for v in range(256) if not used[v]]
                rnd.shuffle(fv)
                # prefer S[j] values that make the output slot already hold an allowed value, or be free
                good = []
                for v in fv:
                    t = (cur[i2] + v) % 256
                    if t == i2: ov = v
                    elif t == j2: ov = cur[i2]
                    else: ov = cur[t]
                    if ov is not None and ov in pref: good.insert(0, v)
                    elif ov is None: good.insert(len([g for g in good]) if False else len(good), v)
                    elif ov in allowed: good.append(v)
                    if len(good) >= width: break
                cands_j = good
            for vj in cands_j:
                newj = cur[j2] is None
                if newj:
                    assign(j2, vj)
                cur[i2], cur[j2] = cur[j2], cur[i2]; orig[i2], orig[j2] = orig[j2], orig[i2]
                t = (cur[i2] + cur[j2]) % 256
                if cur[t] is not None:
                    outs = [cur[t]]
                else:
                    pv = [v for v in pref if not used[v]]; rnd.shuffle(pv)
                    av = [v for v in allowed if not used[v] and v not in pref]; rnd.shuffle(av)
                    outs = (pv + av)[:width]
                for ov in outs:
                    if ov not in allowed:
                        continue
                    newt = cur[t] is None
                    if newt:
                        assign(t, ov)
                    k = step % 3
                    pushed = False
                    isdummy = False
                    if k == 0:
                        cands = [idx for idx, tf in enumerate(tufts) if idx not in chosen and tf[0] == ov]
                        if cands:
                            order.append(cands[0]); chosen.add(cands[0])
                        else:
                            order.append(None); dummies[0] += 1; isdummy = True
                        pushed = True
                    if rec(step + 1, i2, j2, chosen):
                        return True
                    if pushed:
                        q = order.pop()
                        if q is None:
                            dummies[0] -= 1
                        else:
                            chosen.discard(q)
                    if newt:
                        unassign(t, ov)
                cur[i2], cur[j2] = cur[j2], cur[i2]; orig[i2], orig[j2] = orig[j2], orig[i2]
                if newj:
                    unassign(j2, vj)
            if newi:
                unassign(i2, vi)
        return False

    if not rec(0, 0, 0, set()):
        if best[0] is None:
            return None, None
        S0p, order, _ = best[0]
        S0[:] = S0p
        used[:] = [False] * 256
        for v in S0:
            if v is not None:
                used[v] = True
        order = list(order)
    rest = [v for v in range(256) if not used[v]]; rnd.shuffle(rest)
    for k in range(256):
        if S0[k] is None:
            S0[k] = rest.pop()
    return S0, [None if i is None else tufts[i] for i in order]
