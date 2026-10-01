// DNA -> RNA executor (spec section 3).
// DNA is an implicit treap of pieces referring into an append-only arena.
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <random>
#include <chrono>

using u64 = uint64_t;

static std::vector<char> arena;

struct Node {
    u64 off, sum;
    uint32_t len, pri;
    int l, r;
};
static std::vector<Node> pool(1);
static std::vector<int> freeList;
static std::mt19937 rng(12345);

static inline u64 S(int t) { return t ? pool[t].sum : 0; }
static inline void upd(int t) { pool[t].sum = S(pool[t].l) + S(pool[t].r) + pool[t].len; }

static int newNode(u64 off, u64 len) {
    int t;
    if (!freeList.empty()) { t = freeList.back(); freeList.pop_back(); }
    else { t = (int)pool.size(); pool.push_back(Node()); }
    Node& n = pool[t];
    n.off = off; n.len = (uint32_t)len; n.sum = len; n.pri = rng(); n.l = n.r = 0;
    return t;
}

static void freeTree(int t) {
    std::vector<int> st;
    if (t) st.push_back(t);
    while (!st.empty()) {
        int x = st.back(); st.pop_back();
        if (pool[x].l) st.push_back(pool[x].l);
        if (pool[x].r) st.push_back(pool[x].r);
        freeList.push_back(x);
    }
}

static int merge(int a, int b) {
    if (!a) return b;
    if (!b) return a;
    if (pool[a].pri > pool[b].pri) {
        int m = merge(pool[a].r, b);
        pool[a].r = m; upd(a); return a;
    } else {
        int m = merge(a, pool[b].l);
        pool[b].l = m; upd(b); return b;
    }
}

// a gets first k bases, b the rest
static void split(int t, u64 k, int& a, int& b) {
    if (!t) { a = b = 0; return; }
    u64 ls = S(pool[t].l), len = pool[t].len;
    if (k <= ls) {
        int x, y; split(pool[t].l, k, x, y);
        pool[t].l = y; upd(t); a = x; b = t;
    } else if (k >= ls + len) {
        int x, y; split(pool[t].r, k - ls - len, x, y);
        pool[t].r = x; upd(t); a = t; b = y;
    } else {
        u64 cut = k - ls;
        int nn = newNode(pool[t].off + cut, len - cut);
        pool[t].len = (uint32_t)cut;
        int rr = pool[t].r; pool[t].r = 0; upd(t);
        a = t; b = merge(nn, rr);
    }
}

struct Piece { u64 off, len; };

// collect pieces covering [lo, hi) of tree t
static void collect(int t, u64 lo, u64 hi, std::vector<Piece>& out) {
    while (t && lo < hi) {
        u64 ls = S(pool[t].l), len = pool[t].len;
        if (lo < ls) collect(pool[t].l, lo, std::min(hi, ls), out);
        if (hi > ls && lo < ls + len) {
            u64 a = lo > ls ? lo - ls : 0, b = std::min(hi - ls, len);
            out.push_back({pool[t].off + a, b - a});
        }
        if (hi <= ls + len) return;
        u64 sh = ls + len;
        lo = lo > sh ? lo - sh : 0; hi -= sh; t = pool[t].r;
    }
}

// Forward reader over the treap
struct Reader {
    std::vector<int> stk;
    const char *p = nullptr, *e = nullptr;
    int cur = 0;
    u64 pos = 0; // position of *p
    void pushLeft(int t) { while (t) { stk.push_back(t); t = pool[t].l; } }
    void seek(int root, u64 k) {
        stk.clear(); pos = k; p = e = nullptr; cur = 0;
        int t = root;
        while (t) {
            u64 ls = S(pool[t].l), len = pool[t].len;
            if (k < ls) { stk.push_back(t); t = pool[t].l; }
            else if (k < ls + len) {
                cur = t;
                p = arena.data() + pool[t].off + (k - ls);
                e = arena.data() + pool[t].off + len;
                return;
            } else { k -= ls + len; t = pool[t].r; }
        }
    }
    // advance to next span; false at end
    bool nextSpan() {
        if (cur) { pushLeft(pool[cur].r); cur = 0; }
        while (!stk.empty()) {
            int t = stk.back(); stk.pop_back();
            cur = t;
            if (pool[t].len) {
                p = arena.data() + pool[t].off;
                e = p + pool[t].len;
                return true;
            }
            pushLeft(pool[t].r); cur = 0;
        }
        p = e = nullptr;
        return false;
    }
    inline int get() {
        if (p == e && !nextSpan()) return 0;
        pos++;
        return *p++;
    }
};

// Lookahead decoder on top of Reader
struct Dec {
    Reader rd;
    char buf[8]; int nb = 0;
    u64 consumed = 0;
    inline char peek(int k) {
        while (nb <= k) { int c = rd.get(); if (!c) return 0; buf[nb++] = (char)c; }
        return buf[k];
    }
    inline void adv(int n) {
        for (int i = 0; i < n; i++) peek(i);
        memmove(buf, buf + n, nb - n); nb -= n; consumed += n;
    }
};

// ---------------- global state ----------------
static int root = 0;
static std::string rna; // concatenated commands
static u64 rnaCount = 0;
static u64 cost = 0;
static u64 iters = 0;
static u64 maxLen = 0;
static bool finished = false;

struct Finish {};

enum PType { PB, PSKIP, PSEARCH, POPEN, PCLOSE };
struct PItem { PType t; char b; u64 n; std::string s; };
enum TType { TB, TREF, TLEN };
struct TItem { TType t; char b; u64 n, l; };

static Dec dec;
static FILE* logOut = nullptr;
static u64 origLen = 0, prefixLen = 0;
// arena offset of position k in tree t
static u64 originAt(int t, u64 k) {
    while (t) {
        u64 ls = S(pool[t].l), len = pool[t].len;
        if (k < ls) t = pool[t].l;
        else if (k < ls + len) return pool[t].off + (k - ls);
        else { k -= ls + len; t = pool[t].r; }
    }
    return (u64)-1;
}
static u64 traceFrom = 1, traceTo = 0;
static FILE* traceOut = stderr;

static std::string showPat(const std::vector<PItem>& p) {
    std::string o;
    for (auto& it : p) {
        if (it.t == PB) o += it.b;
        else if (it.t == PSKIP) o += "!" + std::to_string(it.n);
        else if (it.t == PSEARCH) o += "?<" + (it.s.size() > 40 ? it.s.substr(0, 40) + "..(" + std::to_string(it.s.size()) + ")" : it.s) + ">";
        else if (it.t == POPEN) o += "(";
        else o += ")";
    }
    return o;
}
static std::string showTpl(const std::vector<TItem>& t) {
    std::string o;
    for (auto& it : t) {
        if (it.t == TB) o += it.b;
        else if (it.t == TREF) o += "\\" + std::to_string(it.n) + (it.l ? "^" + std::to_string(it.l) : "");
        else o += "|" + std::to_string(it.n) + "|";
    }
    return o;
}

static void emitRna() {
    // dec at 'III'
    dec.adv(3);
    char c[7]; int k = 0;
    for (; k < 7; k++) { char x = dec.peek(0); if (!x) break; c[k] = x; dec.adv(1); }
    rna.append(c, k);
    rnaCount++;
}

static u64 nat() {
    u64 n = 0, bit = 1;
    for (;;) {
        char c = dec.peek(0);
        if (!c) throw Finish();
        dec.adv(1);
        if (c == 'P') return n;
        if (c == 'C') n |= bit;
        bit <<= 1;
    }
}

static std::string consts() {
    std::string s;
    for (;;) {
        char c = dec.peek(0);
        if (c == 'C') { dec.adv(1); s += 'I'; }
        else if (c == 'F') { dec.adv(1); s += 'C'; }
        else if (c == 'P') { dec.adv(1); s += 'F'; }
        else if (c == 'I' && dec.peek(1) == 'C') { dec.adv(2); s += 'P'; }
        else return s;
    }
}

static void pattern(std::vector<PItem>& p) {
    int lvl = 0;
    for (;;) {
        char c = dec.peek(0);
        if (c == 'C') { dec.adv(1); p.push_back({PB, 'I'}); }
        else if (c == 'F') { dec.adv(1); p.push_back({PB, 'C'}); }
        else if (c == 'P') { dec.adv(1); p.push_back({PB, 'F'}); }
        else if (c == 'I') {
            char c1 = dec.peek(1);
            if (c1 == 'C') { dec.adv(2); p.push_back({PB, 'P'}); }
            else if (c1 == 'P') { dec.adv(2); u64 n = nat(); p.push_back({PSKIP, 0, n}); }
            else if (c1 == 'F') { dec.adv(3); PItem it{PSEARCH}; it.s = consts(); p.push_back(it); }
            else if (c1 == 'I') {
                char c2 = dec.peek(2);
                if (c2 == 'P') { dec.adv(3); lvl++; p.push_back({POPEN}); }
                else if (c2 == 'C' || c2 == 'F') {
                    dec.adv(3);
                    if (lvl == 0) return;
                    lvl--; p.push_back({PCLOSE});
                } else if (c2 == 'I') { emitRna(); }
                else throw Finish();
            } else throw Finish();
        } else throw Finish();
    }
}

static void templ(std::vector<TItem>& t) {
    for (;;) {
        char c = dec.peek(0);
        if (c == 'C') { dec.adv(1); t.push_back({TB, 'I'}); }
        else if (c == 'F') { dec.adv(1); t.push_back({TB, 'C'}); }
        else if (c == 'P') { dec.adv(1); t.push_back({TB, 'F'}); }
        else if (c == 'I') {
            char c1 = dec.peek(1);
            if (c1 == 'C') { dec.adv(2); t.push_back({TB, 'P'}); }
            else if (c1 == 'F' || c1 == 'P') {
                dec.adv(2); u64 l = nat(); u64 n = nat();
                t.push_back({TREF, 0, n, l});
            } else if (c1 == 'I') {
                char c2 = dec.peek(2);
                if (c2 == 'C' || c2 == 'F') { dec.adv(3); return; }
                else if (c2 == 'P') { dec.adv(3); u64 n = nat(); t.push_back({TLEN, 0, n}); }
                else if (c2 == 'I') { emitRna(); }
                else throw Finish();
            } else throw Finish();
        } else throw Finish();
    }
}

// ---------- replacement builder ----------
struct Builder {
    std::vector<Piece> ps;
    static const u64 SMALL = 64;
    void addArena(const char* s, u64 n) {
        if (!n) return;
        u64 off = arena.size();
        arena.insert(arena.end(), s, s + n);
        if (!ps.empty() && ps.back().off + ps.back().len == off) ps.back().len += n;
        else ps.push_back({off, n});
    }
    void addChar(char c) { addArena(&c, 1); }
    void addPiece(Piece p) {
        if (!p.len) return;
        if (!ps.empty() && ps.back().off + ps.back().len == p.off) { ps.back().len += p.len; return; }
        if (p.len < SMALL) {
            // copy (arena may reallocate: copy via temp)
            char tmp[SMALL];
            memcpy(tmp, arena.data() + p.off, p.len);
            addArena(tmp, p.len);
        } else ps.push_back(p);
    }
    int build() {
        int t = 0;
        for (auto& p : ps) {
            // pieces may exceed uint32; split
            u64 off = p.off, len = p.len;
            while (len) {
                u64 k = std::min<u64>(len, 1u << 30);
                t = merge(t, newNode(off, k));
                off += k; len -= k;
            }
        }
        return t;
    }
};

static std::string asnat(u64 n) {
    std::string s;
    while (n) { s += (n & 1) ? 'C' : 'I'; n >>= 1; }
    s += 'P';
    return s;
}

static int protectCostMode = 0; // 0: final length only, 1: every level

static void quote(const std::string& d, std::string& out) {
    out.clear(); out.reserve(d.size() * 2);
    for (char c : d) {
        switch (c) {
            case 'I': out += 'C'; break;
            case 'C': out += 'F'; break;
            case 'F': out += 'P'; break;
            case 'P': out += "IC"; break;
        }
    }
}

static u64 totalLen() { return S(root); }

static void flattenIfNeeded() {
    u64 live = pool.size() - freeList.size();
    if (live < 2000000 && arena.size() < (u64)3 << 30) return;
    // flatten whole dna into a fresh arena
    std::vector<Piece> ps;
    collect(root, 0, S(root), ps);
    std::vector<char> na; na.reserve(S(root) + (64 << 20));
    for (auto& p : ps) na.insert(na.end(), arena.begin() + p.off, arena.begin() + p.off + p.len);
    freeTree(root);
    arena.swap(na);
    root = 0;
    u64 off = 0, len = arena.size();
    while (len) { u64 k = std::min<u64>(len, 1u << 30); root = merge(root, newNode(off, k)); off += k; len -= k; }
}

static bool step() {
    dec.rd.seek(root, 0); dec.nb = 0; dec.consumed = 0;
    std::vector<PItem> pat; std::vector<TItem> tpl;
    pattern(pat);
    templ(tpl);
    u64 p0 = dec.consumed;
    cost += p0;
    u64 total = totalLen();
    u64 dlen = total - p0;
    // matching on dna[p0..]
    u64 i = 0;
    std::vector<u64> c;
    std::vector<std::pair<u64,u64>> env;
    Reader rd; u64 rdPos = (u64)-1;
    bool ok = true;
    for (auto& it : pat) {
        if (it.t == PB) {
            cost += 1;
            if (i >= dlen) { ok = false; break; }
            if (rdPos != i) { rd.seek(root, p0 + i); rdPos = i; }
            int ch = rd.get(); rdPos++;
            if (logOut) {
                u64 o = originAt(root, p0 + i);
                if (o >= prefixLen && o < origLen) fprintf(logOut, "B %llu %llu %c %c\n", (unsigned long long)iters, (unsigned long long)(o - prefixLen), it.b, ch);
            }
            if (ch != it.b) { ok = false; break; }
            i++;
        } else if (it.t == PSKIP) {
            i += it.n;
            if (i > dlen) { ok = false; break; }
        } else if (it.t == PSEARCH) {
            const std::string& s = it.s;
            if (s.empty()) continue;
            // KMP
            size_t m = s.size();
            std::vector<int> pi(m, 0);
            for (size_t k = 1; k < m; k++) {
                int j = pi[k - 1];
                while (j > 0 && s[k] != s[j]) j = pi[j - 1];
                if (s[k] == s[j]) j++;
                pi[k] = j;
            }
            rd.seek(root, p0 + i); rdPos = i;
            u64 pos = i; int j = 0; bool found = false;
            while (pos < dlen) {
                int ch = rd.get(); pos++;
                while (j > 0 && ch != s[j]) j = pi[j - 1];
                if (ch == s[j]) j++;
                if ((size_t)j == m) { found = true; break; }
            }
            rdPos = pos;
            if (found) { cost += pos - i; i = pos; }
            else { cost += dlen - i; ok = false; break; }
        } else if (it.t == POPEN) {
            c.push_back(i);
        } else {
            env.push_back({c.back(), i}); c.pop_back();
        }
    }
    if (iters >= traceFrom && iters < traceTo) {
        fprintf(traceOut, "#%llu len=%llu p=%llu rna=%llu %s i=%llu\n  P: %s\n  T: %s\n", (unsigned long long)iters,
                (unsigned long long)dlen, (unsigned long long)p0, (unsigned long long)rnaCount, ok ? "OK" : "FAIL", (unsigned long long)i,
                showPat(pat).c_str(), showTpl(tpl).c_str());
        if (ok) { fprintf(traceOut, "  env:"); for (auto& e : env) fprintf(traceOut, " %llu", (unsigned long long)(e.second - e.first)); fprintf(traceOut, "\n"); }
    }
    if (!ok) {
        int a, b; split(root, p0, a, b); freeTree(a); root = b;
        return true;
    }
    // build replacement
    Builder bld;
    for (auto& t : tpl) {
        if (t.t == TB) bld.addChar(t.b);
        else if (t.t == TLEN) {
            u64 L = t.n < env.size() ? env[t.n].second - env[t.n].first : 0;
            std::string s = asnat(L); bld.addArena(s.data(), s.size());
        } else {
            if (t.n >= env.size()) continue;
            u64 lo = p0 + env[t.n].first, hi = p0 + env[t.n].second;
            if (logOut && hi > lo && hi - lo < 1000000) {
                std::vector<Piece> ps; collect(root, lo, hi, ps);
                for (auto& p : ps) if (p.off + p.len > prefixLen && p.off < origLen && p.len >= 16)
                    fprintf(logOut, "R %llu %llu %llu %llu\n", (unsigned long long)iters, (unsigned long long)(p.off - prefixLen), (unsigned long long)p.len, (unsigned long long)t.l);
            }
            if (t.l == 0) {
                std::vector<Piece> ps; collect(root, lo, hi, ps);
                for (auto& p : ps) bld.addPiece(p);
            } else {
                std::vector<Piece> ps; collect(root, lo, hi, ps);
                std::string d, q;
                for (auto& p : ps) d.append(arena.data() + p.off, p.len);
                for (u64 k = 0; k < t.l; k++) {
                    quote(d, q); d.swap(q);
                    if (protectCostMode == 1) cost += d.size();
                }
                if (protectCostMode == 0) cost += d.size();
                bld.addArena(d.data(), d.size());
            }
        }
    }
    int a, b; split(root, p0 + i, a, b); freeTree(a);
    root = merge(bld.build(), b);
    flattenIfNeeded();
    return true;
}

static std::string readDna(const char* fn) {
    FILE* f = fopen(fn, "rb");
    if (!f) { fprintf(stderr, "cannot open %s\n", fn); exit(1); }
    std::string s; char buf[1 << 16]; size_t n;
    while ((n = fread(buf, 1, sizeof buf, f)) > 0)
        for (size_t k = 0; k < n; k++) { char c = buf[k]; if (c == 'I' || c == 'C' || c == 'F' || c == 'P') s += c; }
    fclose(f);
    return s;
}

int main(int argc, char** argv) {
    // usage: dna [-p prefixfile | -s prefixstring] [-d dnafile] [-o rnafile] [-q]
    std::string prefix, dnaFile = "data/endo.dna", outFile = "out/out.rna";
    bool quiet = false, dumpDna = false;
    u64 maxIters = 0;
    for (int k = 1; k < argc; k++) {
        std::string a = argv[k];
        if (a == "-p") prefix += readDna(argv[++k]);
        else if (a == "-s") prefix += argv[++k];
        else if (a == "-d") dnaFile = argv[++k];
        else if (a == "-o") outFile = argv[++k];
        else if (a == "-q") quiet = true;
        else if (a == "--dump") dumpDna = true;
        else if (a == "-L") logOut = fopen(argv[++k], "w");
        else if (a == "-t") { traceFrom = strtoull(argv[++k], 0, 10); traceTo = strtoull(argv[++k], 0, 10); }
        else if (a == "-n") maxIters = strtoull(argv[++k], 0, 10);
        else if (a == "--protect-all") protectCostMode = 1;
        else { fprintf(stderr, "unknown arg %s\n", a.c_str()); return 1; }
    }
    std::string d = dnaFile == "-" ? std::string() : readDna(dnaFile.c_str());
    std::string all = prefix + d;
    prefixLen = prefix.size(); origLen = all.size();
    arena.reserve(all.size() + (256 << 20));
    arena.assign(all.begin(), all.end());
    pool.reserve(1 << 20);
    {
        u64 off = 0, len = arena.size();
        while (len) { u64 k = std::min<u64>(len, 1u << 30); root = merge(root, newNode(off, k)); off += k; len -= k; }
    }
    auto t0 = std::chrono::steady_clock::now();
    try {
        for (;;) {
            step();
            iters++;
            u64 L = totalLen(); if (L > maxLen) maxLen = L;
            if (maxIters && iters >= maxIters) break;
        }
    } catch (Finish&) {}
    if (dumpDna) {
        std::vector<Piece> ps; collect(root, 0, S(root), ps);
        for (auto& p : ps) fwrite(arena.data() + p.off, 1, p.len, stdout);
        printf("\n");
    }
    double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    if (logOut) fclose(logOut);
    FILE* f = fopen(outFile.c_str(), "wb");
    fwrite(rna.data(), 1, rna.size(), f); fclose(f);
    if (!quiet)
        fprintf(stderr, "iterations %llu rna %llu cost %llu maxlen %llu remaining %llu time %.2fs\n",
                (unsigned long long)iters, (unsigned long long)rnaCount, (unsigned long long)cost,
                (unsigned long long)maxLen, (unsigned long long)totalLen(), sec);
    return 0;
}
