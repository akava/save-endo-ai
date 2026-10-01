// RNA -> image builder (spec section 4).
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <zlib.h>

static const int W = 600;
struct Px { uint8_t r, g, b, a; };
static inline bool operator==(Px x, Px y) { return x.r == y.r && x.g == y.g && x.b == y.b && x.a == y.a; }
static inline bool operator!=(Px x, Px y) { return !(x == y); }
using Bitmap = std::vector<Px>;

static std::vector<Bitmap> bitmaps;
static uint64_t sr, sg, sb, nrgb, sa, na;
static int px = 0, py = 0, mx = 0, my = 0, dir = 1; // 0 N 1 E 2 S 3 W
static bool cacheValid = false; static Px cache;

static void addRgb(int r, int g, int b) { sr += r; sg += g; sb += b; nrgb++; cacheValid = false; }
static void addA(int a) { sa += a; na++; cacheValid = false; }

static Px currentPixel() {
    if (cacheValid) return cache;
    uint64_t rc = nrgb ? sr / nrgb : 0, gc = nrgb ? sg / nrgb : 0, bc = nrgb ? sb / nrgb : 0;
    uint64_t ac = na ? sa / na : 255;
    cache = {(uint8_t)(rc * ac / 255), (uint8_t)(gc * ac / 255), (uint8_t)(bc * ac / 255), (uint8_t)ac};
    cacheValid = true;
    return cache;
}

static inline Px& at(int x, int y) { return bitmaps[0][y * W + x]; }

static void line(int x0, int y0, int x1, int y1) {
    long dx = x1 - x0, dy = y1 - y0;
    long d = std::max(labs(dx), labs(dy));
    long c = (dx * dy <= 0) ? 1 : 0;
    long x = x0 * d + (d - c) / 2, y = y0 * d + (d - c) / 2;
    Px p = currentPixel();
    for (long k = 0; k < d; k++) { at(x / d, y / d) = p; x += dx; y += dy; }
    at(x1, y1) = p;
}

static void tryfill() {
    Px nw = currentPixel(), old = at(px, py);
    if (nw == old) return;
    std::vector<int> st; st.push_back(py * W + px);
    Bitmap& b = bitmaps[0];
    while (!st.empty()) {
        int k = st.back(); st.pop_back();
        if (b[k] != old) continue;
        b[k] = nw;
        int x = k % W, y = k / W;
        if (x > 0) st.push_back(k - 1);
        if (x < W - 1) st.push_back(k + 1);
        if (y > 0) st.push_back(k - W);
        if (y < W - 1) st.push_back(k + W);
    }
}

static void compose() {
    if (bitmaps.size() < 2) return;
    Bitmap &b0 = bitmaps[0], &b1 = bitmaps[1];
    for (int k = 0; k < W * W; k++) {
        Px p0 = b0[k], p1 = b1[k]; int ia = 255 - p0.a;
        b1[k] = {(uint8_t)(p0.r + p1.r * ia / 255), (uint8_t)(p0.g + p1.g * ia / 255),
                 (uint8_t)(p0.b + p1.b * ia / 255), (uint8_t)(p0.a + p1.a * ia / 255)};
    }
    bitmaps.erase(bitmaps.begin());
}

static void clip() {
    if (bitmaps.size() < 2) return;
    Bitmap &b0 = bitmaps[0], &b1 = bitmaps[1];
    for (int k = 0; k < W * W; k++) {
        Px p1 = b1[k]; int a0 = b0[k].a;
        b1[k] = {(uint8_t)(p1.r * a0 / 255), (uint8_t)(p1.g * a0 / 255), (uint8_t)(p1.b * a0 / 255), (uint8_t)(p1.a * a0 / 255)};
    }
    bitmaps.erase(bitmaps.begin());
}

static void put32(std::vector<uint8_t>& v, uint32_t x) { for (int s = 24; s >= 0; s -= 8) v.push_back((x >> s) & 255); }
static void chunk(FILE* f, const char* type, const std::vector<uint8_t>& data) {
    std::vector<uint8_t> v; put32(v, data.size());
    v.insert(v.end(), type, type + 4); v.insert(v.end(), data.begin(), data.end());
    uint32_t crc = crc32(0, v.data() + 4, v.size() - 4); put32(v, crc);
    fwrite(v.data(), 1, v.size(), f);
}

// mode 0: rgb with alpha forced opaque; mode 1: rgba as is
static void writePng(const char* fn, const Bitmap& b, bool rgba) {
    FILE* f = fopen(fn, "wb");
    if (!f) { fprintf(stderr, "cannot write %s\n", fn); exit(1); }
    static const uint8_t sig[8] = {137, 80, 78, 71, 13, 10, 26, 10};
    fwrite(sig, 1, 8, f);
    std::vector<uint8_t> ih; put32(ih, W); put32(ih, W);
    ih.push_back(8); ih.push_back(rgba ? 6 : 2); ih.push_back(0); ih.push_back(0); ih.push_back(0);
    chunk(f, "IHDR", ih);
    int ch = rgba ? 4 : 3;
    std::vector<uint8_t> raw; raw.reserve(W * (W * ch + 1));
    for (int y = 0; y < W; y++) {
        raw.push_back(0);
        for (int x = 0; x < W; x++) {
            Px p = b[y * W + x];
            raw.push_back(p.r); raw.push_back(p.g); raw.push_back(p.b);
            if (rgba) raw.push_back(p.a);
        }
    }
    uLongf cl = compressBound(raw.size());
    std::vector<uint8_t> comp(cl);
    compress2(comp.data(), &cl, raw.data(), raw.size(), 6);
    comp.resize(cl);
    chunk(f, "IDAT", comp);
    chunk(f, "IEND", {});
    fclose(f);
}

int main(int argc, char** argv) {
    // usage: build in.rna out.png [--layers dir] (dumps all bitmaps at the end) [--steps N dir] (snapshot every N cmds)
    if (argc < 3) { fprintf(stderr, "usage: build in.rna out.png [--layers prefix] [--stop N]\n"); return 1; }
    const char* layers = nullptr; long stopAt = -1;
    for (int k = 3; k < argc; k++) {
        std::string a = argv[k];
        if (a == "--layers") layers = argv[++k];
        else if (a == "--stop") stopAt = atol(argv[++k]);
    }
    FILE* f = fopen(argv[1], "rb");
    if (!f) { fprintf(stderr, "cannot open %s\n", argv[1]); return 1; }
    std::string rna; char buf[1 << 16]; size_t n;
    while ((n = fread(buf, 1, sizeof buf, f)) > 0) rna.append(buf, n);
    fclose(f);
    bitmaps.push_back(Bitmap(W * W, Px{0, 0, 0, 0}));
    long cmds = 0, unknown = 0;
    for (size_t i = 0; i + 7 <= rna.size(); i += 7, cmds++) {
        if (stopAt >= 0 && cmds >= stopAt) break;
        std::string r = rna.substr(i, 7);
        if (r == "PIPIIIC") addRgb(0, 0, 0);
        else if (r == "PIPIIIP") addRgb(255, 0, 0);
        else if (r == "PIPIICC") addRgb(0, 255, 0);
        else if (r == "PIPIICF") addRgb(255, 255, 0);
        else if (r == "PIPIICP") addRgb(0, 0, 255);
        else if (r == "PIPIIFC") addRgb(255, 0, 255);
        else if (r == "PIPIIFF") addRgb(0, 255, 255);
        else if (r == "PIPIIPC") addRgb(255, 255, 255);
        else if (r == "PIPIIPF") addA(0);
        else if (r == "PIPIIPP") addA(255);
        else if (r == "PIIPICP") { sr = sg = sb = nrgb = sa = na = 0; cacheValid = false; }
        else if (r == "PIIIIIP") {
            if (dir == 0) py = (py + W - 1) % W; else if (dir == 1) px = (px + 1) % W;
            else if (dir == 2) py = (py + 1) % W; else px = (px + W - 1) % W;
        }
        else if (r == "PCCCCCP") dir = (dir + 3) % 4;
        else if (r == "PFFFFFP") dir = (dir + 1) % 4;
        else if (r == "PCCIFFP") { mx = px; my = py; }
        else if (r == "PFFICCP") line(px, py, mx, my);
        else if (r == "PIIPIIP") tryfill();
        else if (r == "PCCPFFP") { if (bitmaps.size() < 10) bitmaps.insert(bitmaps.begin(), Bitmap(W * W, Px{0, 0, 0, 0})); }
        else if (r == "PFFPCCP") compose();
        else if (r == "PFFICCF") clip();
        else unknown++;
    }
    writePng(argv[2], bitmaps[0], false);
    if (layers)
        for (size_t k = 0; k < bitmaps.size(); k++) {
            std::string fn = std::string(layers) + "_" + std::to_string(k) + ".png";
            writePng(fn.c_str(), bitmaps[k], true);
        }
    fprintf(stderr, "commands %ld unknown %ld bitmaps %zu\n", cmds, unknown, bitmaps.size());
    return 0;
}
