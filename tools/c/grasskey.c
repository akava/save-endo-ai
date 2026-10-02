// score RC4 keys: count of target tuft triples (b1,b2,t) appearing among the 70 generated tufts
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
static int T[64][3], NT;
static int score(const unsigned char *k, int L) {
    unsigned char S[256]; int i, j = 0, n, sc = 0;
    for (i = 0; i < 256; i++) S[i] = i;
    for (i = 0; i < 256; i++) { j = (j + S[i] + k[i % L]) & 255; unsigned char t = S[i]; S[i] = S[j]; S[j] = t; }
    i = j = 0; int o[210];
    for (n = 0; n < 210; n++) { i = (i + 1) & 255; j = (j + S[i]) & 255; unsigned char t = S[i]; S[i] = S[j]; S[j] = t; o[n] = S[(S[i] + S[j]) & 255]; }
    for (n = 0; n < 70; n++) { int b1 = o[3*n], b2 = o[3*n+1], ty = o[3*n+2] & 3; if (!ty) continue;
        for (int q = 0; q < NT; q++) if (T[q][0] == b1 && T[q][1] == b2 && T[q][2] == ty) { sc++; break; } }
    return sc;
}
int main(int argc, char **argv) {
    // argv[1]: key codes comma sep; stdin: triples
    unsigned char k[300]; int L = 0; char *p = strtok(argv[1], ",");
    while (p) { k[L++] = atoi(p); p = strtok(NULL, ","); }
    while (scanf("%d %d %d", &T[NT][0], &T[NT][1], &T[NT][2]) == 3) NT++;
    int base = score(k, L); printf("base %d of %d\n", base, NT);
    int best = 0;
    for (int pos = 0; pos < L; pos++) for (int v = 0; v < 256; v++) {   // substitution
        unsigned char kk[300]; memcpy(kk, k, L); kk[pos] = v; int s = score(kk, L);
        if (s > best || s >= 4) { if (s > best) best = s; printf("sub pos %d val %d score %d\n", pos, v, s); } }
    for (int pos = 0; pos < L; pos++) {   // deletion
        unsigned char kk[300]; memcpy(kk, k, pos); memcpy(kk + pos, k + pos + 1, L - pos - 1); int s = score(kk, L - 1);
        if (s >= 3) printf("del pos %d score %d\n", pos, s); }
    for (int pos = 0; pos <= L; pos++) for (int v = 0; v < 256; v++) {   // insertion
        unsigned char kk[300]; memcpy(kk, k, pos); kk[pos] = v; memcpy(kk + pos + 1, k + pos, L - pos); int s = score(kk, L + 1);
        if (s >= 3) printf("ins pos %d val %d score %d\n", pos, v, s); }
    printf("best %d\n", best);
}
