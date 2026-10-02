// score many keys (one per line: comma-separated codes) against target tuft triples (file argv[1])
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
    FILE *f = fopen(argv[1], "r");
    while (fscanf(f, "%d %d %d", &T[NT][0], &T[NT][1], &T[NT][2]) == 3) NT++;
    char line[4096]; int ln = 0;
    while (fgets(line, sizeof line, stdin)) {
        ln++; unsigned char k[1024]; int L = 0; char *p = strtok(line, ",\n");
        while (p) { k[L++] = atoi(p); p = strtok(NULL, ",\n"); }
        if (!L) continue;
        int s = score(k, L); if (s >= 2) printf("%d score %d\n", ln, s);
    }
}
