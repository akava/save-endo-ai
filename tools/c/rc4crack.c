// Brute-force a FuunTech crypt (RC4) key from a purchase code: first 6 keystream bytes known.
// usage: rc4crack b0 b1 b2 b3 b4 b5 maxlen charset_codes(comma separated)
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
static unsigned char want[6];
static int cs[256], ncs, maxlen;
static int key[16];
static int test(int len) {
    unsigned char S[256]; int i, j = 0;
    for (i = 0; i < 256; i++) S[i] = i;
    for (i = 0; i < 256; i++) { j = (j + S[i] + key[i % len]) & 255; unsigned char t = S[i]; S[i] = S[j]; S[j] = t; }
    i = j = 0;
    for (int k = 0; k < 6; k++) {
        i = (i + 1) & 255; j = (j + S[i]) & 255; unsigned char t = S[i]; S[i] = S[j]; S[j] = t;
        if (S[(S[i] + S[j]) & 255] != want[k]) return 0;
    }
    return 1;
}
static int rec(int pos, int len) {
    if (pos == len) {
        if (test(len)) { printf("FOUND len %d:", len); for (int k = 0; k < len; k++) printf(" %d", key[k]); printf("\n"); fflush(stdout); return 1; }
        return 0;
    }
    for (int c = 0; c < ncs; c++) { key[pos] = cs[c]; if (rec(pos + 1, len)) return 1; }
    return 0;
}
int main(int argc, char** argv) {
    for (int k = 0; k < 6; k++) want[k] = atoi(argv[1 + k]);
    maxlen = atoi(argv[7]);
    char* s = strdup(argv[8]);
    for (char* p = strtok(s, ","); p; p = strtok(NULL, ",")) cs[ncs++] = atoi(p);
    for (int len = 1; len <= maxlen; len++) {
        fprintf(stderr, "len %d\n", len);
        if (rec(0, len)) return 0;
    }
    printf("not found\n");
    return 0;
}
