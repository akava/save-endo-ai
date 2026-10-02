// RC4 key search with partially known keystream bytes.
// usage: rc4crack2 maxlen charset "i:val:mask,i:val:mask,..." [firstchar_index_from firstchar_index_to]
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
static int nk, kidx[32], kval[32], kmask[32], maxi;
static int cs[256], ncs, key[16];
static int test(int len) {
    unsigned char S[256]; int i, j = 0;
    for (i = 0; i < 256; i++) S[i] = i;
    for (i = 0; i < 256; i++) { j = (j + S[i] + key[i % len]) & 255; unsigned char t = S[i]; S[i] = S[j]; S[j] = t; }
    i = j = 0; int q = 0;
    for (int k = 0; k <= maxi; k++) {
        i = (i + 1) & 255; j = (j + S[i]) & 255; unsigned char t = S[i]; S[i] = S[j]; S[j] = t;
        int b = S[(S[i] + S[j]) & 255];
        while (q < nk && kidx[q] == k) { if ((b & kmask[q]) != kval[q]) return 0; q++; }
    }
    return 1;
}
static int lo, hi;
static int rec(int pos, int len) {
    if (pos == len) {
        if (test(len)) { printf("FOUND len %d:", len); for (int k = 0; k < len; k++) printf(" %d", key[k]); printf("\n"); fflush(stdout); }
        return 0;
    }
    int a = pos == 0 ? lo : 0, b = pos == 0 ? hi : ncs;
    for (int c = a; c < b; c++) { key[pos] = cs[c]; rec(pos + 1, len); }
    return 0;
}
int main(int argc, char** argv) {
    int maxlen = atoi(argv[1]);
    char* s = strdup(argv[2]);
    for (char* p = strtok(s, ","); p; p = strtok(NULL, ",")) cs[ncs++] = atoi(p);
    char* t = strdup(argv[3]);
    for (char* p = strtok(t, ","); p; p = strtok(NULL, ",")) { sscanf(p, "%d:%d:%d", &kidx[nk], &kval[nk], &kmask[nk]); if (kidx[nk] > maxi) maxi = kidx[nk]; nk++; }
    lo = argc > 4 ? atoi(argv[4]) : 0; hi = argc > 5 ? atoi(argv[5]) : ncs;
    for (int len = 1; len <= maxlen; len++) { fprintf(stderr, "len %d\n", len); rec(0, len); }
    fprintf(stderr, "done\n");
    return 0;
}
