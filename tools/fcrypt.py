"""FuunTech crypt = RC4 over acids: key = EBCDIC-64 char codes; each keystream byte XORs 4 acids (2 bits each, LSB first)."""
V = {'I': 0, 'C': 1, 'F': 2, 'P': 3}
A = 'ICFP'


def keycodes(s):
    return [ch.encode('cp037')[0] - 64 for ch in s]


def crypt(key, data):
    S = list(range(256)); j = 0; k = keycodes(key) if isinstance(key, str) else key
    for i in range(256):
        j = (j + S[i] + k[i % len(k)]) % 256; S[i], S[j] = S[j], S[i]
    i = j = 0; out = []
    for n in range(0, len(data), 4):
        i = (i + 1) % 256; j = (j + S[i]) % 256; S[i], S[j] = S[j], S[i]
        b = S[(S[i] + S[j]) % 256]
        for t, c in enumerate(data[n:n + 4]):
            out.append(A[V[c] ^ ((b >> (2 * t)) & 3)])
    return ''.join(out)
