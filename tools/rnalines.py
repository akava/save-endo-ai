"""Extract drawn line segments from an RNA file (with bucket composition) for a command index range."""
import sys

COL = {'PIPIIIC': 'k', 'PIPIIIP': 'r', 'PIPIICC': 'g', 'PIPIICF': 'y', 'PIPIICP': 'b', 'PIPIIFC': 'm', 'PIPIIFF': 'c',
       'PIPIIPC': 'w', 'PIPIIPF': 'T', 'PIPIIPP': 'O'}


def lines(path, a=0, b=None):
    s = open(path).read()
    n = len(s) // 7
    b = n if b is None else b
    x = y = 0; d = 0; mx = my = 0
    bucket = []; out = []
    DIRS = [(1, 0), (0, 1), (-1, 0), (0, -1)]  # E, S, W, N
    for i in range(0, b):
        c = s[7 * i:7 * i + 7]
        if c in COL:
            bucket.append(COL[c])
        elif c == 'PIIPICP':
            bucket = []
        elif c == 'PIIIIIP':
            x = (x + DIRS[d][0]) % 600; y = (y + DIRS[d][1]) % 600
        elif c == 'PCCCCCP':
            d = (d + 3) % 4
        elif c == 'PFFFFFP':
            d = (d + 1) % 4
        elif c == 'PCCIFFP':
            mx, my = x, y
        elif c == 'PFFICCP':
            if i >= a:
                out.append((i, (mx, my), (x, y), ''.join(sorted(bucket))))
    return out


if __name__ == '__main__':
    for l in lines(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]))[:50]:
        print(l)
