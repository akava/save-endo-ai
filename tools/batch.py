"""Run a dict of named prefixes in parallel, print score, build a labeled montage (with target)."""
import os
import sys
from concurrent.futures import ThreadPoolExecutor

from PIL import Image, ImageDraw

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from endo import OUT, ROOT, run, score  # noqa: E402


def batch(t, montage=None, timeout=120, W=300, cols=4, jobs=4):
    def one(k):
        png, info = run(t[k], 'b_' + k, timeout=timeout)
        try:
            sc = score(png) if 'TIMEOUT' not in info else None
        except Exception:
            sc = None
        r = (k, len(t[k]), sc, info[:60])
        print(*r, flush=True)
        return r
    with ThreadPoolExecutor(jobs) as ex:
        res = list(ex.map(one, t))
    if montage:
        ks = list(t) + ['target']
        rows = (len(ks) + cols - 1) // cols
        m = Image.new('RGB', (W * cols, (W + 14) * rows), 'white')
        d = ImageDraw.Draw(m)
        for i, k in enumerate(ks):
            x, y = (i % cols) * W, (i // cols) * (W + 14)
            f = os.path.join(ROOT, 'doc', 'Target-image.png') if k == 'target' else os.path.join(OUT, 'b_%s.png' % k)
            try:
                m.paste(Image.open(f).convert('RGB').resize((W, W)), (x, y + 14))
            except Exception:
                pass
            sc = dict((r[0], r[2]) for r in res).get(k)
            d.text((x + 2, y + 1), '%s %s' % (k, sc if sc is not None else ''), fill='black')
        m.save(montage)
    return res
