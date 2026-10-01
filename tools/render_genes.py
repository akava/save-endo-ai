"""Call given root genes from exit (no-night scene) and build a labeled montage: render_genes.py OUT.png W gene..."""
import json
import os
import sys
from concurrent.futures import ThreadPoolExecutor

from PIL import Image, ImageDraw

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from endo import OUT, ROOT, call_from_exit, combine, no_night, run  # noqa: E402

ends = {s: e for s, e in json.load(open(os.path.join(ROOT, 'analysis', 'genes.json'))) if s is not None}
out, W = sys.argv[1], int(sys.argv[2])
gs = [int(x) for x in sys.argv[3:]]


def one(s):
    png = os.path.join(OUT, 'g_%d.png' % s)
    if not os.path.exists(png):
        run(combine(no_night, lambda o: call_from_exit([(s, ends[s] - s)], o)), 'g_%d' % s, timeout=300)
    return s


with ThreadPoolExecutor(4) as ex:
    list(ex.map(one, gs))
cols = min(len(gs), 3)
rows = (len(gs) + cols - 1) // cols
m = Image.new('RGB', (W * cols, (W + 14) * rows), 'white')
dr = ImageDraw.Draw(m)
for k, s in enumerate(gs):
    try:
        im = Image.open(os.path.join(OUT, 'g_%d.png' % s)).resize((W, W))
    except Exception:
        continue
    x, y = (k % cols) * W, (k // cols) * (W + 14)
    m.paste(im, (x, y + 14))
    dr.text((x + 2, y + 1), 'G+%d' % s, fill='black')
m.save(out)
