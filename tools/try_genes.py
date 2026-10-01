"""Call each gene after main (on top of the no-night scene) and render a montage."""
import json
import os
import sys

from PIL import Image, ImageDraw

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from endo import OUT, ROOT, call_from_exit, combine, no_night, run  # noqa: E402

genes = [(s, e) for s, e in json.load(open(os.path.join(ROOT, 'out', 'genes.json'))) if s is not None]
called = {s for s, e, *_ in json.load(open(os.path.join(ROOT, 'out', 'functions.json')))}
which = sys.argv[1] if len(sys.argv) > 1 else 'uncalled'
todo = [(s, e) for s, e in genes if (which == 'all' or s not in called)]
res = []
for s, e in todo:
    pre = combine(no_night, lambda o, s=s, e=e: call_from_exit([(s, e - s)], o))
    png, info = run(pre, 'g_%d' % s)
    res.append((s, e, info))
    print(s, e - s, info, flush=True)
W, cols = 150, 8
rows = (len(res) + cols - 1) // cols
m = Image.new('RGB', (W * cols, (W + 12) * rows), 'white')
dr = ImageDraw.Draw(m)
for k, (s, e, info) in enumerate(res):
    try:
        im = Image.open(os.path.join(OUT, 'g_%d.png' % s)).resize((W, W))
    except Exception:
        continue
    x, y = (k % cols) * W, (k // cols) * (W + 12)
    m.paste(im, (x, y + 12))
    dr.text((x + 2, y), '%d' % s, fill='black')
m.save(os.path.join(OUT, 'genes_%s.png' % which))
