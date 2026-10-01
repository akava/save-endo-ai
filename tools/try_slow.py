"""Re-run selected genes (called from exit) with a long timeout."""
import json
import os
import sys
from concurrent.futures import ThreadPoolExecutor

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from endo import ROOT, call_from_exit, combine, no_night, run  # noqa: E402

ends = {s: e for s, e in json.load(open(os.path.join(ROOT, 'analysis', 'genes.json'))) if s is not None}
starts = [int(x) for x in sys.argv[2:]]
timeout = int(sys.argv[1])


def one(s):
    e = ends[s]
    pre = combine(no_night, lambda o: call_from_exit([(s, e - s)], o))
    png, info = run(pre, 'slow_%d' % s, timeout=timeout)
    print(s, e - s, info, flush=True)


with ThreadPoolExecutor(4) as ex:
    list(ex.map(one, starts))
