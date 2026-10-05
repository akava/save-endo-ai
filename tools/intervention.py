"""How much does the solution change Endo? A = bases of Endo's DNA changed per patch, B = new code injected
(adapters, raw instructions). 'Blueprint' patches carry data that only the target can give (coordinates, polygon, text...)."""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import build_best as bb  # noqa: E402

BLUEPRINT = {'caravan', 'cloudpos', 'whale_pos', 'balloon', 'bubble', 'hill2', 'blades', 'fish', 'tailalpha', 'text'}

if __name__ == '__main__':
    pre = bb.build_merged(list(bb.ORDER))
    st = bb.LAST_STATS
    ch, raw, ad = st['changed'], st['raw'], st['adapters']
    print('%-12s %8s %8s %s' % ('patch', 'A(DNA)', 'B(new)', 'kind'))
    for n in bb.ORDER:
        print('%-12s %8d %8d %s' % (n, ch.get(n, 0), raw.get(n, 0), 'blueprint' if n in BLUEPRINT else 'repair'))
    for a, l in ad.items():
        print('%-12s %8d %8d %s' % ('ADAPTER_' + a, 0, l, 'new call'))
    A = sum(ch.values()); B = sum(raw.values()) + sum(ad.values())
    Af = sum(v for k, v in ch.items() if k in BLUEPRINT)
    print('TOTAL A=%d (blueprint %d, repair %d)  B=%d  prefix=%d' % (A, Af, A - Af, B, len(pre)))
