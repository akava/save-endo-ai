import sys,json,itertools; sys.path.insert(0, __import__("os").path.dirname(__import__("os").path.abspath(__file__)))
from endo import *
from batch import batch
dna=open('data/endo.dna').read()
ecc=''.join(push_arg(word(a)) for a in [890971,893863,2868])+adapter_call(5995507,59614)
M=(1<<23)-1
def wpatch(instr_at, old, new):
    p=dna.index(lit(word(old&M)),G+instr_at); L=len(lit(word(old&M)))
    return lambda o: write_at(p,lit(word(new&M)),L,o)
def build(base, patches):
    rest=base[len(ecc):]
    return ecc+combine(*patches, lambda o: rest)
def tune(base, params, grids, tag, jobs=4):
    """params: list of (instr_at, oldval); grids: list of value lists"""
    t={}
    for vals in itertools.product(*grids):
        ps=[wpatch(a,o,v) for (a,o),v in zip(params,vals) if v!=o]
        t[tag+'_'.join(map(str,vals))]=build(base,ps)
    r=batch(t,None,timeout=900,jobs=jobs)
    return sorted(r,key=lambda q:q[2])
if __name__=='__main__':
    base=open(sys.argv[1]).read(); params=json.loads(sys.argv[2]); grids=json.loads(sys.argv[3])
    r=tune(base,params,grids,sys.argv[4])
    print('BEST',r[:3])
