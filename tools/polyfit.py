import numpy as np
from collections import deque
def line(p0,p1):
    x0,y0=p0;x1,y1=p1
    dx=x1-x0;dy=y1-y0;d=max(abs(dx),abs(dy))
    if d==0: return [(x1,y1)]
    c=1 if dx*dy<=0 else 0
    x=x0*d+(d-c)//2; y=y0*d+(d-c)//2
    out=[]
    for _ in range(d):
        out.append((x//d,y//d)); x+=dx; y+=dy
    out.append((x1,y1))
    return out
def raster(pts, fillpt, shape):
    H,W=shape
    m=np.zeros((H,W),bool)
    for a,b in zip(pts,pts[1:]+pts[:1]):
        for x,y in line(a,b):
            if 0<=x<W and 0<=y<H: m[y,x]=1
    q=deque([fillpt]); 
    if m[fillpt[1],fillpt[0]]: return m
    seen=np.zeros_like(m); seen[fillpt[1],fillpt[0]]=1
    while q:
        x,y=q.popleft(); m[y,x]=1
        for nx,ny in ((x+1,y),(x-1,y),(x,y+1),(x,y-1)):
            if 0<=nx<W and 0<=ny<H and not seen[ny,nx] and not m[ny,nx]:
                seen[ny,nx]=1; q.append((nx,ny))
    return m
def trace(M):
    """Moore-neighbour boundary trace, clockwise in image coords; returns list of (x,y)"""
    H,W=M.shape
    ys,xs=np.nonzero(M)
    i=np.lexsort((xs,ys))[0]; start=(xs[i],ys[i])
    # directions clockwise starting East (image coords, y down)
    D=[(1,0),(1,1),(0,1),(-1,1),(-1,0),(-1,-1),(0,-1),(1,-1)]
    def inside(x,y): return 0<=x<W and 0<=y<H and M[y,x]
    cur=start; back=4  # came from west
    path=[cur]
    while True:
        for k in range(8):
            dd=(back+1+k)%8
            nx,ny=cur[0]+D[dd][0],cur[1]+D[dd][1]
            if inside(nx,ny):
                back=(dd+4)%8
                cur=(nx,ny); break
        if cur==start and len(path)>1: break
        path.append(cur)
        if len(path)>100000: raise Exception('loop')
    return path
def simplify(path):
    n=len(path); verts=[0]; i=0
    while i<n:
        best=i+1
        j=i+1
        while j<=n:
            seg=[path[k%n] for k in range(i,j+1)]
            if line(seg[0],seg[-1])==seg: best=j; j+=1
            else:
                j+=1
                if j-best>6: break
        i=best
        verts.append(i)
    verts=verts[:-1]
    return [path[v%n] for v in verts]
