import csv, sys, collections, statistics
ox, oy, L = 77715, 20614, 6
W = {s: 128 >> s for s in range(L)}; H_ = {s: 96 >> s for s in range(L)}
size = {}
for r in csv.DictReader(open(sys.argv[1])):
    x,y,c,s,k,z = int(r['x']),int(r['y']),int(r['color']),int(r['scale']),int(r['coef']),int(r['size'])
    off = (3<<s)-2 + ((1<<s) if k==0 else 0)
    tx = (x-ox-off)//(256<<s); ty = (y-oy-off)//(256<<s)
    size[(c,s,k,tx,ty)] = z
def chunk(s,tx,ty):
    g=s//3; top=3*g+2; sh=top-s
    return (g, tx>>sh, ty>>sh)
def nb(s,tiles):
    out=set()
    for tx,ty in tiles:
        for dx in(-1,0,1):
            for dy in(-1,0,1):
                a,b=tx+dx,ty+dy
                if 0<=a<W[s] and 0<=b<H_[s]: out.add((a,b))
    return out
def closure(tx,ty,full):
    S={(tx,ty)}; chunks=set(); blocks=0; byts=0
    for s in range(L):
        N=nb(s,S) if (full or s==0) else set()
        for (a,b) in N:
            chunks.add(chunk(s,a,b))
            for c in range(3):
                blocks+=1; byts+=size.get((c,s,1,a,b),0)
        if not full: break
        if s==L-1:
            for (a,b) in N:
                chunks.add(chunk(s,a,b))   # top LL blocks live in the same top chunk
                for c in range(3):
                    blocks+=1; byts+=size.get((c,s,0,a,b),0)
        else:
            S={(a//2,b//2) for a,b in N}
    return chunks,blocks,byts
chunk_bytes=collections.Counter()
for (c,s,k,tx,ty),z in size.items(): chunk_bytes[chunk(s,tx,ty)]+=z
for full in (False,True):
    cnt=collections.Counter(); g0=collections.Counter(); bl=[]; by=[]; cb=[]
    for ty in range(96):
        for tx in range(128):
            ch,b,y=closure(tx,ty,full)
            cnt[len(ch)]+=1; bl.append(b); by.append(y); cb.append(sum(chunk_bytes[c] for c in ch))
            g0[sum(1 for c in ch if c[0]==0)]+=1
    print('FULL closure (H + LL chain to top)' if full else 'SAME-SCALE H neighbours only (scale 0)')
    print(' distinct chunks per tile:',dict(sorted(cnt.items())))
    print(' of which scale0-2 chunks:',dict(sorted(g0.items())))
    print(' codeblocks (3 colours): min %d median %d max %d'%(min(bl),statistics.median(bl),max(bl)))
    print(' codeblock bytes: median %d max %d | whole-chunk bytes: median %d max %d'%(statistics.median(by),max(by),statistics.median(cb),max(cb)))
