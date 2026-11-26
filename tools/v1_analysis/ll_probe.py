import sys, csv, collections
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hulsken_v1 import decode_block_v1, unpack_values
path, csvp = sys.argv[1], sys.argv[2]
ox, oy = 77715, 20614; s = 5; off = (3<<s)-2+(1<<s)
f = open(path,'rb'); planes = {0:{},1:{},2:{}}
for r in csv.DictReader(open(csvp)):
    if int(r['coef'])!=0 or int(r['scale'])!=5: continue
    c=int(r['color']); tx=(int(r['x'])-ox-off)//(256<<s); ty=(int(r['y'])-oy-off)//(256<<s)
    f.seek(int(r['offset'])); data=f.read(int(r['size']))
    masks,pl,info=decode_block_v1(data,1); assert info['len_ok']
    v=unpack_values(masks,pl,1)[0]
    out=[0]*16384
    for i in range(16384):
        area=i//16; w=i%16
        out[((area//32)*4+w//4)*128+(area%32)*4+w%4]=v[i]
    planes[c][(tx,ty)]=[(-(x&0x7fff) if x>>15 else x) for x in out]
print('LL blocks per colour:',[len(planes[c]) for c in range(3)])
W=4*128;H=3*128
img=[[ [0]*(W) for _ in range(H)] for c in range(3)]
for c in range(3):
    for (tx,ty),b in planes[c].items():
        for y in range(128):
            img[c][ty*128+y][tx*128:(tx+1)*128]=b[y*128:(y+1)*128]
Y,Co,Cg=img
pts=[(Y[y][x],Co[y][x],Cg[y][x]) for y in range(H) for x in range(W) if Y[y][x]!=0]
print('LL grid',W,'x',H,' nonzero-Y samples',len(pts))
for n,k in (('Y',0),('Co',1),('Cg',2)):
    vs=[p[k] for p in pts]; print(' %s min %d max %d'%(n,min(vs),max(vs)))
fix=lambda x: x//2 if x>=0 else -((-x)//2)
fl=lambda x: x>>1
def inv(Y,Co,Cg,h):
    t=Y-h(Cg); G=t+Cg; B=t-h(Co); return B+Co,G,B
n=len(pts)
print('\nparity of (Co,Cg): counts', dict(collections.Counter((p[1]&1,p[2]&1) for p in pts)))
print('joint (Y,Co,Cg) mod 2:', dict(sorted(collections.Counter((p[0]&1,p[1]&1,p[2]&1) for p in pts).items())))
dif=sum(1 for p in pts if inv(*p,fix)!=inv(*p,fl)); print('\nsamples where trunc and floor inverse differ:',dif,'of',n)
for nm,h in (('trunc',fix),('floor',fl)):
    lo=sum(1 for p in pts if min(inv(*p,h))<0); hi=sum(1 for p in pts if max(inv(*p,h))>255)
    print(' %s inverse: RGB<0 in %d samples, >255 in %d samples'%(nm,lo,hi))
print('non-gray samples (Co or Cg !=0):',sum(1 for p in pts if p[1] or p[2]))
