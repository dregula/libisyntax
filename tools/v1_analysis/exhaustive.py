import sys, os, csv, collections, random, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hulsken_v1 import decode_block_v1, unpack_values
# usage: python3 -I exhaustive.py <file.isyntax> <blocks.csv from parse_v1.py>
path=sys.argv[1]
rows=list(csv.DictReader(open(sys.argv[2])))
f=open(path,'rb')
# Philips Q2 table: (Q, DZ) by (color-group, scale)
Q2={('Y',0):(1,2),('Y',1):(0,0),('C',0):(4,16),('C',1):(2,8)}
def qdz(c,s): return Q2.get(('Y' if c==0 else 'C', s),(0,0))
stats=collections.defaultdict(lambda: dict(n=0,empty8=0,decoded=0,badlen=0,viol=0,unionbits=[0]*3,minnz=None,notmult=0,nzcount=0,maxabs=0,hh_zero=0,sampled=0, bit_set_blocks=collections.Counter()))
random.seed(7)
per_group_val=6
byg=collections.defaultdict(list)
for r in rows: byg[(int(r['color']),int(r['scale']),'LL' if int(r['coef'])==0 else 'H')].append(r)
valsample=set()
for k,v in byg.items():
    nz=[r for r in v if int(r['size'])>8]
    for r in random.sample(nz,min(per_group_val,len(nz))): valsample.add(id(r))
fail=0
for idx,r in enumerate(rows):
    c,s=int(r['color']),int(r['scale']); band='LL' if int(r['coef'])==0 else 'H'; ncoef=1 if band=='LL' else 3
    st=stats[(c,s,band)]; st['n']+=1
    sz=int(r['size'])
    if sz<=8: st['empty8']+=1; continue
    f.seek(int(r['offset'])); data=f.read(sz)
    try: masks,planes,info=decode_block_v1(data,ncoef)
    except Exception as e: fail+=1; st['badlen']+=1; continue
    st['decoded']+=1
    if not info['len_ok']: st['badlen']+=1
    q,dz=qdz(c,s) if band=='H' else (0,0)
    forced=sum(1<<i for i in range(1,q+1))   # v1 mask bits 1..Q must be clear
    for k in range(ncoef):
        st['unionbits'][k]|=masks[k]
        if masks[k] & forced: st['viol']+=1
        for i in range(16):
            if (masks[k]>>i)&1: st['bit_set_blocks'][i]+=1
    if id(r) in valsample:
        vals=unpack_values(masks,planes,ncoef)
        st['sampled']+=1
        for k in range(ncoef):
            for v in vals[k]:
                if v==0: continue
                a=v&0x7FFF; st['nzcount']+=1
                if st['minnz'] is None or a<st['minnz']: st['minnz']=a
                st['maxabs']=max(st['maxabs'],a)
                if q and (a % (1<<q)): st['notmult']+=1
    if idx%5000==0: print('..',idx,file=sys.stderr)
print('failed decodes:',fail)
print('group        | blocks  <=8B  decoded  badlen | Q,DZ | blocks violating "mask bits 1..Q clear" | value check on sampled blocks: sampled nonzero-coefs  min|c|  max|c|  not-multiple-of-2^Q')
tot_viol=0; tot_bad=0
for key in sorted(stats):
    st=stats[key]; c,s,b=key; q,dz=qdz(c,s) if b=='H' else (0,0)
    tot_viol+=st['viol']; tot_bad+=st['badlen']
    print(f"c{c} s{s} {b:2s}    | {st['n']:6d} {st['empty8']:5d} {st['decoded']:7d} {st['badlen']:6d} | {q},{dz:<2d} | {st['viol']:4d} | {st['sampled']:2d} blocks {st['nzcount']:8d} {str(st['minnz']):>6s} {st['maxabs']:6d} {st['notmult']:6d}")
print('TOTAL violations of forced-clear mask bits:',tot_viol,' blocks with length mismatch:',tot_bad)
