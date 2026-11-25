# Compare the v1 seektable (parse_v1.py CSV) with the v2 cluster header table of the matching .fic (header-only export).
# usage: compare_v1_v2.py <blocks.csv> <file.fic>
import sys, csv, re, base64, struct, collections
import xml.etree.ElementTree as ET
rows=[dict(x=int(r['x']),y=int(r['y']),c=int(r['color']),s=int(r['scale']),k=int(r['coef']),o=int(r['offset']),z=int(r['size'])) for r in csv.DictReader(open(sys.argv[1]))]
raw=open(sys.argv[2],'rb').read().decode('utf-8','replace'); root=ET.fromstring(raw)
def attrs(n): return {a.get('Element'):(a.text or '').strip() for a in n.findall('Attribute')}
def sub(n,el):
    for a in n.findall('Attribute'):
        if a.get('Element')==el: return a
def dos(a): return [] if a is None or a.find('Array') is None else a.find('Array').findall('DataObject')
def rng(d): return [tuple(int(x) for x in attrs(r)['0x200b'].split()) for r in dos(sub(d,'0x200a'))]
wsi=[im for im in dos(sub(root,'0x1003')) if attrs(im).get('0x1004')=='WSI'][0]
tpls=[(rng(t)[2][0],rng(t)[3][0],rng(t)[4][0]==3) for t in dos(sub(wsi,'0x2009'))]
ctpl=[[int(attrs(h)['0x2012']) for h in dos(sub(c,'0x200d'))] for c in dos(sub(wsi,'0x2016'))]
b=base64.b64decode(re.search(r'Element="0x201f"[^>]*>([^<]*)<',raw).group(1)+'===')
assert struct.unpack_from('<I',b,0)[0]==len(b)-4
pos=4; v2=[]; nclusters=0; ctids=collections.Counter()
while pos<len(b):
    g,e,sz=struct.unpack_from('<HHI',b,pos); end=pos+8+sz; q=pos+8; el={}
    while q<end:
        g,e,s=struct.unpack_from('<HHI',b,q); el[e]=b[q+8:q+8+s]; q+=8+s
    pos=end; nclusters+=1
    cx,cy,csc=struct.unpack('<3I',el[0x2019]); ct=struct.unpack('<I',el[0x201a])[0]; ctids[ct]+=1
    off=struct.unpack('<Q',el[0x201b])[0]; size=struct.unpack('<Q',el[0x201c])[0]
    boff=struct.unpack('<%dI'%(len(el[0x201d])//4),el[0x201d]); bsz=struct.unpack('<%dI'%(len(el[0x201e])//4),el[0x201e])
    for j,(o,s) in enumerate(zip(boff,bsz)):
        col,sc,ll=tpls[ctpl[ct][j]]
        v2.append((col,sc,0 if ll else 1,off+o,s))
print('v2 clusters',nclusters,dict(ctids),'v2 blocks',len(v2))
# v1 keyed by (colour, scale, coef) -> multiset of (offset,size); v2 same (coordinates not compared: only offset/size/identity)
m1=collections.Counter((r['c'],r['s'],r['k'],r['o'],r['z']) for r in rows)
m2=collections.Counter(v2)
print('identical (colour,scale,coef,offset,size) multisets:',m1==m2, '| only in v1:',sum((m1-m2).values()),'only in v2:',sum((m2-m1).values()))
