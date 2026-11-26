# Parse a v1 (data model <100) .isyntax: block header table (0x2014) + seektable (0x2015) -> CSV.
# Follows Philips extract_block_header.py (read-only; reads the file, writes only to /tmp).
import sys, re, base64, struct, math, csv
path=sys.argv[1]; out=sys.argv[2]
f=open(path,'rb')
head=f.read(8*1024*1024)
while head.find(b'\r\n\x04')<0:          # the XML header of the larger files exceeds 8 MB: keep reading until the EOT marker
    more=f.read(8*1024*1024)
    if not more: raise SystemExit('no EOT marker found')
    head+=more
eot=head.find(b'\r\n\x04'); xml=head[:eot].decode('utf-8')
data_pos=eot+3
# --- WSI general info from XML (regex; no base64 decode of images)
wsi_start=xml.index('WSI')
def attr(name, s=xml, n=0):
    return re.findall(r'Name="%s"[^>]*>([^<]*)<'%name, s)
tbl=re.search(r'Name="UFS_IMAGE_BLOCK_HEADER_TABLE"[^>]*>([^<]*)<',xml).group(1)
bht=base64.b64decode(tbl)
nblocks=int(attr('UFS_IMAGE_NUMBER_OF_BLOCKS')[0])
ranges=[tuple(map(int,r.split())) for r in attr('UFS_IMAGE_DIMENSION_RANGE')]
# first 5 ranges after the general-header start belong to x,y,component,scale,waveletcoef of the WSI general header
gh=xml[xml.index('UFS_IMAGE_GENERAL_HEADERS'):]
gr=[tuple(map(int,r.split())) for r in re.findall(r'Name="UFS_IMAGE_DIMENSION_RANGE"[^>]*>([^<]*)<',gh)][:5]
nx=(gr[0][2]-gr[0][0]+gr[0][1])//gr[0][1]; ny=(gr[1][2]-gr[1][0]+gr[1][1])//gr[1][1]
orig_x,orig_y=gr[0][0],gr[1][0]
L=(gr[3][2]-gr[3][0]+gr[3][1])//gr[3][1]; max_scale=gr[3][2]
# codeblock dimension from template 0
tp=xml[xml.index('UFS_IMAGE_BLOCK_HEADER_TEMPLATES'):]
tr=[tuple(map(int,r.split())) for r in re.findall(r'Name="UFS_IMAGE_DIMENSION_RANGE"[^>]*>([^<]*)<',tp)][:2]
cbw=(tr[0][2]-tr[0][0]+tr[0][1])//tr[0][1]; cbh=(tr[1][2]-tr[1][0]+tr[1][1])//tr[1][1]
print('padded',nx,ny,'levels',L,'orig',orig_x,orig_y,'codeblock',cbw,cbh,'blocks',nblocks,'bht bytes',len(bht))
# --- block headers
pos=0; total=struct.unpack_from('<I',bht,pos)[0]; pos+=4
hdrs=[]
for _ in range(nblocks):
    tag,size=struct.unpack_from('<4sI',bht,pos); assert tag==b'\xfe\xff\x00\xe0',(pos,tag); pos+=8
    end=pos+size
    t,s=struct.unpack_from('<4sI',bht,pos); assert t==b'\x1d\x30\x0e\x20' and s==20; pos+=8
    x,y,c,sc,co=struct.unpack_from('<5I',bht,pos); pos+=20
    off=sz=0
    if size==72:
        pos+=8; off=struct.unpack_from('<Q',bht,pos)[0]; pos+=8
        pos+=8; sz=struct.unpack_from('<Q',bht,pos)[0]; pos+=8
    t,s=struct.unpack_from('<4sI',bht,pos); assert t==b'\x1d\x30\x12\x20' and s==4; pos+=8
    tid=struct.unpack_from('<I',bht,pos)[0]; pos+=4
    assert pos==end,(pos,end)
    hdrs.append([x,y,c,sc,co,off,sz,tid])
partial = all(h[5]==0 for h in hdrs)
print('headers parsed',len(hdrs),'partial (seektable needed):',partial)
# --- seektable
if partial:
    f.seek(data_pos); tag,size=struct.unpack('<4sI',f.read(8)); assert tag==b'\x1d\x30\x15\x20',tag
    n=size//40; st=[]
    raw=f.read(size)
    for i in range(n):
        p=i*40; t,s=struct.unpack_from('<4sI',raw,p); assert t==b'\xfe\xff\x00\xe0' and s==32
        o=struct.unpack_from('<Q',raw,p+16)[0]; z=struct.unpack_from('<Q',raw,p+32)[0]; st.append((o,z))
    print('seektable entries',n,'size field',size)
    base=data_pos+8+size   # data offsets in the seektable may be absolute; check below
    # block id as in Philips calculate_block_id
    tbw=int((nx+(cbw<<L)-1)/(cbw<<L))<<(L-1); tbh=int((ny+(cbh<<L)-1)/(cbh<<L))<<(L-1)
    coeff=sum((tbw*tbh)>>(2*i) for i in range(L)); ll=(tbw*tbh)>>(2*(L-1)); per_color=coeff+ll
    def fv(s): return (3<<s)-2
    for h in hdrs:
        x,y,c,sc,co=h[:5]; is_ll=(co==0)
        bid=sum((tbw>>s)*(tbh>>s) for s in range(sc))
        if is_ll: bid+=(tbw>>max_scale)*(tbh>>max_scale)
        off=fv(sc)+((1<<sc) if is_ll else 0)
        tx=int((x-orig_x-off)/((cbw*2)<<sc)); ty=int((y-orig_y-off)/((cbh*2)<<sc))
        bid+=tx+ty*(tbw>>sc)+c*per_color
        h[5],h[6]=st[bid]
    print('per-color tile count',per_color,'grid',tbw,tbh)
rows=[h for h in hdrs if h[5]!=0]
print('blocks with data:',len(rows),'of',len(hdrs))
w=csv.writer(open(out,'w')); w.writerow(['x','y','color','scale','coef','offset','size','tpl']); w.writerows(rows)
