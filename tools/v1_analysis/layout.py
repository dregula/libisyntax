# Layout of codeblocks on disk: chunk (cluster) structure from real offsets. usage: layout.py <blocks.csv>
import sys, csv, collections
rows=[dict(x=int(r['x']),y=int(r['y']),c=int(r['color']),s=int(r['scale']),k=int(r['coef']),o=int(r['offset']),z=int(r['size'])) for r in csv.DictReader(open(sys.argv[1]))]
rows.sort(key=lambda r:r['o'])
# contiguity: every block starts where the previous one ended?
gaps=[rows[i+1]['o']-(rows[i]['o']+rows[i]['z']) for i in range(len(rows)-1)]
print('blocks',len(rows),'| gaps between consecutive blocks (bytes): min',min(gaps),'max',max(gaps),'count nonzero',sum(1 for g in gaps if g))
print('gap histogram (first 8):',collections.Counter(gaps).most_common(8))
# chunks = runs separated by a nonzero gap (or by the 16-byte DICOM item header pattern)
chunks=[[rows[0]]]
for i,g in enumerate(gaps):
    if g!=0 and g!=8 and g!=16: chunks.append([])
    chunks[-1].append(rows[i+1])
print('runs split on gaps other than 0/8/16:',len(chunks))
# composition of blocks per scale in consecutive groups of 63/66 (check chunk shape)
def comp(ch): return tuple(sorted(collections.Counter((r['s'],'LL' if r['k']==0 else 'H') for r in ch).items()))
cc=collections.Counter(comp(c) for c in chunks)
for k,v in cc.most_common(6): print(v,'chunks with composition',k)
# order inside a chunk: print the first full chunk in file order
c=chunks[len(chunks)//2]
print('\nexample chunk in FILE ORDER (x y color scale coef size):')
for r in c[:26]: print(' ',r['x'],r['y'],r['c'],r['s'],r['k'],r['z'])
print('  ... (%d blocks in this chunk)'%len(c))
