# Statistics of the Hulsken v1 streams of every non-empty codeblock (needs numpy).
# usage: stream_stats.py <file.isyntax> <blocks.csv>
import sys, csv, collections, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hulsken_v1 import decode_block_v1
path, csvp = sys.argv[1], sys.argv[2]
rows = list(csv.DictReader(open(csvp)))
f = open(path, 'rb')
POP = np.array([bin(i).count('1') for i in range(256)], dtype=np.uint8)
G = collections.defaultdict(lambda: dict(n=0, comp=0, ser=0, planes=0, planes_full=0, zero_bytes=0, ones=collections.Counter(), stored=collections.Counter(), cnt=collections.Counter(), sym=np.zeros(256, dtype=np.int64)))
leaves = collections.Counter(); zrb = collections.Counter(); zrs = collections.Counter(); maxlen = collections.Counter()
empty_n = 0
for r in rows:
    c, s, k, z = int(r['color']), int(r['scale']), int(r['coef']), int(r['size'])
    if z <= 8: empty_n += 1; continue
    ncoef = 1 if k == 0 else 3
    f.seek(int(r['offset'])); data = f.read(z)
    masks, planes, info = decode_block_v1(data, ncoef)
    g = G[(c, s, 'LL' if k == 0 else 'H')]
    g['n'] += 1; g['comp'] += z; g['ser'] += info['ser_len']
    arr = np.frombuffer(planes, dtype=np.uint8)
    g['zero_bytes'] += int((arr == 0).sum()); g['sym'] += np.bincount(arr, minlength=256)
    leaves[info['tree_leaves']] += 1; zrb[info['zr_bits']] += 1; zrs[info['zr_sym']] += 1
    maxlen[max(d for _, d in info['codes'].values())] += 1
    p = 0
    for kk in range(ncoef):
        for i in range(16):
            g['cnt'][i] += 1                                    # chance for plane i
            if (masks[kk] >> i) & 1:
                pl = arr[p * 2048:(p + 1) * 2048]; p += 1
                g['ones'][i] += int(POP[pl].sum()); g['planes'] += 1; g['stored'][i] += 1
print('non-empty blocks decoded:', sum(g['n'] for g in G.values()), '| empty (<=8 B):', empty_n)
print('\nHuffman/zero-run parameters over all non-empty blocks:')
print('  tree leaves  min/median/max:', min(leaves), sorted(leaves.elements())[len(list(leaves.elements()))//2], max(leaves))
print('  max code length histogram:', dict(sorted(maxlen.items())))
print('  zero-run counter bits histogram:', dict(sorted(zrb.items())))
print('  zero-run symbol (most common 6):', zrs.most_common(6))
print('\nper group: blocks | raw 16-bit size MB | compressed MB | ratio raw/comp | serialized MB (after plane dropping) | zero bytes in serialized %  | entropy stage ratio ser/comp')
tot_raw = tot_c = tot_s = 0
for key in sorted(G):
    g = G[key]; c, s, kind = key; ncoef = 1 if kind == 'LL' else 3
    raw = g['n'] * ncoef * 128 * 128 * 2
    tot_raw += raw; tot_c += g['comp']; tot_s += g['ser']
    zb = g['zero_bytes'] / max(1, g['planes'] * 2048)
    print(f"c{c} s{s} {kind:2s} | {g['n']:5d} | {raw/1e6:8.1f} | {g['comp']/1e6:8.2f} | {raw/g['comp']:6.1f} | {g['ser']/1e6:8.2f} | {100*zb:5.1f}% | {g['ser']/g['comp']:5.2f}")
print(f"TOTAL non-empty: raw {tot_raw/1e6:.0f} MB -> compressed {tot_c/1e6:.1f} MB ({tot_raw/tot_c:.1f}:1); after plane dropping {tot_s/1e6:.1f} MB; entropy stage {tot_s/tot_c:.2f}:1")
print('\nbit density per stored plane (fraction of the 16384 coefficients with that bit set, averaged over blocks/subbands where the plane is stored)')
print('  sign = mask bit 0, bK = magnitude bit K (mask bit K+1); "." = plane never stored in this group; n = share of subband-blocks that store the plane')
print('group     ' + ' '.join(f'{"sign" if b < 0 else "b%d" % b:>7s}' for b in range(-1, 9)))
for key in sorted(G):
    c, s_, kind = key
    if kind == 'LL': continue
    g = G[key]; nsb = g['n'] * 3
    cells = []
    for b in range(-1, 9):
        i = b + 1
        cells.append('      .' if not g['stored'][i] else f"{g['ones'][i] / (g['stored'][i] * 16384):7.4f}")
    print(f"c{c} s{s_}   " + ' '.join(cells))
print('stored-plane share (fraction of subband-blocks that store the plane):')
for key in sorted(G):
    c, s_, kind = key
    if kind == 'LL': continue
    g = G[key]; nsb = g['n'] * 3
    print(f"c{c} s{s_}   " + ' '.join(f"{g['stored'][b + 1] / nsb:7.2f}" for b in range(-1, 9)))
