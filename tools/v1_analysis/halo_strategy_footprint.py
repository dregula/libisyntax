# Index-only comparison of halo-correction strategies for the iSyntax -> DICOM (JPEG 2000 5/3) transcoder. No pixel or coefficient decoding:
# it replays the exact window arithmetic of windowed_idwt_halo.py (build_windowed / Windowed) and only records which codeblocks are touched and how many
# samples the lifting and the strip transforms would process, then totals them over EVERY frame of every layer using the block table (blocks.csv).
#   python3 -B -I halo_strategy_footprint.py <file.isyntax> <blocks.csv> [F=256] [D=5] [h=4]
# Strategies (per colour; all three colours touch the same blocks, empty blocks (size <= 8 B) cost nothing):
#   THIN   frame-independent task, thin ring strips (width 2h+4) from a windowed top-down IDWT, interior H copied from the stored blocks.
#          'no cache' = every frame fetches and decodes its own footprint; 'LRU C' = same with a shared cache of C blocks, frames in raster order.
#   RECRUIT same interior copy, but the ring strips come from WHOLE neighbouring iSyntax tiles as libisyntax_tile_read() does it: at every scale from the
#          frame's own down to the top, the tiles under the frame footprint (ancestors above the frame's own scales) plus their 8 neighbours are loaded,
#          and every footprint tile is inverse-transformed in full (2 x (128 + 8))^2 samples per tile and colour (ISYNTAX_IDWT_PAD_L/R = 4).
#   STRIP  the planned pipeline: decode every block once; 'per layer' = each layer runs its own sliding window (a block of scale s is decoded again for each
#          layer that uses it), 'shared' = all layers advance together and every block is decoded once.
import sys, os, csv, time, collections
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from windowed_idwt_halo import open_file
from dicom_tile_tests import axis_grid, tile_record
need = lambda n, s: -(-n // (1 << (s + 1))) + 2
path, csvp = sys.argv[1], sys.argv[2]
F = int(sys.argv[3]) if len(sys.argv) > 3 else 256; D = int(sys.argv[4]) if len(sys.argv) > 4 else 5; h = int(sys.argv[5]) if len(sys.argv) > 5 else 4
g, rows = open_file(path, csvp); top, pad, Hp, Wp = g['top'], g['pad'], g['Hp'], g['Wp']; Wv, Hv = Wp - 2 * pad, Hp - 2 * pad
SZ = {}                                              # (c, s, kind, tx, ty) -> compressed bytes (kind 0 = LL, 1 = H)
for r in rows:
    off = (3 << r['s']) - 2 + ((1 << r['s']) if r['k'] == 0 else 0)
    SZ[(r['c'], r['s'], r['k'], (r['x'] - g['ox'] - off) // (256 << r['s']), (r['y'] - g['oy'] - off) // (256 << r['s']))] = r['z']
NB = {s: (-(-need(Hp, s) // 128), -(-need(Wp, s) // 128)) for s in range(top + 1)}   # block rows, cols per scale
def cost(blk):                                       # (decodes, bytes) over the three colours for one touched block (s, kind, bc, br)
    s, kind, bc, br = blk; n = b = 0
    for c in range(3):
        z = SZ.get((c, s, 1 if kind == 'H' else 0, bc, br), 0)
        if z > 8: n += 1; b += z
    return n, b
class Sim:
    def __init__(s_, top, Hp, Wp):
        s_.top = top; s_.lim = {s: (need(Hp, s), need(Wp, s)) for s in range(top + 1)}; s_.touched = set(); s_.lift = 0; s_.rects = []
    def touch(s_, s, kind, r0, r1, c0, c1):
        for br in range(r0 // 128, (r1 - 1) // 128 + 1):
            for bc in range(c0 // 128, (c1 - 1) // 128 + 1): s_.touched.add((s, kind, bc, br))
    def cur_region(s_, s, r0, r1, c0, c1):
        assert 0 <= r0 < r1 <= s_.lim[s][0] and 0 <= c0 < c1 <= s_.lim[s][1]
        if s == s_.top: s_.touch(s, 'LL', r0, r1, c0, c1); return
        pr = 1 if r0 == 0 else 0; pc = 1 if c0 == 0 else 0
        if r1 - pr <= 0 or c1 - pc <= 0: return
        s_.out_region(s + 1, r0 - 1 + pr, r1 - 1, c0 - 1 + pc, c1 - 1)
    def out_region(s_, s, r0, r1, c0, c1):
        R, C = s_.lim[s]; assert 0 <= r0 < r1 <= 2 * R and 0 <= c0 < c1 <= 2 * C
        nr0 = max(0, (r0 >> 1) - 1); nr1 = min(R, ((r1 - 1) >> 1) + 3); nc0 = max(0, (c0 >> 1) - 1); nc1 = min(C, ((c1 - 1) >> 1) + 3)
        s_.cur_region(s, nr0, nr1, nc0, nc1); s_.touch(s, 'H', nr0, nr1, nc0, nc1); s_.lift += 4 * (nr1 - nr0) * (nc1 - nc0); s_.rects.append((s, nr0, nr1, nc0, nc1))
def rect_fill(n, getter, W):
    if n <= 2 * W: getter(0, n, 0, n); return
    getter(0, n, 0, W); getter(0, n, n - W, n); getter(0, W, 0, n); getter(n - W, n, 0, n)
def thin(dec, k, i0y, i0x):
    """replays build_windowed(): returns forward-transform samples of the strips; dec.touched / dec.lift accumulate"""
    W = 2 * h + 4; qy = i0y + 1; qx = i0x + 1; fwd = 0
    rect_fill(F, lambda r0, r1, c0, c1: dec.out_region(k, i0y + r0, i0y + r1, i0x + c0, i0x + c1), W); n = F
    for j in range(D):
        s = k + j; nn = n // 2
        if s > dec.top: fwd += n * n; n = nn; continue
        hx = ((qx + 2) >> (j + 1)) - 1; hy = ((qy + 2) >> (j + 1)) - 1
        dec.touch(s, 'H', hy - 1, hy + nn - 1, hx, hx + nn); dec.touch(s, 'H', hy, hy + nn, hx - 1, hx + nn - 1); dec.touch(s, 'H', hy, hy + nn, hx, hx + nn)
        by = ((qy + 2) >> (j + 1)) - 2; bx = ((qx + 2) >> (j + 1)) - 2
        if j == D - 1 or s == dec.top or nn <= 2 * W: dec.cur_region(s, by, by + nn, bx, bx + nn)
        else: rect_fill(nn, lambda r0, r1, c0, c1: dec.cur_region(s, by + r0, by + r1, bx + c0, bx + c1), W)
        fwd += n * n if nn <= 2 * h else 4 * W * n
        n = nn
    return fwd
def recruit(k, i0y, i0x):
    """whole-tile variant: tiles under the frame footprint at every scale k..top (+ their 8 neighbours); returns (blocks touched, tiles inverse-transformed)"""
    qy = i0y + 1; qx = i0x + 1; blocks = set(); tiles = 0
    lo_y = ((qy + 2) >> 1) - 2; hi_y = lo_y + F // 2 + 2; lo_x = ((qx + 2) >> 1) - 2; hi_x = lo_x + F // 2 + 2     # coefficient range at the frame's own scale k
    for s in range(k, top + 1):
        if s > k: lo_y >>= 1; hi_y = (hi_y >> 1) + 1; lo_x >>= 1; hi_x = (hi_x >> 1) + 1                         # ancestors
        lo_y = max(lo_y, 0); lo_x = max(lo_x, 0); nbr, nbc = NB[s]
        ty0, ty1, tx0, tx1 = lo_y // 128, min(nbr - 1, (hi_y - 1) // 128), lo_x // 128, min(nbc - 1, (hi_x - 1) // 128)
        tiles += (ty1 - ty0 + 1) * (tx1 - tx0 + 1)                                                               # tiles whose IDWT is run
        for ty in range(max(0, ty0 - 1), min(nbr - 1, ty1 + 1) + 1):
            for tx in range(max(0, tx0 - 1), min(nbc - 1, tx1 + 1) + 1):
                blocks.add((s, 'H', tx, ty))
                if s == top: blocks.add((s, 'LL', tx, ty))
    return blocks, tiles
class LRU:
    def __init__(s_, cap): s_.cap = cap; s_.d = collections.OrderedDict(); s_.dec = 0; s_.bytes = 0
    def get(s_, blk, nb):
        n, b = nb
        if blk in s_.d: s_.d.move_to_end(blk); return
        s_.dec += n; s_.bytes += b
        if s_.cap > 0: s_.d[blk] = 1
        while len(s_.d) > s_.cap: s_.d.popitem(last=False)
caps = [0, 64, 256, 1024, 4096]
print('file %s: padded %dx%d, %d levels (top scale %d), F=%d D=%d h=%d' % (os.path.basename(path), Wp, Hp, g['L'], top, F, D, h), flush=True)
tot = collections.defaultdict(float); uniq_all = set(); t0 = time.time()
hdr = 'layer  frames  | THIN per frame: blocks(med/p95) decodes bytes(KB) lifted-Ksamples fwd-Ksamples | RECRUIT per frame: blocks decodes bytes(KB) tiles | unique decodes (layer) | decodes with LRU ' + ' '.join('C=%d' % c for c in caps)
print(hdr)
for k in range(top + 1):
    gy, gx = axis_grid(k, Hv, pad, D, F), axis_grid(k, Wv, pad, D, F); dec0 = Sim(top, Hp, Wp)
    R, C = 2 * dec0.lim[k][0], 2 * dec0.lim[k][1]
    UNION = os.environ.get('UNION', '1') == '1'; allrects = collections.defaultdict(list); thin_n = []; thin_b = []; thin_l = []; thin_f = []; rec_n = []; rec_b = []; rec_t = []; skipped = 0; uniq = set(); lrus = [LRU(c) for c in caps]; rlru = [LRU(c) for c in caps]
    for ty in range(gy['tiles']):
        for tx in range(gx['tiles']):
            rec = tile_record(k, ty, tx, F, gy, gx, g['ox'], g['oy'])
            if rec['i0y'] < 8 or rec['i0x'] < 8 or rec['i0y'] + F > R - 8 or rec['i0x'] + F > C - 8: skipped += 1; continue
            d = Sim(top, Hp, Wp)
            try: fwd = thin(d, k, rec['i0y'], rec['i0x'])
            except AssertionError: skipped += 1; continue
            if UNION:
                for (rs, a0, a1, b0, b1) in d.rects: allrects[rs].append((a0, a1, b0, b1))
            nbs = {b: cost(b) for b in d.touched}; thin_n.append(len(d.touched)); thin_b.append(sum(v[1] for v in nbs.values())); thin_l.append(d.lift); thin_f.append(fwd)
            tn = sum(v[0] for v in nbs.values()); thin_n[-1] = tn
            for b in sorted(d.touched, key=lambda b: (b[3], b[2], b[0])):
                uniq.add(b); [l.get(b, nbs[b]) for l in lrus]
            rb, rt = recruit(k, rec['i0y'], rec['i0x']); rn = {b: cost(b) for b in rb}; rec_n.append(sum(v[0] for v in rn.values())); rec_b.append(sum(v[1] for v in rn.values())); rec_t.append(rt)
            for b in sorted(rb, key=lambda b: (b[3], b[2], b[0])): [l.get(b, rn[b]) for l in rlru]
    uniq_all |= {(k,) + b for b in uniq} if False else set(uniq)
    ulift = 0
    if UNION:
        for rs, rl in allrects.items():
            R_, C_ = dec0.lim[rs]; diff = np.zeros((R_ + 2, C_ + 2), dtype=np.int32)
            for a0, a1, b0, b1 in rl: diff[a0, b0] += 1; diff[a0, b1] -= 1; diff[a1, b0] -= 1; diff[a1, b1] += 1
            ulift += 4 * int((np.cumsum(np.cumsum(diff, 0), 1) > 0).sum()); del diff
        tot['union_lift'] += ulift
        print('   layer %d: lifted samples per colour, summed over frames %.1f M, union (each cell once) %.1f M  -> sharing saves %.0f%%' % (k, sum(thin_l) / 1e6, ulift / 1e6, 100 * (1 - ulift / max(1, sum(thin_l)))), flush=True)
    nf = len(thin_n); u_dec = sum(cost(b)[0] for b in uniq); u_b = sum(cost(b)[1] for b in uniq)
    if nf == 0: print(k, 'no interior frames'); continue
    tn = np.array(thin_n); print('%d  %6d (%d edge frames skipped) | THIN %5.1f / %5.1f  %5.1f dec  %7.0f KB  %6.0f  %6.0f | RECRUIT %6.1f dec %8.0f KB %5.1f tiles | unique %6d dec %6.2f MB | THIN-LRU %s | RECRUIT-LRU %s'
          % (k, nf, skipped, np.median(tn), np.percentile(tn, 95), tn.mean(), np.mean(thin_b) / 1e3, np.mean(thin_l) / 1e3, np.mean(thin_f) / 1e3, np.mean(rec_n), np.mean(rec_b) / 1e3, np.mean(rec_t), u_dec, u_b / 1e6,
             ' '.join('%d' % l.dec for l in lrus), ' '.join('%d' % l.dec for l in rlru)), flush=True)
    tot['frames'] += nf; tot['thin_dec'] += sum(thin_n); tot['thin_bytes'] += sum(thin_b); tot['thin_lift'] += sum(thin_l); tot['thin_fwd'] += sum(thin_f)
    tot['rec_dec'] += sum(rec_n); tot['rec_bytes'] += sum(rec_b); tot['rec_tiles'] += sum(rec_t); tot['uniq_dec_layer'] += u_dec; tot['uniq_bytes_layer'] += u_b
    for i, c in enumerate(caps): tot['thin_lru%d' % c] += lrus[i].dec; tot['thin_lru_b%d' % c] += lrus[i].bytes; tot['rec_lru%d' % c] += rlru[i].dec; tot['rec_lru_b%d' % c] += rlru[i].bytes
    uniq_all |= {(k,) + b for b in []}
allb = {(s, 'H' if kk == 1 else 'LL', tx, ty) for (c, s, kk, tx, ty), z in SZ.items() if z > 8 for c in [c]}      # every non-empty (scale, kind, block) once
shared_dec = sum(1 for (c, s, kk, tx, ty), z in SZ.items() if z > 8); shared_bytes = sum(z for z in SZ.values() if z > 8)
print('\nTOTALS over %d interior frames of all %d layers (decodes = non-empty codeblock decodes over the three colours):' % (tot['frames'], top + 1))
print('  STRIP shared (all layers in lockstep, every block once): %d decodes, %.1f MB read (= the file)' % (shared_dec, shared_bytes / 1e6))
print('  STRIP per layer (each layer own window)                : %d decodes, %.1f MB read  (%.2fx shared)' % (tot['uniq_dec_layer'], tot['uniq_bytes_layer'] / 1e6, tot['uniq_dec_layer'] / shared_dec))
print('  THIN, per-frame task, no cache                           : %d decodes (%.1fx shared), %.1f MB read (%.1fx file)' % (tot['thin_dec'], tot['thin_dec'] / shared_dec, tot['thin_bytes'] / 1e6, tot['thin_bytes'] / shared_bytes))
for c in caps[1:]: print('  THIN with LRU cache of %5d blocks                       : %d decodes (%.2fx shared), %.1f MB read (%.2fx file)' % (c, tot['thin_lru%d' % c], tot['thin_lru%d' % c] / shared_dec, tot['thin_lru_b%d' % c] / 1e6, tot['thin_lru_b%d' % c] / shared_bytes))
print('  RECRUIT whole tiles, no cache                            : %d decodes (%.1fx shared), %.1f MB read (%.1fx file)' % (tot['rec_dec'], tot['rec_dec'] / shared_dec, tot['rec_bytes'] / 1e6, tot['rec_bytes'] / shared_bytes))
for c in caps[1:]: print('  RECRUIT with LRU cache of %5d blocks                    : %d decodes (%.2fx shared), %.1f MB read (%.2fx file)' % (c, tot['rec_lru%d' % c], tot['rec_lru%d' % c] / shared_dec, tot['rec_lru_b%d' % c] / 1e6, tot['rec_lru_b%d' % c] / shared_bytes))
print('  transform work per colour, THIN: lifted %.1f M output samples + strip forward %.1f M samples ; RECRUIT: %d tile inverse transforms = %.1f M samples (+ the frame forward transform)' % (tot['thin_lift'] / 1e6, tot['thin_fwd'] / 1e6, tot['rec_tiles'], tot['rec_tiles'] * 2 * 136 * 2 * 136 / 1e6))
if UNION: print('  transform work per colour, THIN summed %.1f M lifted samples vs shared (union) %.1f M (saving %.0f%%)' % (tot['thin_lift'] / 1e6, tot['union_lift'] / 1e6, 100 * (1 - tot['union_lift'] / tot['thin_lift'])))
print('total %.0fs' % (time.time() - t0))
