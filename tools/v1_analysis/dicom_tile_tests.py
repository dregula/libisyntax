# DICOM-WSI tile test harness for a v1 iSyntax file: tile grids, absolute slide coordinates, per-tile exactness and
# independence, seamless stitching, sparse (blank) tile classification, for frame sizes 256, 512 and 1024.
# Read-only on the input; needs numpy (repo .venv) and ~14 GB RAM (all three channels are decoded for ground truth).
#   python3 -I dicom_tile_tests.py <file.isyntax> <blocks.csv> [D=5] [sizes=256,512,1024]
# Conventions (see iSyntax_to_DICOM_WSI_Transcoding_Claude.md):
#  * DICOM layer k = pixel plane out_k of the iSyntax decoder; sample i sits at padded level-0 position 2^k (i+3) - 2.
#  * Absolute position = ox + that (header origin, units of 0.25 um = DIMENSION_SCALE_FACTOR). Pixel spacing of layer k = 0.25 um * 2^k.
#  * Tile origin i0 = i_base + n*F with (i_base + 3) % 2^D == 0 so that every frame sits on the JPEG 2000 grid.
#  * A tile is "interior" if its pixels lie inside the decoded array: coefficients = copied interior + recomputed ring (strip algorithm).
#    Otherwise it is an "edge" tile: pixels outside the decoded array are filler (0 = black padding) and the frame is transformed locally.
import sys, os, time, csv, math, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hulsken_v1 import decode_block_v1
from windowed_idwt_halo import (open_file, unpack, need, lift_axis, fwd2, inv2, full_pyramid, inv_pyramid, Windowed, build_windowed)
UM = 0.25                                            # micrometre per level-0 unit (header: DIMENSION_SCALE_FACTOR)

def decode_channel(path, rows, g, c):
    Hp, Wp, top = g['Hp'], g['Wp'], g['top']
    H = {s: np.zeros((3, need(Hp, s), need(Wp, s)), dtype=np.int16) for s in range(g['L'])}
    LLtop = np.zeros((need(Hp, top), need(Wp, top)), dtype=np.int16); f = open(path, 'rb')
    for r in rows:
        if r['c'] != c or r['z'] <= 8: continue
        s, k = r['s'], r['k']; off = (3 << s) - 2 + ((1 << s) if k == 0 else 0)
        tx = (r['x'] - g['ox'] - off) // (256 << s); ty = (r['y'] - g['oy'] - off) // (256 << s)
        f.seek(r['o']); masks, planes, info = decode_block_v1(f.read(r['z']), 1 if k == 0 else 3)
        v = unpack(masks, planes, 1 if k == 0 else 3); tgt = LLtop[None] if k == 0 else H[s]
        y0, x0 = ty * 128, tx * 128; hh = min(128, tgt.shape[1] - y0); ww = min(128, tgt.shape[2] - x0)
        if hh > 0 and ww > 0: tgt[:, y0:y0 + hh, x0:x0 + ww] = v[:, :hh, :ww]
    cur = LLtop.astype(np.int32); outs = {}                       # full decode = ground truth only
    for s in range(top, -1, -1):
        n_r, n_c = need(Hp, s), need(Wp, s); cur = cur[:n_r, :n_c]
        out = lift_axis(lift_axis(cur, H[s][0][:n_r, :n_c], 1), lift_axis(H[s][1][:n_r, :n_c], H[s][2][:n_r, :n_c], 1), 0)
        outs[s] = out; nxt = np.zeros_like(out); nxt[1:, 1:] = out[:-1, :-1]; cur = nxt
    return H, LLtop, outs

# ------------------------------------------------------------------ geometry
def axis_grid(k, n_valid, pad, D, F):
    """Layer-k sample range overlapping the valid image on one axis, the aligned tile origin and the tile count.
    centre of sample i (padded level-0 units): p(i) = 2^k (i+3) - 2; its footprint is p +- 2^k/2; the image is [pad-0.5, pad+n_valid-0.5]."""
    p = lambda i: (1 << k) * (i + 3) - 2; half = (1 << k) / 2
    i_lo = 0
    while p(i_lo) + half <= pad - 0.5: i_lo += 1
    i_hi = i_lo
    while p(i_hi + 1) - half < pad + n_valid - 0.5: i_hi += 1
    i_base = ((i_lo + 3) // (1 << D)) * (1 << D) - 3
    cols = i_hi - i_base + 1
    return dict(i_lo=i_lo, i_hi=i_hi, i_base=i_base, n=cols, tiles=-(-cols // F))
def tile_record(k, ty, tx, F, gy, gx, ox, oy):
    i0y = gy['i_base'] + ty * F; i0x = gx['i_base'] + tx * F
    cx = ox + (1 << k) * (i0x + 3) - 2; cy = oy + (1 << k) * (i0y + 3) - 2          # centre of the top-left sample, 0.25 um units
    half = (1 << k) / 2
    return dict(i0y=i0y, i0x=i0x, col=1 + tx * F, row=1 + ty * F, cx_mm=cx * UM / 1000, cy_mm=cy * UM / 1000,
                cornerx_mm=(cx - half) * UM / 1000, cornery_mm=(cy - half) * UM / 1000, spacing_um=UM * (1 << k))

def is_interior(dec, k, rec, F):
    R, C = dec.shape(k); return rec['i0y'] >= 0 and rec['i0x'] >= 0 and rec['i0y'] + F <= R and rec['i0x'] + F <= C

def frame_from_decode(dec, k, rec, F):
    """Pixels of an edge tile: decoded part (windowed) + zero filler outside the decoded array."""
    R, C = dec.shape(k); fr = np.zeros((F, F), dtype=np.int32)
    y0, x0 = max(rec['i0y'], 0), max(rec['i0x'], 0); y1, x1 = min(rec['i0y'] + F, R), min(rec['i0x'] + F, C)
    if y1 > y0 and x1 > x0: fr[y0 - rec['i0y']:y1 - rec['i0y'], x0 - rec['i0x']:x1 - rec['i0x']] = dec.out_region(k, y0, y1, x0, x1)
    return fr, (F * F - (y1 - y0) * (x1 - x0)) / (F * F)

def tile_coeffs(decs, outs, k, rec, F, D, rng, h=2):
    """Per channel: JPEG 2000 coefficients of one tile, from the windowed decoder only. Returns [(LLD, levels, pixels_or_None)]."""
    res = []
    for c in range(3):
        if is_interior(decs[c], k, rec, F):
            LL, lv, _ = build_windowed(decs[c], k, rec['i0y'], rec['i0x'], F, D, h, rng); res.append((LL, lv))
        else:
            fr, _ = frame_from_decode(decs[c], k, rec, F); LL, lv = full_pyramid(fr, D); res.append((LL, lv))
    return res
def tile_pixels_truth(decs, outs, k, rec, F):
    return [frame_from_decode(decs[c], k, rec, F)[0] if not is_interior(decs[c], k, rec, F) else outs[c][k][rec['i0y']:rec['i0y'] + F, rec['i0x']:rec['i0x'] + F] for c in range(3)]

# ------------------------------------------------------------------ sparse classification
def classify(decs, k, rec, F, D):
    """'edge' (partly outside the decoded array), 'present', or 'blank?' (all stored detail coefficients of the tile are zero
    in all channels and the last stored LL band is constant: candidate for omission, verified separately)."""
    if not all(is_interior(decs[c], k, rec, F) for c in range(3)): return 'edge', None
    qy, qx = rec['i0y'] + 1, rec['i0x'] + 1; bg = []
    for c in range(3):
        dec = decs[c]; last = None
        for j in range(D):
            s = k + j
            if s > dec.top: break
            nn = F >> (j + 1); hx = ((qx + 2) >> (j + 1)) - 1; hy = ((qy + 2) >> (j + 1)) - 1; last = (j, s, nn, hy, hx)
            if dec.H[s][0][hy - 1:hy - 1 + nn, hx:hx + nn].any() or dec.H[s][1][hy:hy + nn, hx - 1:hx - 1 + nn].any() or dec.H[s][2][hy:hy + nn, hx:hx + nn].any():
                return 'present', None
        if last is None: return 'present', None
        j, s, nn, hy, hx = last
        ll = dec.cur_region(s, hy - 1, hy - 1 + nn, hx - 1, hx - 1 + nn)
        if ll.min() != ll.max(): return 'present', None
        bg.append(int(ll.flat[0]))
    return 'blank?', tuple(bg)

# ------------------------------------------------------------------ main
if __name__ == '__main__':
    path, csvp = sys.argv[1], sys.argv[2]; D = int(sys.argv[3]) if len(sys.argv) > 3 else 5
    SIZES = [int(v) for v in sys.argv[4].split(',')] if len(sys.argv) > 4 else [256, 512, 1024]
    t0 = time.time(); g, rows = open_file(path, csvp); top, pad, ox, oy = g['top'], g['pad'], g['ox'], g['oy']
    Wv, Hv = g['Wp'] - 2 * pad, g['Hp'] - 2 * pad
    blockx = {}                                       # header coordinates of every block, for the anchor test
    for r in rows:
        off = (3 << r['s']) - 2 + ((1 << r['s']) if r['k'] == 0 else 0)
        blockx[(r['c'], r['s'], r['k'], (r['x'] - ox - off) // (256 << r['s']), (r['y'] - oy - off) // (256 << r['s']))] = (r['x'], r['y'])
    Hs, LLs, outs, decs = {}, {}, {}, {}
    for c in range(3):
        Hs[c], LLs[c], o = decode_channel(path, rows, g, c); outs[c] = o; decs[c] = Windowed(Hs[c], LLs[c], top, g['Hp'], g['Wp'])
        print('decoded channel %d (%.0fs)' % (c, time.time() - t0), flush=True)
    print('valid image %d x %d level-0 px (%.2f x %.2f mm), origin (%.4f, %.4f) mm, level-0 spacing %.2f um' % (
        Wv, Hv, Wv * UM / 1000, Hv * UM / 1000, ox * UM / 1000, oy * UM / 1000, UM))
    print('padding region of out_0: Y at padded (5,5)=%d, (100,100)=%d, first valid pixel (%d,%d)=%d' % (
        outs[0][0][4, 4], outs[0][0][99, 99], pad, pad, outs[0][0][pad - 1, pad - 1]))
    rng = np.random.default_rng(11); allok = True
    # level-0 padded positions of the categories used for sampling (from the earlier tests: tile (60,40), (35,35), (10,10), (97,50))
    cats = {'tissue': (15358 + 128, 10238 + 128), 'tissue edge': (8958 + 128, 8958 + 128), 'blank glass': (2558 + 128, 2558 + 128),
            'scan edge': (24830 + 128, 12798 + 128)}
    for F in SIZES:
        print('\n' + '=' * 110 + '\nFRAME / TILE SIZE %d x %d  (D = %d decomposition levels, 2^D = %d)' % (F, F, D, 1 << D))
        for k in range(top + 1):
            gy = axis_grid(k, Hv, pad, D, F); gx = axis_grid(k, Wv, pad, D, F); R, C = decs[0].shape(k)
            o = tile_record(k, 0, 0, F, gy, gx, ox, oy)
            print('\n-- layer %d (spacing %.2f um): total pixel matrix %d x %d (cols x rows), %d x %d tiles = %d, matrix origin index (%d,%d), '
                  'tile(0,0) centre (%.4f, %.4f) mm' % (k, o['spacing_um'], gx['n'], gy['n'], gx['tiles'], gy['tiles'], gx['tiles'] * gy['tiles'],
                  gx['i_base'], gy['i_base'], o['cx_mm'], o['cy_mm']))
            # ---- sparse classification of all tiles
            cnt = {'edge': 0, 'present': 0, 'blank?': 0}; blanks = []; bgs = set()
            for ty in range(gy['tiles']):
                for tx in range(gx['tiles']):
                    cl, bg = classify(decs, k, tile_record(k, ty, tx, F, gy, gx, ox, oy), F, D); cnt[cl] += 1
                    if cl == 'blank?': blanks.append((ty, tx)); bgs.add(bg)
            ver = 0; nver = 0
            for ty, tx in [blanks[i] for i in rng.permutation(len(blanks))[:12]]:
                rec = tile_record(k, ty, tx, F, gy, gx, ox, oy); nver += 1
                ver += all((lambda a: a.min() == a.max())(decs[c].out_region(k, rec['i0y'], rec['i0y'] + F, rec['i0x'], rec['i0x'] + F)) for c in range(3))
            print('   sparse: %d edge, %d present, %d blank-candidates (all details zero, LL constant); %d/%d sampled candidates are pixel-constant; '
                  'distinct (Y,Co,Cg) LL values among candidates: %s' % (cnt['edge'], cnt['present'], cnt['blank?'], ver, nver, sorted(bgs)[:4]))
            # ---- sampled tiles
            picks = {}
            def tile_of(pp):
                i = ((pp[0] + 2) >> k) - 3; j = ((pp[1] + 2) >> k) - 3
                return min(max((j - gy['i_base']) // F, 0), gy['tiles'] - 1), min(max((i - gx['i_base']) // F, 0), gx['tiles'] - 1)
            picks['origin tile'] = (0, 0); picks['last tile'] = (gy['tiles'] - 1, gx['tiles'] - 1)
            picks['top edge'] = (0, gx['tiles'] // 2); picks['left edge'] = (gy['tiles'] // 2, 0)
            for nm, pp in cats.items(): picks[nm] = tile_of(pp)
            seen = {}
            for nm, t in picks.items(): seen.setdefault(t, nm)
            for (ty, tx), nm in seen.items():
                rec = tile_record(k, ty, tx, F, gy, gx, ox, oy); interior = all(is_interior(decs[c], k, rec, F) for c in range(3))
                truth = tile_pixels_truth(decs, outs, k, rec, F); line = []
                for h in ((2, 4) if interior else (None,)):
                    cf = tile_coeffs(decs, outs, k, rec, F, D, rng, h if h else 2); ok = True; nbad = 0
                    for c in range(3):
                        refLL, refLV = full_pyramid(truth[c], D); LL, lv = cf[c]
                        nbad += int((LL != refLL).sum()) + sum(int((a != b).sum()) for l1, l2 in zip(lv, refLV) for a, b in zip(l1, l2))
                        ok &= bool((inv_pyramid(LL, lv) == truth[c]).all())          # decoded from the tile's OWN coefficients only
                    line.append(('h=%d: ' % h if h else '') + '%d coeff mismatches, standalone decode exact %s' % (nbad, ok))
                    if nbad or not ok: allok = False
                # absolute coordinates: anchor to the file's block headers + cross-layer check
                anc = 'n/a'; xl = 'n/a'
                if interior:
                    qy, qx = rec['i0y'] + 1, rec['i0x'] + 1; nh_x = ((qx + 2) >> 1) - 1; nh_y = ((qy + 2) >> 1) - 1     # first high-pass coeff, level 1
                    key = (0, k, 1, nh_x // 128, nh_y // 128); hx_, hy_ = blockx[key]
                    px = ox + (3 << k) - 2 + (1 << (k + 1)) * nh_x; py = oy + (3 << k) - 2 + (1 << (k + 1)) * nh_y        # abs. position of that coefficient
                    anc = 'ok' if (hx_ + (1 << (k + 1)) * (nh_x - 128 * (nh_x // 128)) == px and hy_ + (1 << (k + 1)) * (nh_y - 128 * (nh_y // 128)) == py
                                   and px == ox + (1 << k) * (2 * nh_x + 3) - 2) else 'MISMATCH'
                    if anc != 'ok': allok = False
                    if k < top:
                        LL1 = fwd2(truth[0])[0]; i2y = (rec['i0y'] - 3) // 2; i2x = (rec['i0x'] - 3) // 2; m = F // 2; d = 4
                        bad = int((LL1[d:m - d, d:m - d] != outs[0][k + 1][i2y + d:i2y + m - d, i2x + d:i2x + m - d]).sum())
                        xl = 'ok' if bad == 0 else 'MISMATCH(%d)' % bad
                        if bad: allok = False
                print('   tile (%d,%d) [%s]%s col,row=(%d,%d) top-left sample centre (%.5f, %.5f) mm, corner (%.5f, %.5f) mm | %s | header-anchor %s, cross-layer LL %s' % (
                    ty, tx, nm, '' if interior else ' EDGE (%.0f%% filler)' % (100 * frame_from_decode(decs[0], k, rec, F)[1]), rec['col'], rec['row'],
                    rec['cx_mm'], rec['cy_mm'], rec['cornerx_mm'], rec['cornery_mm'], ' | '.join(line), anc, xl), flush=True)
            # ---- seamless stitching of independently decoded tiles (2 x 2 around the tissue position)
            ty0, tx0 = tile_of(cats['tissue']); ty0 = min(ty0, gy['tiles'] - 2); tx0 = min(tx0, gx['tiles'] - 2)
            if ty0 >= 0 and tx0 >= 0 and gy['tiles'] >= 2 and gx['tiles'] >= 2:
                bad = 0
                for c in range(3):
                    canvas = np.zeros((2 * F, 2 * F), dtype=np.int32)
                    for dy in (0, 1):
                        for dx in (0, 1):
                            rec = tile_record(k, ty0 + dy, tx0 + dx, F, gy, gx, ox, oy); LL, lv = tile_coeffs(decs, outs, k, rec, F, D, rng)[c]
                            canvas[dy * F:(dy + 1) * F, dx * F:(dx + 1) * F] = inv_pyramid(LL, lv)           # each tile decoded from its own data
                    r0 = tile_record(k, ty0, tx0, F, gy, gx, ox, oy)
                    if all(is_interior(decs[c], k, r0, 2 * F) for _ in (0,)): bad += int((canvas != outs[c][k][r0['i0y']:r0['i0y'] + 2 * F, r0['i0x']:r0['i0x'] + 2 * F]).sum())
                    else: bad += -1000000
                print('   2x2 independently decoded tiles stitched vs iSyntax decode: %s' % ('%d mismatching samples (3 channels)' % bad if bad >= 0 else 'skipped (touches the array edge)'))
                if bad > 0: allok = False
    print('\nALL OK:', allok, ' total %.0fs' % (time.time() - t0))
