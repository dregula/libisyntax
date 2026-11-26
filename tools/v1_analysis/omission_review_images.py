# Visual-review material for tile-omission strategies (read-only on the iSyntax file; needs numpy + Pillow; reads the CSVs of omission_size_estimate.py).
#   python3 -B -I omission_review_images.py <file.isyntax> <blocks.csv> [dir=../transcoding_review] [F=256]
# Strategies (layer-0 tile grid of frame size F):
#   A      : omit only exact filler tiles (YCoCg = 255, 0, 0)
#   R<T>   : omit filler + glass tiles with robust RGB range (99.5 %-0.5 % percentile, worst channel) <= T          (default T = 4)
#   S<S>   : omit filler + glass tiles with strict RGB range (max - min, worst channel) <= S, i.e. no speck beyond S levels (default S = 16)
# Output (in dir): overview_original_L3.png (layer-3 image = 1/8 resolution), masks_A.png, masks_R4.png, masks_S16.png (omitted tiles tinted),
# overview_R4_singlefill.png / overview_S16_singlefill.png (omitted tiles replaced by ONE fill colour, as a viewer would show them) with matching
# diff_*_x16.png (|difference| amplified 16x), and crops_*.png: 3 x 3 tile crops at full resolution (original | omitted | |diff| x 16) for borderline cases.
# The single fill colour is the median colour of the omitted GLASS tiles; the filler bands therefore also take that colour in the singlefill images.
import sys, os, csv, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from PIL import Image
from colour_ebcot_probe import decode_coeffs, inv_rgb, fl
from windowed_idwt_halo import open_file, need, lift_axis, Windowed
from dicom_tile_tests import axis_grid, tile_record, frame_from_decode
T_R, S_S = 4, 16
def rgb_img(Y, Co, Cg): return np.clip(np.stack(inv_rgb(Y, Co, Cg, fl), -1), 0, 255).astype(np.uint8)
def save(a, name, d): Image.fromarray(a).save(os.path.join(d, name)); print('wrote', name, a.shape, flush=True)
if __name__ == '__main__':
    path, csvp = sys.argv[1], sys.argv[2]; d = sys.argv[3] if len(sys.argv) > 3 else '../transcoding_review'; F = int(sys.argv[4]) if len(sys.argv) > 4 else 256
    g, rows = open_file(path, csvp); top, pad, Hp, Wp = g['top'], g['pad'], g['Hp'], g['Wp']; Hv, Wv = Hp - 2 * pad, Wp - 2 * pad
    Hs, LLs = {}, {}
    for c in range(3): Hs[c], LLs[c] = decode_coeffs(path, rows, g, c)
    # layer-3 planes: run the inverse only for scales top..3
    L3 = []
    for c in range(3):
        cur = LLs[c].astype(np.int32)
        for s in range(top, 2, -1):
            n_r, n_c = need(Hp, s), need(Wp, s); cur = cur[:n_r, :n_c]
            out = lift_axis(lift_axis(cur, Hs[c][s][0][:n_r, :n_c], 1), lift_axis(Hs[c][s][1][:n_r, :n_c], Hs[c][s][2][:n_r, :n_c], 1), 0)
            nxt = np.zeros_like(out); nxt[1:, 1:] = out[:-1, :-1]; cur = nxt
        L3.append(out)
    orig = rgb_img(*L3); save(orig, 'overview_original_L3.png', d)
    gy = axis_grid(0, Hv, pad, 5, F); gx = axis_grid(0, Wv, pad, 5, F)
    r = list(csv.DictReader(open(os.path.join(d, 'tile_stats_F%d_L0.csv' % F))))
    def col(k, t=float): return np.array([t(x[k]) for x in r])
    ty, tx = col('ty', int), col('tx', int); ex = col('exact_filler', int).astype(bool); rob = col('rgb_range_robust'); st = col('rgb_range'); outf = col('outside_frac')
    med = np.stack([col('med_r'), col('med_g'), col('med_b')], 1)
    masks = {'A': ex, 'R%d' % T_R: ex | ((rob <= T_R) & (outf == 0)), 'S%d' % S_S: ex | ((st <= S_S) & (outf == 0))}
    def blk(t, y, x):                                                   # layer-3 block of layer-0 tile (ty, tx): 8 level-0 px per layer-3 px
        p0y = gy['i_base'] + t[0] * F + 1; p0x = gx['i_base'] + t[1] * F + 1                 # padded level-0 position of the tile's first sample
        i3y = (p0y + 2) // 8 - 3; i3x = (p0x + 2) // 8 - 3; n = F // 8; return slice(max(i3y, 0), min(i3y + n, orig.shape[0])), slice(max(i3x, 0), min(i3x + n, orig.shape[1]))
    for nm, m in masks.items():
        ov = orig.astype(np.float32).copy(); tint = {'A': (255, 0, 255)}.get(nm, (255, 140, 0) if nm.startswith('R') else (0, 200, 255))
        for i in np.nonzero(m)[0]:
            sy, sx = blk((ty[i], tx[i]), 0, 0); ov[sy, sx] = 0.45 * ov[sy, sx] + 0.55 * np.array(tint)
        save(ov.astype(np.uint8), 'masks_%s.png' % nm, d); print('   %s: %d of %d tiles omitted (layer 0, F=%d)' % (nm, m.sum(), len(m), F))
    for nm in ('R%d' % T_R, 'S%d' % S_S):
        m = masks[nm]; glass = m & ~ex; fill = np.median(med[glass], axis=0).round().astype(np.uint8); out = orig.copy()
        for i in np.nonzero(m)[0]:
            sy, sx = blk((ty[i], tx[i]), 0, 0); out[sy, sx] = fill
        save(out, 'overview_%s_singlefill.png' % nm, d); save(np.clip(np.abs(out.astype(int) - orig.astype(int)) * 16, 0, 255).astype(np.uint8), 'diff_%s_x16.png' % nm, d)
        print('   %s single fill colour (median of omitted glass tiles): RGB %s' % (nm, tuple(int(v) for v in fill)))
    # full-resolution crops for borderline cases
    decs = [Windowed(Hs[c], LLs[c], top, Hp, Wp) for c in range(3)]; fillR = np.median(med[masks['R%d' % T_R] & ~ex], axis=0).round().astype(np.uint8)
    interior = (outf == 0)
    cases = {'speck_inside_robust_tile': np.nonzero(masks['R%d' % T_R] & ~ex & (st > 60) & interior)[0], 'glass_next_to_tissue': None, 'borderline_robust_5_to_8': np.nonzero((rob > T_R) & (rob <= 8) & interior)[0]}
    # a glass tile (omitted by R4) with the most non-omitted tiles among its 8 neighbours
    key = {(a, b): i for i, (a, b) in enumerate(zip(ty, tx))}; best = None; bn = -1
    for i in np.nonzero(masks['R%d' % T_R] & ~ex & interior)[0]:
        n = sum(1 for dy in (-1, 0, 1) for dx in (-1, 0, 1) if (dy or dx) and (ty[i] + dy, tx[i] + dx) in key and not masks['R%d' % T_R][key[(ty[i] + dy, tx[i] + dx)]])
        if n > bn: bn, best = n, i
    cases['glass_next_to_tissue'] = np.array([best])
    rng = np.random.default_rng(0)
    for nm, idx in cases.items():
        if len(idx) == 0: print('   no case for', nm); continue
        for j, i in enumerate(rng.permutation(idx)[:2] if nm != 'glass_next_to_tissue' else idx):
            cy, cx = int(ty[i]), int(tx[i]); y0 = gy['i_base'] + (cy - 1) * F; x0 = gx['i_base'] + (cx - 1) * F
            if y0 < 0 or x0 < 0: continue
            pix = [dd.out_region(0, y0, y0 + 3 * F, x0, x0 + 3 * F) for dd in decs]; a = rgb_img(*pix); o = a.copy()
            for dy in (-1, 0, 1):
                for dx in (-1, 0, 1):
                    k2 = key.get((cy + dy, cx + dx))
                    if k2 is not None and masks['R%d' % T_R][k2]: o[(dy + 1) * F:(dy + 2) * F, (dx + 1) * F:(dx + 2) * F] = fillR
            df = np.clip(np.abs(o.astype(int) - a.astype(int)) * 16, 0, 255).astype(np.uint8); sep = np.full((3 * F, 8, 3), 255, np.uint8)
            save(np.concatenate([a, sep, o, sep, df], 1), 'crops_%s_%d_tile_%d_%d.png' % (nm, j, cy, cx), d)
    print('done')
