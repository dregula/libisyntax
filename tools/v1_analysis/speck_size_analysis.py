# How big are the specks that the 'robust range' omission rule (8h, strategy B) would delete?  Read-only; numpy only.
#   python3 -B -I speck_size_analysis.py <file.isyntax> <blocks.csv> [stats_dir=../transcoding_review] [T=4] [F=256] [out.csv]
# Takes the layer-0 tiles that rule B omits (robust RGB range <= T, not exact filler, from tile_stats_F<F>_L0.csv), decodes each tile
# (exact pixel-domain RGB, floor inverse, as omission_size_estimate.py), and finds pixels that deviate from the tile's median colour (= the fill colour
# that would replace the tile) by more than DEV levels in any channel (noise on glass is about 4 levels, so DEV = 8 and 16 are real objects).
# 8-connected components of those pixels are the specks: area in pixels, equivalent circular diameter and longest bounding-box side in um
# (layer 0: 0.25 um per pixel, from the header's UFS_IMAGE_DIMENSION_SCALE_FACTOR; layer-0 pixel centres are 0.25 um apart).
import sys, os, csv, numpy as np, multiprocessing as mp
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from colour_ebcot_probe import decode_coeffs, inv_rgb, fl
from windowed_idwt_halo import open_file, Windowed
from dicom_tile_tests import axis_grid, tile_record, frame_from_decode
UM = 0.25
G = {}
def components(mask):
    ys, xs = np.nonzero(mask); pts = {(int(y), int(x)) for y, x in zip(ys, xs)}; out = []
    while pts:
        seed = pts.pop(); stack = [seed]; comp = [seed]
        while stack:
            y, x = stack.pop()
            for dy in (-1, 0, 1):
                for dx in (-1, 0, 1):
                    q = (y + dy, x + dx)
                    if q in pts: pts.remove(q); stack.append(q); comp.append(q)
        a = np.array(comp); out.append((len(comp), int(a[:, 0].max() - a[:, 0].min() + 1), int(a[:, 1].max() - a[:, 1].min() + 1)))
    return out
def work(args):
    ty, tx, F = args; rec = tile_record(0, ty, tx, F, G['gy'], G['gx'], G['ox'], G['oy']); pix = []
    for c in range(3): fr, _ = frame_from_decode(G['decs'][c], 0, rec, F); pix.append(fr)
    rgb = np.clip(np.stack(inv_rgb(*pix, fl), -1), 0, 255).astype(np.int16); med = np.median(rgb.reshape(-1, 3), axis=0)
    dev = np.abs(rgb - med).max(-1); res = {}
    for d in (8, 16): res[d] = components(dev > d)
    return ty, tx, int(dev.max()), res
if __name__ == '__main__':
    path, csvp = sys.argv[1], sys.argv[2]; sd = sys.argv[3] if len(sys.argv) > 3 else '../transcoding_review'
    T = float(sys.argv[4]) if len(sys.argv) > 4 else 4; F = int(sys.argv[5]) if len(sys.argv) > 5 else 256; outp = sys.argv[6] if len(sys.argv) > 6 else None
    rows_ = list(csv.DictReader(open(os.path.join(sd, 'tile_stats_F%d_L0.csv' % F))))
    sel = [(int(r['ty']), int(r['tx']), F) for r in rows_ if float(r['rgb_range_robust']) <= T and int(r['exact_filler']) == 0]
    print('layer 0, F=%d: %d tiles omitted by robust rule T=%g (of %d)' % (F, len(sel), T, len(rows_)), flush=True)
    g, rows = open_file(path, csvp); top, pad = g['top'], g['pad']; Hs, LLs = {}, {}
    for c in range(3): Hs[c], LLs[c] = decode_coeffs(path, rows, g, c)
    Wv, Hv = g['Wp'] - 2 * pad, g['Hp'] - 2 * pad
    G.update(decs=[Windowed(Hs[c], LLs[c], top, g['Hp'], g['Wp']) for c in range(3)], ox=g['ox'], oy=g['oy'], gy=axis_grid(0, Hv, pad, 5, F), gx=axis_grid(0, Wv, pad, 5, F))
    with mp.get_context('fork').Pool(16) as pool: res = pool.map(work, sel, chunksize=16)
    if outp:
        with open(outp, 'w', newline='') as f:
            w = csv.writer(f); w.writerow(['ty', 'tx', 'dev', 'area_px', 'bbox_h', 'bbox_w'])
            for ty, tx, mx, r in res:
                for a, h, ww in r[8]: w.writerow([ty, tx, mx, a, h, ww])
    for d in (8, 16):
        comps = [c for _, _, _, r in res for c in r[d]]; A = np.array([c[0] for c in comps]) if comps else np.zeros(0); L = np.array([max(c[1], c[2]) for c in comps]) if comps else np.zeros(0)
        tiles_with = sum(1 for _, _, _, r in res if r[d])
        print('\n== deviation from the fill colour > %d levels (any channel): %d specks in %d of %d omitted tiles (%.1f %%)' % (d, len(comps), tiles_with, len(res), 100 * tiles_with / len(res)))
        if not len(comps): continue
        eq = 2 * np.sqrt(A / np.pi) * UM; ext = L * UM
        pc = lambda v, q: np.percentile(v, q)
        print('  area (px):                 median %d  p90 %d  p99 %d  max %d' % (pc(A, 50), pc(A, 90), pc(A, 99), A.max()))
        print('  equivalent diameter (um):  median %.2f  p90 %.2f  p99 %.2f  max %.2f' % (pc(eq, 50), pc(eq, 90), pc(eq, 99), eq.max()))
        print('  longest bbox side (um):    median %.2f  p90 %.2f  p99 %.2f  max %.2f' % (pc(ext, 50), pc(ext, 90), pc(ext, 99), ext.max()))
        for lim in (2, 5, 10, 20):
            print('  specks with equivalent diameter < %2d um: %5.1f %% ; longest side < %2d um: %5.1f %%' % (lim, 100 * np.mean(eq < lim), lim, 100 * np.mean(ext < lim)))
        big = sorted(comps, reverse=True)[:5]; print('  five largest (area px, bbox h x w px):', ', '.join('%d (%dx%d)' % c for c in big))
    mx = np.array([m for _, _, m, _ in res]); print('\nlargest single-pixel deviation per tile: median %d, p95 %d, max %d' % (np.median(mx), np.percentile(mx, 95), mx.max()))
