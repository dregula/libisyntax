# Single-core timing of OpenJPEG (stock Tier-1 is identical to the patched coefficient-input encoder) on real 256-px frames of every layer.
#   python3 -B -I tier1_profile.py <file.isyntax> <blocks.csv> <stats_dir with tile_stats_F256_L*.csv> [out.csv] [tiles_per_layer=24]
# Frames are the W-floor route of omission_size_estimate.py (wavelet-domain colour conversion, decoded to 8-bit pixels so that the encoder sees exactly
# the converted coefficients). Stock encode = forward RCT + 5/3 DWT + Tier-1 + Tier-2 + headers. The patched exact route skips the first two.
# Their cost is estimated from encoding a CONSTANT frame (identical DWT/MCT work, trivial Tier-1), so  Tier-1+2 ~ t(frame) - t(constant).
# Timing: time.perf_counter around the encode call, best of N repeats, one process, nothing else of ours running in it.
# Extrapolation: per-layer least-squares fit  t = a + b * coded_bytes  on the sampled frames, applied to the bytes of EVERY frame of that layer
# (the per-frame bytes come from the existing tile_stats CSVs, so no further decoding is needed for the whole-slide estimate).
import sys, os, io, time, csv, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from PIL import Image
from colour_ebcot_probe import decode_coeffs, inv_rgb, rct_inv, convert_pyramid, fl
from windowed_idwt_halo import open_file, full_pyramid, inv_pyramid, Windowed
from dicom_tile_tests import axis_grid, tile_record, frame_from_decode
from openjpeg_ebcot_sizes import encode, to_u8
path, csvp, sdir = sys.argv[1], sys.argv[2], sys.argv[3]
out = sys.argv[4] if len(sys.argv) > 4 else 'tier1_profile.csv'; NS = int(sys.argv[5]) if len(sys.argv) > 5 else 24
F, D, REP = 256, 5, 5
t0 = time.time(); g, rows = open_file(path, csvp); top, pad = g['top'], g['pad']; Hs, LLs = {}, {}
for c in range(3): Hs[c], LLs[c] = decode_coeffs(path, rows, g, c)
decs = [Windowed(Hs[c], LLs[c], top, g['Hp'], g['Wp']) for c in range(3)]
Wv, Hv = g['Wp'] - 2 * pad, g['Hp'] - 2 * pad; print('decoded (%.0fs)' % (time.time() - t0), flush=True)
def best(a, n=REP):
    ts = []
    for _ in range(n):
        t = time.perf_counter(); d = encode(a, 64, D); ts.append(time.perf_counter() - t)
    return min(ts), len(d)
const = {v: best(np.full((F, F, 3), v, np.uint8), 9)[0] for v in (255, 230, 128)}
tconst = float(np.median(list(const.values()))); print('constant-frame encode (DWT+MCT+headers, trivial Tier-1): ', {k: round(v * 1e3, 2) for k, v in const.items()}, 'ms', flush=True)
res = []; fits = {}
for k in range(top + 1):
    gy, gx = axis_grid(k, Hv, pad, 5, F), axis_grid(k, Wv, pad, 5, F)
    allrows = list(csv.DictReader(open(os.path.join(sdir, 'tile_stats_F256_L%d.csv' % k))))
    nz = sorted([r for r in allrows if r['exact_filler'] == '0'], key=lambda r: int(r['bytes']))
    idx = sorted(set(np.linspace(0, len(nz) - 1, min(NS, len(nz))).astype(int).tolist())) if nz else []
    pick = [nz[i] for i in idx] + [r for r in allrows if r['exact_filler'] == '1'][:2]
    pts = []
    for r in pick:
        ty, tx = int(r['ty']), int(r['tx']); rec = tile_record(k, ty, tx, F, gy, gx, g['ox'], g['oy']); pix = []
        for c in range(3): fr, _ = frame_from_decode(decs[c], k, rec, F); pix.append(fr)
        if r['exact_filler'] == '1': a, _ = to_u8([np.full((F, F), 255)] * 3)
        else:
            py = [full_pyramid(p, D) for p in pix]; (LLy, lvy), cbw, crw = convert_pyramid(py, fl)
            a, _ = to_u8(rct_inv(pix[0], inv_pyramid(*cbw), inv_pyramid(*crw)))
        t, nb = best(a); pts.append((nb, t)); res.append((k, ty, tx, int(r['exact_filler']), nb, t))
    P = np.array([p for p in pts if p[0] > 400] or pts)
    b, a_ = np.polyfit(P[:, 0], P[:, 1], 1) if len(P) > 2 else (0.0, float(P[:, 1].mean()))
    fits[k] = (a_, b); allb = np.array([int(r['bytes']) for r in allrows], dtype=float)
    est = float((a_ + b * allb).sum()); est_t12 = float((np.maximum(a_ + b * allb - tconst, 0)).sum())
    print('layer %d: %6d frames, sampled %2d | frame bytes median %7.0f | encode ms: min %.2f median %.2f max %.2f | fit t = %.2f ms + %.3f us/byte (%.1f MB/s coded) | est. whole layer %.0f core-s (Tier-1+2 only %.0f)'
          % (k, len(allrows), len(pick), np.median(allb), P[:, 1].min() * 1e3, np.median(P[:, 1]) * 1e3, P[:, 1].max() * 1e3, a_ * 1e3, b * 1e6, 1 / b / 1e6 if b > 0 else 0, est, est_t12), flush=True)
with open(out, 'w', newline='') as f:
    w = csv.writer(f); w.writerow(['layer', 'ty', 'tx', 'filler', 'bytes', 'encode_s']); w.writerows(res)
tot = 0.0; tot12 = 0.0; totb = 0.0; nfr = 0
for k in range(top + 1):
    allrows = list(csv.DictReader(open(os.path.join(sdir, 'tile_stats_F256_L%d.csv' % k)))); allb = np.array([int(r['bytes']) for r in allrows], dtype=float); a_, b = fits[k]
    tot += float((a_ + b * allb).sum()); tot12 += float(np.maximum(a_ + b * allb - tconst, 0).sum()); totb += allb.sum(); nfr += len(allb)
print('\nWHOLE PYRAMID (this file, 256-px frames, %d frames, %.1f MB coded): stock encode %.0f core-s; Tier-1+2 only (constant-frame overhead removed) %.0f core-s -> %.1f MB/s coded per core, %.2f ms per frame' % (nfr, totb / 1e6, tot, tot12, totb / 1e6 / tot12, tot12 / nfr * 1e3))
print('overhead share avoided by the exact route (DWT+MCT per frame): %.0f%% of the stock encode time' % (100 * (1 - tot12 / tot)))
print('total %.0fs' % (time.time() - t0))
