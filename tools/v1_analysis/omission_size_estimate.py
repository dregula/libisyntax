# Per-tile statistics and real OpenJPEG sizes for every DICOM tile of every pyramid layer, to estimate what omitting tiles saves.
#   python3 -B -I omission_size_estimate.py <file.isyntax> <blocks.csv> [outdir] [sizes=256,512,1024]
# Route measured: wavelet-domain YCoCg->YBR_RCT (W-floor, see companion note 8d/8g), 8-bit, reversible 5/3 + RCT, D = 5, 64 x 64 code-blocks.
# For every tile: exact-filler flag (YCoCg == 255,0,0 everywhere), strict RGB range, robust RGB range (99.5 % - 0.5 % percentile), median RGB,
# share of pixels outside the decoded array (edge tiles, zero filler) and the encoded size.  One CSV per frame size and layer is written.
# Uses multiprocessing (fork); the decoded coefficient arrays are shared copy-on-write.  Read-only on the input.
import sys, os, io, time, csv, numpy as np
import multiprocessing as mp
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from PIL import Image
from colour_ebcot_probe import decode_coeffs, inv_rgb, rct_inv, convert_pyramid, fl
from windowed_idwt_halo import open_file, full_pyramid, inv_pyramid, Windowed
from dicom_tile_tests import axis_grid, tile_record, frame_from_decode
from openjpeg_ebcot_sizes import encode, to_u8
G = {}                                                  # shared state for the workers (set before the pool is created)
def work(args):
    k, F, ty, tx = args; D = 5; rec = tile_record(k, ty, tx, F, G['gy'][(k, F)], G['gx'][(k, F)], G['ox'], G['oy'])
    pix = []; filler = 0.0
    for c in range(3):
        fr, f = frame_from_decode(G['decs'][c], k, rec, F); pix.append(fr); filler = f
    Y, Co, Cg = pix
    exact = bool(Y.min() == 255 and Y.max() == 255 and Co.min() == 0 and Co.max() == 0 and Cg.min() == 0 and Cg.max() == 0)
    rgb = np.clip(np.stack(inv_rgb(Y, Co, Cg, fl), -1), 0, 255)
    strict = int((rgb.reshape(-1, 3).max(0) - rgb.reshape(-1, 3).min(0)).max())
    lo, hi = np.percentile(rgb.reshape(-1, 3), [0.5, 99.5], axis=0); robust = float((hi - lo).max()); med = np.median(rgb.reshape(-1, 3), axis=0)
    if exact: a, _ = to_u8([np.full((F, F), 255)] * 3)
    else:
        py = [full_pyramid(p, D) for p in pix]; (LLy, lvy), cbw, crw = convert_pyramid(py, fl)
        a, _ = to_u8(rct_inv(Y, inv_pyramid(*cbw), inv_pyramid(*crw)))
    return (k, F, ty, tx, rec['i0y'], rec['i0x'], int(exact), strict, round(robust, 2), round(float(med[0]), 1), round(float(med[1]), 1), round(float(med[2]), 1), round(filler, 4), len(encode(a, 64, D)))
if __name__ == '__main__':
    path, csvp = sys.argv[1], sys.argv[2]; outdir = sys.argv[3] if len(sys.argv) > 3 else '../transcoding_review'
    SIZES = [int(v) for v in sys.argv[4].split(',')] if len(sys.argv) > 4 else [256, 512, 1024]; os.makedirs(outdir, exist_ok=True); t0 = time.time()
    g, rows = open_file(path, csvp); top, pad = g['top'], g['pad']; Hs, LLs = {}, {}
    for c in range(3): Hs[c], LLs[c] = decode_coeffs(path, rows, g, c)
    Wv, Hv = g['Wp'] - 2 * pad, g['Hp'] - 2 * pad
    G.update(decs=[Windowed(Hs[c], LLs[c], top, g['Hp'], g['Wp']) for c in range(3)], ox=g['ox'], oy=g['oy'], gy={}, gx={})
    for k in range(top + 1):
        for F in SIZES: G['gy'][(k, F)] = axis_grid(k, Hv, pad, 5, F); G['gx'][(k, F)] = axis_grid(k, Wv, pad, 5, F)
    print('decoded (%.0fs)' % (time.time() - t0), flush=True)
    with mp.get_context('fork').Pool(16) as pool:      # fork: workers inherit the decoded arrays (Python 3.14 defaults to forkserver)
        for F in SIZES:
            for k in range(top + 1):
                gy, gx = G['gy'][(k, F)], G['gx'][(k, F)]; jobs = [(k, F, ty, tx) for ty in range(gy['tiles']) for tx in range(gx['tiles'])]
                res = pool.map(work, jobs, chunksize=max(1, len(jobs) // 64))
                with open(os.path.join(outdir, 'tile_stats_F%d_L%d.csv' % (F, k)), 'w', newline='') as f:
                    w = csv.writer(f); w.writerow(['layer', 'F', 'ty', 'tx', 'i0y', 'i0x', 'exact_filler', 'rgb_range', 'rgb_range_robust', 'med_r', 'med_g', 'med_b', 'outside_frac', 'bytes']); w.writerows(res)
                tot = sum(r[-1] for r in res); print('F=%4d layer %d: %5d tiles, %.3f MB (%.0fs)' % (F, k, len(res), tot / 1e6, time.time() - t0), flush=True)
    print('done %.0fs' % (time.time() - t0))
