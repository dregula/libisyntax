# Reference dump of DICOM frame coefficients (Python): wavelet-domain colour-converted YBR_RCT coefficients in OpenJPEG's Mallat layout, to compare byte for byte
# with the C++ strip transcoder (src/transcode/strip_transcoder.cpp --dump).
#   python3 -B -I dump_frame_coeffs.py <file.isyntax> <blocks.csv> <outdir> <floor|trunc> <layer> <ty> <tx> [<layer> <ty> <tx> ...]    (F = 256, D = 5, h = 4)
# Interior frames use build_windowed (ring recomputation from windowed strips), edge frames the local transform of the zero-padded decoded pixels; both as in dicom_tile_tests.py.
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from colour_ebcot_probe import decode_coeffs, convert_pyramid, fl, tr
from windowed_idwt_halo import open_file, full_pyramid, Windowed, build_windowed
from dicom_tile_tests import axis_grid, tile_record, frame_from_decode, is_interior
from make_coeff_test_frames import mallat
path, csvp, outdir, mode = sys.argv[1:5]; q = fl if mode == 'floor' else tr; D, F = 5, 256; rest = [int(v) for v in sys.argv[5:]]
g, rows = open_file(path, csvp); top, pad = g['top'], g['pad']; Hs, LLs = {}, {}
for c in range(3): Hs[c], LLs[c] = decode_coeffs(path, rows, g, c)
decs = [Windowed(Hs[c], LLs[c], top, g['Hp'], g['Wp']) for c in range(3)]; Wv, Hv = g['Wp'] - 2 * pad, g['Hp'] - 2 * pad; rng = np.random.default_rng(1)
for k, ty, tx in zip(rest[0::3], rest[1::3], rest[2::3]):
    gy, gx = axis_grid(k, Hv, pad, D, F), axis_grid(k, Wv, pad, D, F); rec = tile_record(k, ty, tx, F, gy, gx, g['ox'], g['oy']); py = []
    for c in range(3):
        if rec['i0y'] >= 8 and rec['i0x'] >= 8 and rec['i0y'] + F <= 2 * decs[c].lim[k][0] - 8 and rec['i0x'] + F <= 2 * decs[c].lim[k][1] - 8:
            LL, lv, _ = build_windowed(decs[c], k, rec['i0y'], rec['i0x'], F, D, 4, rng); py.append((LL, lv))
        else:
            fr, _ = frame_from_decode(decs[c], k, rec, F); py.append(full_pyramid(fr, D))
    (LLy, lvy), (LLb, lvb), (LLr, lvr) = convert_pyramid(py, q); LLy = LLy - 128
    A = np.stack([mallat(LLy, lvy), mallat(LLb, lvb), mallat(LLr, lvr)]).astype('<i4')
    A.tofile(os.path.join(outdir, 'py_L%d_%d_%d.coef' % (k, ty, tx))); print('layer', k, 'frame', ty, tx, 'origin', rec['i0y'], rec['i0x'], flush=True)
