# Test data for the OpenJPEG coefficient-input patch (ext/openjpeg-2.5.4, patches/openjpeg-2.5.4-coefficient-input.patch). Read-only on the input.
#   python3 -B -I make_coeff_test_frames.py <file.isyntax> <blocks.csv> <outdir>
# For three real 256 x 256 layer-0 frames writes   <name>_rgb.raw   (256*256*3 uint8, interleaved R G B; truncating inverse, clamped)
#                                                  <name>_coef.raw  (3 components x 256*256 int32, little endian): the YBR_RCT coefficients that
# OpenJPEG's own encoder would hold after DC level shift, RCT and the 5-level reversible 5/3 DWT, in OpenJPEG's Mallat layout
# (LL5 top-left; per resolution HL to the right, LH below, HH diagonal; row stride 256), computed here with our NumPy transforms.
# The C test (tools/transcode_tests/coeff_input_test.c) encodes the RGB with stock OpenJPEG and the coefficients through the patched
# coefficient-input path; the two codestreams must be byte-identical.
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from colour_ebcot_probe import decode_coeffs, inv_rgb, tr
from windowed_idwt_halo import open_file, full_pyramid, Windowed
D = 5; F = 256
def mallat(LL, lv):
    A = np.zeros((F, F), dtype=np.int32); h = F >> D; A[:h, :h] = LL
    for j in range(D - 1, -1, -1):                 # lv[j] = (HL, LH, HH) of the transform of the (F >> j) image; coarsest first when going down
        HL, LH, HH = lv[j]; n = F >> (j + 1)       # sub-band size
        A[:n, n:2 * n] = HL; A[n:2 * n, :n] = LH; A[n:2 * n, n:2 * n] = HH
    return A
if __name__ == '__main__':
    path, csvp, outdir = sys.argv[1:4]; os.makedirs(outdir, exist_ok=True)
    g, rows = open_file(path, csvp); top, pad = g['top'], g['pad']; Hs, LLs = {}, {}
    for c in range(3): Hs[c], LLs[c] = decode_coeffs(path, rows, g, c)
    decs = [Windowed(Hs[c], LLs[c], top, g['Hp'], g['Wp']) for c in range(3)]
    for name, i0y, i0x in (('tissue', 64 * 157 + pad, 64 * 237 + pad), ('edge', 64 * 137 + pad, 64 * 137 + pad), ('glass', 64 * 37 + pad, 64 * 37 + pad)):
        pix = [d.out_region(0, i0y, i0y + F, i0x, i0x + F) for d in decs]
        R, G, B = [np.clip(p, 0, 255).astype(np.int32) for p in inv_rgb(pix[0], pix[1], pix[2], tr)]
        rgb = np.stack([R, G, B], -1).astype(np.uint8); rgb.tofile(os.path.join(outdir, name + '_rgb.raw'))
        Rs, Gs, Bs = R - 128, G - 128, B - 128                       # DC level shift first (T.800 G.1), then the RCT
        Y = (Rs + 2 * Gs + Bs) >> 2; Cb = Bs - Gs; Cr = Rs - Gs
        comps = [mallat(*full_pyramid(p, D)) for p in (Y, Cb, Cr)]
        np.stack(comps).astype('<i4').tofile(os.path.join(outdir, name + '_coef.raw')); print(name, 'written', flush=True)
