# Real EBCOT (OpenJPEG via Pillow) sizes for the DICOM transcoding routes on real regions of a v1 iSyntax file.
#   python3 -B -I openjpeg_ebcot_sizes.py <file.isyntax> <blocks.csv> [region=2048]      (needs numpy + Pillow with OpenJPEG; ~10 GB RAM)
# Each frame is encoded as an independent raw codestream: reversible 5/3, RCT (mct=1), 5 decomposition levels, one JPEG 2000 tile.
# Routes (all carry the SAME iSyntax decode, differing only in how YCoCg becomes YBR_RCT):
#   E-floor / E-fix : exact pixel-domain conversion; RGB = Philips inverse with floor (>>1) or fix (truncation), then OpenJPEG's RCT.
#   W-floor         : wavelet-domain conversion Cb = -(Co>>1) - Cg, Cr = Cb + Co, Y' = Y on each frame's own pyramid, decoded back
#                     to pixels (so the encoder sees exactly the sparse converted coefficients; RCT is a bijection on integers).
# Pixels outside 0..255 cannot be fed to an 8-bit encoder; they are clipped and counted (this changes the image slightly there).
# The iSyntax size of the same region is the area-weighted compressed size of the codeblocks overlapping it.
import sys, os, io, time, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from PIL import Image
from colour_ebcot_probe import (decode_coeffs, inv_rgb, rct_inv, convert_pyramid, fl, tr)
from windowed_idwt_halo import open_file, full_pyramid, inv_pyramid, Windowed

def encode(rgb, cb=64, D=5):
    b = io.BytesIO(); Image.fromarray(rgb).save(b, format='JPEG2000', no_jp2=True, irreversible=False, mct=1, num_resolutions=D + 1, codeblock_size=(cb, cb))
    return b.getvalue()
def to_u8(planes):
    a = np.stack(planes, -1); bad = int(((a < 0) | (a > 255)).any(-1).sum()); return np.clip(a, 0, 255).astype(np.uint8), bad
def isyntax_region_bytes(rows, g, x0p, y0p, w, h):
    """area-weighted compressed bytes of all codeblocks (all scales, colours) overlapping the padded level-0 rectangle."""
    tot = 0.0
    for r in rows:
        span = 256 << r['s']                                                         # a codeblock represents about span x span level-0 px
        x0 = r['x'] - g['ox']; y0 = r['y'] - g['oy']                                 # its first coefficient, padded level-0 units (error < 2^s px)
        ox_ = max(0, min(x0 + span, x0p + w) - max(x0, x0p)); oy_ = max(0, min(y0 + span, y0p + h) - max(y0, y0p))
        tot += r['z'] * (ox_ * oy_) / (span * span)
    return tot
if __name__ == '__main__':
    path, csvp = sys.argv[1], sys.argv[2]; R = int(sys.argv[3]) if len(sys.argv) > 3 else 2048; D = 5; t0 = time.time()
    g, rows = open_file(path, csvp); top, pad = g['top'], g['pad']; Hs, LLs = {}, {}
    for c in range(3): Hs[c], LLs[c] = decode_coeffs(path, rows, g, c)
    decs = [Windowed(Hs[c], LLs[c], top, g['Hp'], g['Wp']) for c in range(3)]; print('decoded (%.0fs)' % (time.time() - t0), flush=True)
    regions = [('tissue', 64 * 157 + pad, 64 * 237 + pad), ('tissue edge', 64 * 137 + pad, 64 * 137 + pad), ('blank glass', 64 * 37 + pad, 64 * 37 + pad)]
    for rname, y0, x0 in regions:
        isx = isyntax_region_bytes(rows, g, x0 + 1, y0 + 1, R, R)                   # out index i <-> padded position i + 1 at scale 0
        print('\n== region %-12s %dx%d px at level 0: iSyntax (area-weighted) %.3f MB = %.3f bits/px' % (rname, R, R, isx / 1e6, isx * 8 / R / R), flush=True)
        for F in (256, 512, 1024, 2048):
            if F > R: continue
            res = {'E-floor': 0, 'E-fix': 0, 'W-floor': 0}; clip = {'E-floor': 0, 'E-fix': 0, 'W-floor': 0}; res32 = {'W-floor': 0}; ok = True; te = 0.0
            for ty in range(R // F):
                for tx in range(R // F):
                    iy, ix = y0 + ty * F, x0 + tx * F
                    pix = [d.out_region(0, iy, iy + F, ix, ix + F) for d in decs]
                    for nm, h in (('E-floor', fl), ('E-fix', tr)):
                        a, bad = to_u8(inv_rgb(pix[0], pix[1], pix[2], h)); clip[nm] += bad; t1 = time.time(); d = encode(a, 64, D); te += time.time() - t1; res[nm] += len(d)
                        if tx == 0 and ty == 0: ok &= bool((np.array(Image.open(io.BytesIO(d))) == a).all())
                    py = [full_pyramid(p, D) for p in pix]; (LLy, lvy), cbw, crw = convert_pyramid(py, fl)
                    a, bad = to_u8(rct_inv(pix[0], inv_pyramid(*cbw), inv_pyramid(*crw))); clip['W-floor'] += bad
                    t1 = time.time(); d = encode(a, 64, D); te += time.time() - t1; res['W-floor'] += len(d)
                    if tx == 0 and ty == 0: ok &= bool((np.array(Image.open(io.BytesIO(d))) == a).all())
                    res32['W-floor'] += len(encode(a, 32, D))
            n = R * R
            print(' frame %4d x %4d (%3d frames): ' % (F, F, (R // F) ** 2) + ' | '.join('%s %.3f MB (%.3f b/px, %.2fx iSyntax, %d clipped)' % (k, v / 1e6, v * 8 / n, v / isx, clip[k]) for k, v in res.items())
                  + ' | W-floor with 32x32 code-blocks %.3f MB | round trips lossless: %s | encode %.0fs' % (res32['W-floor'] / 1e6, ok, te), flush=True)
    print('\ntotal %.0fs' % (time.time() - t0))
