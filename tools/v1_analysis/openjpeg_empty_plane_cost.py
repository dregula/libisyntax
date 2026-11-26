# What do the always-zero least-significant planes cost in a real EBCOT stream?  Coefficients are coded MSB plane first, so the empty LSB
# planes are coded last; removing them (dividing the stored multiples of 2^Q by 2^Q) leaves the earlier planes' coding identical.  The size
# difference is therefore the cost of coding those planes (upper bound on what plane deletion with r = 0 could save).
# Ring (halo) coefficients are arbitrary integers and are left unscaled, i.e. the estimate is slightly conservative.  Needs numpy + Pillow.
#   python3 -B -I openjpeg_empty_plane_cost.py <file.isyntax> <blocks.csv>
import sys, os, io, time, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from PIL import Image
from colour_ebcot_probe import decode_coeffs, rct_inv, convert_pyramid, fl
from windowed_idwt_halo import open_file, full_pyramid, inv_pyramid, Windowed
from openjpeg_ebcot_sizes import encode, to_u8
def scale(pyr, Q):
    LL, lv = pyr; out = []
    for j, bands in enumerate(lv):
        q = Q[j]; nb = []
        for a in bands:
            if q: m = (a & ((1 << q) - 1)) == 0; a = np.where(m, a >> q, a)
            nb.append(a)
        out.append(tuple(nb))
    return LL, out
if __name__ == '__main__':
    path, csvp = sys.argv[1], sys.argv[2]; D = 5; t0 = time.time()
    g, rows = open_file(path, csvp); top, pad = g['top'], g['pad']; Hs, LLs = {}, {}
    for c in range(3): Hs[c], LLs[c] = decode_coeffs(path, rows, g, c)
    decs = [Windowed(Hs[c], LLs[c], top, g['Hp'], g['Wp']) for c in range(3)]
    QY = [1, 0, 0, 0, 0]; QC = [3, 1, 0, 0, 0]; Z = [0] * 5
    for rname, y0, x0 in (('tissue', 64 * 157 + pad, 64 * 237 + pad), ('blank glass', 64 * 37 + pad, 64 * 37 + pad)):
        R = 2048
        for F in (256, 512, 1024):
            tot = {'baseline': 0, 'no empty Y planes': 0, 'no empty chroma planes': 0, 'neither': 0}; clip = 0
            for ty in range(R // F):
                for tx in range(R // F):
                    iy, ix = y0 + ty * F, x0 + tx * F
                    pix = [d.out_region(0, iy, iy + F, ix, ix + F) for d in decs]; py = [full_pyramid(p, D) for p in pix]
                    ypy, cbw, crw = convert_pyramid(py, fl)
                    for name, qy, qc in (('baseline', Z, Z), ('no empty Y planes', QY, Z), ('no empty chroma planes', Z, QC), ('neither', QY, QC)):
                        Y = inv_pyramid(*scale(ypy, qy)); Cb = inv_pyramid(*scale(cbw, qc)); Cr = inv_pyramid(*scale(crw, qc))
                        a, bad = to_u8(rct_inv(Y, Cb, Cr)); clip += bad; tot[name] += len(encode(a, 64, D))
            b = tot['baseline']
            print('%-11s frames %4d: baseline %.3f MB | no empty Y planes %.3f MB (%.2f%% smaller) | no empty chroma planes %.3f MB (%.2f%%) | neither %.3f MB (%.2f%%)' % (
                rname, F, b / 1e6, tot['no empty Y planes'] / 1e6, 100 * (1 - tot['no empty Y planes'] / b), tot['no empty chroma planes'] / 1e6,
                100 * (1 - tot['no empty chroma planes'] / b), tot['neither'] / 1e6, 100 * (1 - tot['neither'] / b)), flush=True)
    print('total %.0fs' % (time.time() - t0))
