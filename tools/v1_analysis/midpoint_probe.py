# What if the structurally-zero least-significant planes are not coded (J2K coding passes dropped) and the decoder reconstructs
# with reconstruction parameter r = 1/2 (OpenJPEG) instead of r = 0?  T.800 E.1.2.2 (E-8): R = q +- r * 2^(Mb-Nb) for Nb < Mb.
# Read-only; needs numpy.  python3 -B -I midpoint_probe.py <file.isyntax> <blocks.csv>
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from colour_ebcot_probe import *
path, csvp = sys.argv[1], sys.argv[2]
g, rows = open_file(path, csvp); top, pad = g['top'], g['pad']; Hs, LLs = {}, {}
for c in range(3): Hs[c], LLs[c] = decode_coeffs(path, rows, g, c)
decs = [Windowed(Hs[c], LLs[c], top, g['Hp'], g['Wp']) for c in range(3)]
D, CB = 5, 64
Q = {'Y': [1, 0, 0, 0, 0], 'Cb': [3, 1, 0, 0, 0], 'Cr': [3, 1, 0, 0, 0]}          # empty LSB planes per J2K level (measured in colour_ebcot_probe.py)
def midpoint(lv, j, F, q):
    """add +-2^(q-1) to non-zero coefficients of the code-blocks that do not touch the frame border (ring-free blocks)."""
    if q == 0: return lv[j], 0
    out = []; touched = 0
    for a in lv[j]:
        a = a.copy(); n = a.shape[0]; m = n // CB
        if m >= 3:
            sl = slice(CB, n - CB); sub = a[sl, sl]; nz = sub != 0; touched += int(nz.sum())
            sub[nz] += np.sign(sub[nz]) * (1 << (q - 1))
        out.append(a)
    return tuple(out), touched
for name, i0y, i0x, F in (('tissue', 64 * 157 + pad, 64 * 237 + pad, 1024), ('tissue', 64 * 157 + pad, 64 * 237 + pad, 512), ('tissue', 64 * 157 + pad, 64 * 237 + pad, 256)):
    pix = [d.out_region(0, i0y, i0y + F, i0x, i0x + F) for d in decs]; py = [full_pyramid(p, D) for p in pix]
    (LLy, lvy), (LLb, lvb), (LLr, lvr) = convert_pyramid(py, fl)
    base = [inv_pyramid(LLy, lvy), inv_pyramid(LLb, lvb), inv_pyramid(LLr, lvr)]
    mod = {}; info = []
    for nm, (LL, lv) in zip(('Y', 'Cb', 'Cr'), ((LLy, lvy), (LLb, lvb), (LLr, lvr))):
        lv2 = list(lv); cnt = 0
        for j in range(D):
            lv2[j], t = midpoint(lv, j, F, Q[nm][j]); cnt += t
        mod[nm] = inv_pyramid(LL, lv2); info.append('%s: %d coefficients shifted' % (nm, cnt))
    RGBb = rct_inv(*base); RGBm = rct_inv(mod['Y'], mod['Cb'], mod['Cr'])
    dRGB = np.stack([a - b for a, b in zip(RGBm, RGBb)])
    print('%s %dx%d (code-blocks %dx%d, ring-free ones only) | %s' % (name, F, F, CB, CB, ', '.join(info)))
    print('   r=1/2 vs r=0 pixels: Y differs %.2f%% (max %d) | Cb differs %.2f%% (max %d) | Cr differs %.2f%% (max %d) | RGB differs %.2f%% (max |d| %d)' % (
        100 * np.mean(mod['Y'] != base[0]), np.abs(mod['Y'] - base[0]).max(), 100 * np.mean(mod['Cb'] != base[1]), np.abs(mod['Cb'] - base[1]).max(),
        100 * np.mean(mod['Cr'] != base[2]), np.abs(mod['Cr'] - base[2]).max(), 100 * np.mean((dRGB != 0).any(0)), np.abs(dRGB).max()))
# ring-free share of code-blocks, by frame size and code-block size, for J2K levels 1 and 2 (the only ones with a quantizer structure)
print('\nshare of code-blocks that touch no frame border (ring-free), per level (n = F / 2^level coefficients per side):')
for cb in (64, 32):
    for F in (256, 512, 1024):
        s = []
        for j in (1, 2):
            n = F >> j; m = n // cb; s.append('level %d: %d blocks/side -> %s' % (j, m, ('%.0f%%' % (100 * (m - 2) ** 2 / m ** 2)) if m >= 2 else 'single block, 0%'))
        print('  code-block %2d, frame %4d: %s' % (cb, F, ' | '.join(s)))
