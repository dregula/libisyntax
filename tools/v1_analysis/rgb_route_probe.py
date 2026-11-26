# YCoCg-R -> DICOM Photometric Interpretation RGB (no component transform) on a v1 iSyntax file. Read-only; numpy + Pillow(OpenJPEG), ~10 GB RAM.
#   python3 -B -I rgb_route_probe.py <file.isyntax> <blocks.csv> [A] [B] [C]      (parts; default all)
# Counterpart of colour_ebcot_probe.py / colour_sparsity_cost.py / openjpeg_ebcot_sizes.py, which target YBR_RCT.
# With RGB the codestream has SGcod = 0 (mct=0): R, G, B are coded as three independent components, each with its own 5/3 pyramid.
#   Part A  wavelet-domain conversion: the YCoCg-R inverse applied coefficient by coefficient (t = Y - h(Cg); G = t + Cg;
#           B = t - h(Co); R = B + Co) on LL and every detail band of the frame's own pyramid, then decoded, versus the exact
#           pixel-domain RGB (the libisyntax output). Error of R, G, B per number of decomposition levels D.
#   Part B  sparsity of the detail coefficients of R, G, B per level: wavelet-domain (stays sparse?) versus exact pixel-domain.
#   Part C  real EBCOT sizes (OpenJPEG): RGB with mct=0 for routes E-floor / E-fix / W-floor, and for comparison the same
#           frames as YBR_RCT (mct=1 on the exact RGB; and the earlier wavelet-domain RCT route), versus iSyntax bytes.
import sys, os, io, time, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from PIL import Image
from colour_ebcot_probe import decode_coeffs, inv_rgb, rct_inv, convert_pyramid, fl, tr, rct_fwd
from windowed_idwt_halo import open_file, full_pyramid, inv_pyramid, Windowed
from openjpeg_ebcot_sizes import isyntax_region_bytes, to_u8

def convert_rgb_pyramid(py, h):
    """py = [(LL, levels)] for Y, Co, Cg  ->  three pyramids (R, G, B), coefficient-wise YCoCg-R inverse with halving h."""
    (LLy, lvy), (LLo, lvo), (LLg, lvg) = py
    f = lambda Y, Co, Cg: inv_rgb(Y, Co, Cg, h)               # returns (R, G, B)
    LL = f(LLy, LLo, LLg); lv = [[], [], []]
    for ty, to, tg in zip(lvy, lvo, lvg):
        bands = [f(ty[b], to[b], tg[b]) for b in range(3)]
        for ch in range(3): lv[ch].append(tuple(bands[b][ch] for b in range(3)))
    return [(LL[ch], lv[ch]) for ch in range(3)]

def encode(rgb, mct, cb=64, D=5):
    b = io.BytesIO(); Image.fromarray(rgb).save(b, format='JPEG2000', no_jp2=True, irreversible=False, mct=mct, num_resolutions=D + 1, codeblock_size=(cb, cb))
    return b.getvalue()

def h0(a):
    v, cnt = np.unique(a.ravel(), return_counts=True); p = cnt / a.size; return float(-(p * np.log2(p)).sum())

def part_a(decs, frames):
    print('\n== PART A: wavelet-domain YCoCg-R -> RGB versus exact pixel-domain RGB (per-frame mirror pyramids)')
    for name, i0y, i0x, F in frames:
        pix = [d.out_region(0, i0y, i0y + F, i0x, i0x + F) for d in decs]
        for hn, h in (('floor', fl), ('trunc', tr)):
            ref = inv_rgb(pix[0], pix[1], pix[2], h); print(' frame %-12s %dx%d, halving = %s (reference: same halving in the pixel domain)' % (name, F, F, hn))
            for D in (1, 2, 3, 4, 5):
                py = [full_pyramid(p, D) for p in pix]; pr = convert_rgb_pyramid(py, h)
                out = [inv_pyramid(*p) for p in pr]; e = [o - r for o, r in zip(out, ref)]
                anyd = np.stack([x != 0 for x in e]).any(0); mx = max(int(np.abs(x).max()) for x in e)
                print('   D=%d  mean err R %+.2f G %+.2f B %+.2f | std %.2f %.2f %.2f | exact R %5.1f%% G %5.1f%% B %5.1f%% | any channel differs %5.1f%% (max |d| %d)' % (
                    D, e[0].mean(), e[1].mean(), e[2].mean(), e[0].std(), e[1].std(), e[2].std(),
                    *[100 * np.mean(x == 0) for x in e], 100 * anyd.mean(), mx), flush=True)

def lowbit(a):
    o = int(np.bitwise_or.reduce(np.abs(a).ravel())); return (o & -o).bit_length() - 1 if o else '-'

def part_b(decs, frames, D=5):
    print('\n== PART B: detail-coefficient sparsity of R, G, B (non-zero %% / zero-order entropy bits per coefficient / empty LSB planes)')
    for name, i0y, i0x, F in frames:
        pix = [d.out_region(0, i0y, i0y + F, i0x, i0x + F) for d in decs]; py = [full_pyramid(p, D) for p in pix]
        ex = [full_pyramid(c, D) for c in inv_rgb(pix[0], pix[1], pix[2], fl)]; wd = convert_rgb_pyramid(py, fl)
        print('\n%s %dx%d, halving floor. level 1 = finest. Stored Y, Co, Cg for reference.' % (name, F, F))
        print(' level | Y stored         Co stored        Cg stored        | G wavelet-dom    G exact          | R wavelet-dom    R exact          | B wavelet-dom    B exact')
        def st(lv, j): a = np.concatenate([x.ravel() for x in lv[j]]); return '%6.2f%% %5.2fb L%s' % (100 * np.mean(a != 0), h0(a), lowbit(a))
        for j in range(D):
            print(' %5d | %-16s %-16s %-16s | %-16s %-16s | %-16s %-16s | %-16s %-16s' % (j + 1, st(py[0][1], j), st(py[1][1], j), st(py[2][1], j),
                  st(wd[1][1], j), st(ex[1][1], j), st(wd[0][1], j), st(ex[0][1], j), st(wd[2][1], j), st(ex[2][1], j)), flush=True)

def part_c(decs, rows, g, R=2048, D=5):
    print('\n== PART C: OpenJPEG sizes, %dx%d regions, independent frames, reversible 5/3, %d levels, 64x64 code-blocks' % (R, R, D))
    pad = g['pad']
    regions = [('tissue', 64 * 157 + pad, 64 * 237 + pad), ('tissue edge', 64 * 137 + pad, 64 * 137 + pad), ('blank glass', 64 * 37 + pad, 64 * 37 + pad)]
    for rname, y0, x0 in regions:
        isx = isyntax_region_bytes(rows, g, x0 + 1, y0 + 1, R, R)
        print('\n region %-12s iSyntax %.3f MB = %.3f bits/px' % (rname, isx / 1e6, isx * 8 / R / R), flush=True)
        for F in (256, 512, 1024):
            keys = ['RGB E-floor', 'RGB E-fix', 'RGB W-floor', 'RCT E-floor', 'RCT W-floor']; res = {k: 0 for k in keys}; clip = {k: 0 for k in keys}; ok = True
            for ty in range(R // F):
                for tx in range(R // F):
                    iy, ix = y0 + ty * F, x0 + tx * F; first = tx == 0 and ty == 0
                    pix = [d.out_region(0, iy, iy + F, ix, ix + F) for d in decs]; arrs = {}
                    for nm, h in (('E-floor', fl), ('E-fix', tr)): arrs['RGB ' + nm], clip['RGB ' + nm] = to_u8(inv_rgb(pix[0], pix[1], pix[2], h))
                    py = [full_pyramid(p, D) for p in pix]
                    arrs['RGB W-floor'], clip['RGB W-floor'] = to_u8([inv_pyramid(*p) for p in convert_rgb_pyramid(py, fl)])
                    arrs['RCT E-floor'] = arrs['RGB E-floor']
                    (LLy, lvy), cbw, crw = convert_pyramid(py, fl)
                    arrs['RCT W-floor'], clip['RCT W-floor'] = to_u8(rct_inv(pix[0], inv_pyramid(*cbw), inv_pyramid(*crw)))
                    for k in keys:
                        d = encode(arrs[k], 0 if k.startswith('RGB') else 1, 64, D); res[k] += len(d)
                        if first: ok &= bool((np.array(Image.open(io.BytesIO(d))) == arrs[k]).all())
            n = R * R
            print('  frame %4d: ' % F + ' | '.join('%s %.2fx (%.3f b/px%s)' % (k, v / isx, v * 8 / n, ', %d clipped' % clip[k] if clip.get(k) else '') for k, v in res.items()) + ' | round trip lossless: %s' % ok, flush=True)

if __name__ == '__main__':
    path, csvp = sys.argv[1], sys.argv[2]; parts = set(sys.argv[3:]) or {'A', 'B', 'C'}; t0 = time.time()
    g, rows = open_file(path, csvp); top, pad = g['top'], g['pad']; Hs, LLs = {}, {}
    for c in range(3): Hs[c], LLs[c] = decode_coeffs(path, rows, g, c)
    decs = [Windowed(Hs[c], LLs[c], top, g['Hp'], g['Wp']) for c in range(3)]; print('decoded (%.0fs)' % (time.time() - t0), flush=True)
    if 'A' in parts: part_a(decs, [('tissue', 64 * 157 + pad, 64 * 237 + pad, 256), ('tissue edge', 64 * 137 + pad, 64 * 137 + pad, 256), ('glass', 64 * 37 + pad, 64 * 37 + pad, 256)])
    if 'B' in parts: part_b(decs, [('tissue', 64 * 157 + pad, 64 * 237 + pad, 1024), ('tissue edge', 64 * 137 + pad, 64 * 137 + pad, 1024)])
    if 'C' in parts: part_c(decs, rows, g)
    print('\ntotal %.0fs' % (time.time() - t0))
