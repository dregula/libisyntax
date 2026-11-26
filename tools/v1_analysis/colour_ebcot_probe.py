# Colour-transform and EBCOT probe for a v1 iSyntax file (read-only; needs numpy, ~10 GB RAM).
#   python3 -B -I colour_ebcot_probe.py <file.isyntax> <blocks.csv>
# Part A: wavelet-domain YCoCg -> YBR_RCT conversion (Cb = -q(Co) - Cg, Cr = Cb + Co, Y' = Y) on real frames, versus the
#         pixel-domain conversion, for both rounding conventions of the Philips inverse (fix / floor) and both q.
# Part B: whole-file statistics per scale of the stored coefficients and of the converted Cb/Cr: non-zero fraction, empty
#         least-significant planes (quantizer), most-significant bit count, zero-order entropy (a crude code-length proxy,
#         NOT an EBCOT size), and the share of all-zero 64 x 64 blocks (what EBCOT codes almost for free).
import sys, os, time, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hulsken_v1 import decode_block_v1
from windowed_idwt_halo import open_file, unpack, need, full_pyramid, inv_pyramid, Windowed

fl = lambda x: x >> 1
def tr(x): return np.where(x >= 0, x >> 1, -((-x) >> 1))
def decode_coeffs(path, rows, g, c):
    Hp, Wp, top = g['Hp'], g['Wp'], g['top']
    H = {s: np.zeros((3, need(Hp, s), need(Wp, s)), dtype=np.int16) for s in range(g['L'])}
    LLtop = np.zeros((need(Hp, top), need(Wp, top)), dtype=np.int16); f = open(path, 'rb')
    for r in rows:
        if r['c'] != c or r['z'] <= 8: continue
        s, k = r['s'], r['k']; off = (3 << s) - 2 + ((1 << s) if k == 0 else 0)
        tx = (r['x'] - g['ox'] - off) // (256 << s); ty = (r['y'] - g['oy'] - off) // (256 << s)
        f.seek(r['o']); masks, planes, info = decode_block_v1(f.read(r['z']), 1 if k == 0 else 3)
        v = unpack(masks, planes, 1 if k == 0 else 3); tgt = LLtop[None] if k == 0 else H[s]
        y0, x0 = ty * 128, tx * 128; hh = min(128, tgt.shape[1] - y0); ww = min(128, tgt.shape[2] - x0)
        if hh > 0 and ww > 0: tgt[:, y0:y0 + hh, x0:x0 + ww] = v[:, :hh, :ww]
    return H, LLtop

def inv_rgb(Y, Co, Cg, h):
    t = Y - h(Cg); G = t + Cg; B = t - h(Co); return B + Co, G, B
def rct_fwd(R, G, B): return (R + 2 * G + B) >> 2, B - G, R - G
def rct_inv(Y, Cb, Cr):
    G = Y - ((Cb + Cr) >> 2); return Cr + G, G, Cb + G

def convert_pyramid(py, q):                          # py = [(LL, levels)] for Y, Co, Cg ; returns pyramids of Y', Cb, Cr
    (LLy, lvy), (LLo, lvo), (LLg, lvg) = py
    cb = lambda Co, Cg: -q(Co) - Cg
    LLb = cb(LLo, LLg); LLr = LLb + LLo; lvb = []; lvr = []
    for (a, b, c), (d, e, f_) in zip(lvo, lvg):
        bb = [cb(x, y) for x, y in zip((a, b, c), (d, e, f_))]; lvb.append(tuple(bb)); lvr.append(tuple(bb[i] + (a, b, c)[i] for i in range(3)))
    return (LLy, lvy), (LLb, lvb), (LLr, lvr)

def part_a(decs, frames, D=5):
    print('\n== PART A: wavelet-domain colour conversion on real frames (per-frame mirror pyramids, D = %d)' % D)
    for name, i0y, i0x, F in frames:
        pix = [d.out_region(0, i0y, i0y + F, i0x, i0x + F) for d in decs]
        py = [full_pyramid(p, D) for p in pix]
        odd = []                                       # where does q() actually round? fraction of odd Co per level (detail bands)
        for j, (a, b, c) in enumerate(py[1][1]): odd.append('%.1f%%' % (100 * np.mean(np.concatenate([(x & 1).ravel() for x in (a, b, c)]))))
        print(' frame %-12s %dx%d | odd Co among detail coefficients, levels 1..%d: %s | LL%d: %.1f%% odd' % (
            name, F, F, D, ', '.join(odd), D, 100 * np.mean(py[1][0] & 1)))
        for qn, q in (('q=floor', fl), ('q=trunc', tr)):
            (LLy, lvy), cbp, crp = convert_pyramid(py, q)
            Yp = inv_pyramid(LLy, lvy); Cbp = inv_pyramid(*cbp); Crp = inv_pyramid(*crp); RGBb = rct_inv(Yp, Cbp, Crp)
            for cn, h in (('inverse=fix', tr), ('inverse=floor', fl)):
                R, G, B = inv_rgb(pix[0], pix[1], pix[2], h); Yr, Cbr, Crr = rct_fwd(R, G, B)
                dRGB = np.stack([a - b for a, b in zip(RGBb, (R, G, B))]); anyd = (dRGB != 0).any(0)
                print('   %s vs pixel-domain %-13s: Y\' differs %5.1f%% | Cb differs %5.1f%% (max %d) | Cr differs %5.1f%% (max %d) | RGB differs %5.1f%% (max |d| %d)' % (
                    qn, cn, 100 * np.mean(Yp != Yr), 100 * np.mean(Cbp != Cbr), np.abs(Cbp - Cbr).max(), 100 * np.mean(Crp != Crr),
                    np.abs(Crp - Crr).max(), 100 * anyd.mean(), np.abs(dRGB).max()))

def h0_bits(arr, hist):
    n = arr.size; p = hist[hist > 0] / n; return float(-(p * np.log2(p)).sum() * n)
def band_stats(arr):
    a = arr.astype(np.int32); nz = int((a != 0).sum()); ab = np.abs(a)
    orr = int(np.bitwise_or.reduce(ab.ravel())) if nz else 0
    low = (orr & -orr).bit_length() - 1 if orr else 0; msb = int(ab.max()).bit_length()
    hist = np.bincount((a + 40000).ravel(), minlength=80001)
    r64, c64 = a.shape[0] // 64, a.shape[1] // 64
    blk = (a[:r64 * 64, :c64 * 64] != 0).reshape(r64, 64, c64, 64).any(axis=(1, 3)); zb = int((~blk).sum()); nb = blk.size
    return dict(n=a.size, nz=nz, low=low, msb=msb, bits=h0_bits(a, hist), zb=zb, nb=nb)

def part_b(Hs, LLs, top):
    print('\n== PART B: whole-file statistics per scale (detail bands HL+LH+HH together; scale 5 also lists the top LL)')
    names = ['Y', 'Co', 'Cg', 'Cb(q=floor)', 'Cr(q=floor)', 'Cb(q=trunc)', 'Cr(q=trunc)']
    tot = {n: 0.0 for n in names}
    print(' %-5s %-12s %8s %9s %7s %5s %12s %14s' % ('scale', 'channel', 'nonzero%', 'empty LSB', 'MSB bits', '', 'H0 bits/coef', 'all-zero 64x64'))
    for s in list(range(top + 1)) + ['LL']:
        acc = {n: dict(n=0, nz=0, low=99, msb=0, bits=0.0, zb=0, nb=0) for n in names}
        bands = [(s, b) for b in range(3)] if s != 'LL' else [('LL', 0)]
        for sc, b in bands:
            Y, Co, Cg = [(Hs[c][sc][b] if sc != 'LL' else LLs[c]).astype(np.int32) for c in range(3)]
            Cbf = -(Co >> 1) - Cg; Cbt = -tr(Co) - Cg
            arrs = {'Y': Y, 'Co': Co, 'Cg': Cg, 'Cb(q=floor)': Cbf, 'Cr(q=floor)': Cbf + Co, 'Cb(q=trunc)': Cbt, 'Cr(q=trunc)': Cbt + Co}
            for n in names:
                st = band_stats(arrs[n]); a = acc[n]
                for k in ('n', 'nz', 'bits', 'zb', 'nb'): a[k] += st[k]
                a['low'] = min(a['low'], st['low'] if st['nz'] else 99); a['msb'] = max(a['msb'], st['msb'])
            del Y, Co, Cg, Cbf, Cbt, arrs
        for n in names:
            a = acc[n]; tot[n] += a['bits']
            print(' %-5s %-12s %8.2f %9s %7d %5s %12.3f %13.1f%%' % (s, n, 100 * a['nz'] / a['n'], a['low'] if a['low'] != 99 else '-', a['msb'], '', a['bits'] / a['n'], 100 * a['zb'] / a['nb']), flush=True)
    print('\n zero-order entropy total over all stored coefficients (MB; crude proxy, no context modelling, includes the padded margin):')
    for n in names: print('   %-12s %8.2f MB' % (n, tot[n] / 8e6))
    print('   Co + Cg = %.2f MB | Cb + Cr (q=floor) = %.2f MB | Cb + Cr (q=trunc) = %.2f MB' % ((tot['Co'] + tot['Cg']) / 8e6,
          (tot['Cb(q=floor)'] + tot['Cr(q=floor)']) / 8e6, (tot['Cb(q=trunc)'] + tot['Cr(q=trunc)']) / 8e6))

if __name__ == '__main__':
    path, csvp = sys.argv[1], sys.argv[2]; t0 = time.time()
    g, rows = open_file(path, csvp); top, pad = g['top'], g['pad']
    Hs, LLs = {}, {}
    for c in range(3):
        Hs[c], LLs[c] = decode_coeffs(path, rows, g, c); print('decoded channel %d (%.0fs)' % (c, time.time() - t0), flush=True)
    decs = [Windowed(Hs[c], LLs[c], top, g['Hp'], g['Wp']) for c in range(3)]
    frames = [('tissue', 64 * 157 + pad, 64 * 237 + pad, 256), ('tissue edge', 64 * 137 + pad, 64 * 137 + pad, 256), ('tissue 1024', 64 * 157 + pad, 64 * 237 + pad, 1024)]
    part_a(decs, frames)
    part_b(Hs, {c: LLs[c] for c in range(3)}, top)
    print('\ntotal %.0fs' % (time.time() - t0))
