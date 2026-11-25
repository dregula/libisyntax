# Independent numpy reconstruction of a v1 .isyntax WSI (read-only; writes only to the output directory given).
# usage: reconstruct_np.py <file.isyntax> <blocks.csv> <outdir> [downsample=32]
# Geometry (from Philips MATLAB generate_struct / validwaveletcoefficientposition, verified in Philips_Compression_Code_Samples_Claude.md):
#   scale s: H coefficient n sits at padded-image position (3<<s)-2 + 2^(s+1) n ; LL coefficient n at that + 2^s
#   per level: out[m]: even m <- H (inverse of d-step), odd m <- L ; L_{s-1}[n] = out[n-1]
import sys, os, csv, re, struct, zlib, time, collections
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hulsken_v1 import decode_block_v1
path, csvp, outdir = sys.argv[1], sys.argv[2], sys.argv[3]
DS = int(sys.argv[4]) if len(sys.argv) > 4 else 32
CROP = tuple(int(v) for v in sys.argv[5:9]) if len(sys.argv) >= 9 else None   # x y w h in valid-region pixels
os.makedirs(outdir, exist_ok=True)
t0 = time.time()
head = open(path, 'rb').read(8 * 1024 * 1024); eot = head.find(b'\r\n\x04'); xml = head[:eot].decode('utf-8')
gh = xml[xml.index('UFS_IMAGE_GENERAL_HEADERS'):]
gr = [tuple(map(int, r.split())) for r in re.findall(r'Name="UFS_IMAGE_DIMENSION_RANGE"[^>]*>([^<]*)<', gh)][:5]
Wp = (gr[0][2] - gr[0][0] + gr[0][1]) // gr[0][1]; Hp = (gr[1][2] - gr[1][0] + gr[1][1]) // gr[1][1]
ox, oy = gr[0][0], gr[1][0]; L = (gr[3][2] - gr[3][0] + gr[3][1]) // gr[3][1]; top = L - 1
pad = (3 << L) - 3
print(f'padded {Wp}x{Hp} levels {L} origin ({ox},{oy}) pad {pad} valid {Wp-2*pad}x{Hp-2*pad}')
rows = [dict(x=int(r['x']), y=int(r['y']), c=int(r['color']), s=int(r['scale']), k=int(r['coef']), o=int(r['offset']), z=int(r['size'])) for r in csv.DictReader(open(csvp))]
# permutation: serialized (4x2 cell) order -> raster 128x128
P = np.arange(16384); area = P // 16; within = P % 16
perm = ((area // 32) * 4 + within // 4) * 128 + (area % 32) * 4 + within % 4      # out_flat[perm[i]] = buf[i]
def unpack(masks, planes, ncoef):
    arr = np.frombuffer(planes, dtype=np.uint8); res = np.zeros((ncoef, 16384), dtype=np.uint16); p = 0
    for k in range(ncoef):
        for i in range(16):
            if (masks[k] >> i) & 1:
                bits = np.unpackbits(arr[p * 2048:(p + 1) * 2048], bitorder='little').astype(np.uint16); p += 1
                res[k] |= bits << (15 if i == 0 else i - 1)
    out = np.zeros_like(res)
    out[:, perm] = res
    mag = (out & 0x7FFF).astype(np.int32); sgn = (out >> 15).astype(bool)
    return np.where(sgn, -mag, mag).astype(np.int16).reshape(ncoef, 128, 128)
def need(n, s): return -(-n // (1 << (s + 1))) + 2          # coefficient count needed at scale s
f = open(path, 'rb')
H = {}; LLtop = {}
for c in range(3):
    for s in range(L):
        H[(c, s)] = np.zeros((3, need(Hp, s), need(Wp, s)), dtype=np.int16)
    LLtop[c] = np.zeros((need(Hp, top), need(Wp, top)), dtype=np.int16)
nblk = 0
for r in rows:
    c, s, k = r['c'], r['s'], r['k']; ncoef = 1 if k == 0 else 3
    off = (3 << s) - 2 + ((1 << s) if k == 0 else 0)
    tx = (r['x'] - ox - off) // (256 << s); ty = (r['y'] - oy - off) // (256 << s)
    if r['z'] <= 8: continue
    f.seek(r['o']); data = f.read(r['z'])
    masks, planes, info = decode_block_v1(data, ncoef)
    v = unpack(masks, planes, ncoef); nblk += 1
    tgt = LLtop[c][None] if k == 0 else H[(c, s)]
    y0, x0 = ty * 128, tx * 128
    hh = min(128, tgt.shape[1] - y0); ww = min(128, tgt.shape[2] - x0)
    if hh > 0 and ww > 0: tgt[:, y0:y0 + hh, x0:x0 + ww] = v[:, :hh, :ww]
print(f'decoded+placed {nblk} blocks in {time.time()-t0:.0f}s')
stats = dict(maxabs=0, overflow16=0)
def lift_axis(Lc, Hc, axis):
    """inverse 5/3 along `axis`; Lc, Hc same shape (..., N, ...). returns array with 2N samples along axis; even idx<-H, odd<-L."""
    Lc = np.moveaxis(Lc, axis, -1).astype(np.int32); Hc = np.moveaxis(Hc, axis, -1).astype(np.int32)
    N = Lc.shape[-1]
    Hn = np.concatenate([Hc[..., 1:], Hc[..., -1:]], axis=-1)                   # H[n+1], mirrored at the end
    xL = Lc - ((Hc + Hn + 2) >> 2)
    xLp = np.concatenate([xL[..., :1], xL[..., :-1]], axis=-1)                   # xL[n-1], mirrored at the start
    xH = Hc + ((xLp + xL) >> 1)
    out = np.empty(Lc.shape[:-1] + (2 * N,), dtype=np.int32)
    out[..., 0::2] = xH; out[..., 1::2] = xL
    stats['maxabs'] = max(stats['maxabs'], int(np.abs(out).max()))
    return np.moveaxis(out, -1, axis)
planes_full = []
for c in range(3):
    cur = LLtop[c].astype(np.int32)
    for s in range(top, -1, -1):
        HL, LH, HH = H[(c, s)][0], H[(c, s)][1], H[(c, s)][2]
        n_r, n_c = need(Hp, s), need(Wp, s)
        cur = cur[:n_r, :n_c]
        lowrows = lift_axis(cur, HL[:n_r, :n_c], axis=1)            # horizontal on LL|HL  (vertical-low rows)
        highrows = lift_axis(LH[:n_r, :n_c], HH[:n_r, :n_c], axis=1)  # horizontal on LH|HH  (vertical-high rows)
        out = lift_axis(lowrows, highrows, axis=0)                  # vertical: first arg = low (odd rows), second = high (even rows)
        # next L: L[n] = out[n-1]
        nxt = np.zeros_like(out); nxt[1:, 1:] = out[:-1, :-1]
        cur = nxt
        print(f'  colour {c} scale {s}: out {out.shape} t={time.time()-t0:.0f}s', flush=True)
        if s == 0: final = out
    # pixels: pixel[p] = out[p-1]
    pix = np.zeros((Hp, Wp), dtype=np.int32); pix[1:, 1:] = final[:Hp - 1, :Wp - 1]
    planes_full.append(pix[pad:Hp - pad, pad:Wp - pad].copy()); del final, pix, out, cur
Y, Co, Cg = planes_full
h, w = Y.shape
if os.environ.get('SAVE_PLANES'):
    for nm, arr in (('Y', Y), ('Co', Co), ('Cg', Cg)): np.save(os.path.join(outdir, nm + '.npy'), arr.astype(np.int16))
    print('saved Y/Co/Cg planes (int16) to', outdir)
print('valid region', w, 'x', h, '| max |intermediate| over all levels:', stats['maxabs'], '(int16 limit 32767)')
def png(path, arr):
    hh, ww, _ = arr.shape; raw = b''.join(b'\x00' + arr[i].tobytes() for i in range(hh))
    def chunk(t, d): return struct.pack('>I', len(d)) + t + d + struct.pack('>I', zlib.crc32(t + d) & 0xffffffff)
    open(path, 'wb').write(b'\x89PNG\r\n\x1a\n' + chunk(b'IHDR', struct.pack('>IIBBBBB', ww, hh, 8, 2, 0, 0, 0)) + chunk(b'IDAT', zlib.compress(raw, 6)) + chunk(b'IEND', b''))
def rgb(Yv, Cov, Cgv, mode):
    sh = (lambda a: a >> 1) if mode == 'floor' else (lambda a: np.trunc(a / 2).astype(np.int32))
    tmp = Yv - sh(Cgv); G = tmp + Cgv; B = tmp - sh(Cov); R = B + Cov
    return R, G, B
acc = collections.Counter(); sums = np.zeros(3); dmin = np.array([9, 9, 9]); dmax = np.array([-9, -9, -9])
band = 1024; prev_rows = []; crop = None
ph = (h // DS) * DS; pw = (w // DS) * DS
for y0 in range(0, h, band):
    sl = slice(y0, min(h, y0 + band)); Yb, Cob, Cgb = Y[sl], Co[sl], Cg[sl]
    acc['pixels'] += Yb.size; acc['Y_negative'] += int((Yb < 0).sum())
    Ya = np.abs(Yb)
    Rf, Gf, Bf = rgb(Ya, Cob, Cgb, 'floor'); Rt, Gt, Bt = rgb(Ya, Cob, Cgb, 'trunc'); Rn, Gn, Bn = rgb(Yb, Cob, Cgb, 'floor')
    dd = [Rf - Rt, Gf - Gt, Bf - Bt]
    acc['floor_vs_trunc_differ'] += int(((dd[0] != 0) | (dd[1] != 0) | (dd[2] != 0)).sum())
    for k in range(3): sums[k] += np.abs(dd[k]).sum(); dmin[k] = min(dmin[k], dd[k].min()); dmax[k] = max(dmax[k], dd[k].max())
    acc['rgb_below0'] += int(((Rf < 0) | (Gf < 0) | (Bf < 0)).sum()); acc['rgb_above255'] += int(((Rf > 255) | (Gf > 255) | (Bf > 255)).sum())
    acc['abs_Y_changes_rgb'] += int(((Rn != Rf) | (Gn != Gf) | (Bn != Bf)).sum())
    acc['neg_odd_Cg'] += int(((Cgb < 0) & ((Cgb & 1) == 1)).sum()); acc['neg_odd_Co'] += int(((Cob < 0) & ((Cob & 1) == 1)).sum())
    acc['Cg_neg'] += int((Cgb < 0).sum()); acc['Co_neg'] += int((Cob < 0).sum())
    # preview rows (box average) for complete DS-row groups
    rgbf = np.stack([np.clip(X, 0, 255) for X in (Rf, Gf, Bf)], axis=-1).astype(np.float32)
    if y0 % DS == 0 and y0 + band <= ph:
        pv = rgbf[:, :pw].reshape(band // DS, DS, pw // DS, DS, 3).mean(axis=(1, 3)); prev_rows.append(pv)
    if CROP:
        cx, cy, cw, ch = CROP
        if y0 <= cy < y0 + band:
            part = rgbf[cy - y0:cy - y0 + min(ch, band - (cy - y0)), cx:cx + cw].astype(np.uint8); crop = part if crop is None else np.concatenate([crop, part])
        elif crop is not None and crop.shape[0] < ch:
            part = rgbf[:ch - crop.shape[0], cx:cx + cw].astype(np.uint8); crop = np.concatenate([crop, part])
print('FINDINGS (valid region):')
print(f"  pixels {acc['pixels']}; Y min {int(Y.min())} max {int(Y.max())}; Co range [{int(Co.min())},{int(Co.max())}]; Cg range [{int(Cg.min())},{int(Cg.max())}]")
print(f"  negative Y after IDWT: {acc['Y_negative']}  (abs(Y) changes RGB of {acc['abs_Y_changes_rgb']} pixels)")
print(f"  RGB before clamping: below 0: {acc['rgb_below0']}, above 255: {acc['rgb_above255']} pixels")
print(f"  floor (>>1, SIMD path) vs truncate (/2, scalar path): {acc['floor_vs_trunc_differ']} pixels differ = {100*acc['floor_vs_trunc_differ']/acc['pixels']:.1f}%")
print('   difference floor-minus-trunc range per channel R/G/B:', list(zip(dmin.tolist(), dmax.tolist())), '| mean |diff|:', [round(v / acc['pixels'], 3) for v in sums])
print(f"  fraction with negative Cg {acc['Cg_neg']/acc['pixels']:.3f} (odd {acc['neg_odd_Cg']/acc['pixels']:.3f}), negative Co {acc['Co_neg']/acc['pixels']:.3f} (odd {acc['neg_odd_Co']/acc['pixels']:.3f})")
# tile-seam test: mean |neighbour difference| of Y by position modulo the 256 px tile pitch (padded coordinates)
for axis, name in ((1, 'columns'), (0, 'rows')):
    Ys = Y[::16] if axis == 1 else Y[:, ::16]
    d = np.abs(np.diff(Ys, axis=axis)); n = d.shape[axis]
    idx = (np.arange(n) + 1 + pad) % 256
    d2 = np.moveaxis(d, axis, -1).reshape(-1, n).mean(axis=0)
    prof = np.array([d2[idx == r].mean() for r in range(256)])
    print(f'  seam profile ({name}): mean |dY| by position mod 256: min {prof.min():.4f} median {np.median(prof):.4f} max {prof.max():.4f}; max/median {prof.max()/np.median(prof):.3f} (argmax residue {int(prof.argmax())})')
prev = np.clip(np.concatenate(prev_rows), 0, 255).round().astype(np.uint8)
png(os.path.join(outdir, f'recon_ds{DS}.png'), prev)
print('wrote preview', prev.shape)
if crop is not None: png(os.path.join(outdir, f'crop_{CROP[0]}_{CROP[1]}.png'), crop); print('wrote crop', crop.shape)
print(f'total {time.time()-t0:.0f}s')
