# Non-circular test: stored iSyntax detail coefficients (Y channel, real file) vs a standard-J2K per-frame 5/3 forward
# transform (T.800 order: vertical then horizontal, low at even local coordinates, whole-sample mirror at frame edges)
# of the decoded Y plane. Reads the input read-only; writes nothing.
import sys, re, csv, time, numpy as np
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hulsken_v1 import decode_block_v1
path, csvp = sys.argv[1], sys.argv[2]
t0 = time.time()
head = open(path, 'rb').read(8*1024*1024); eot = head.find(b'\r\n\x04'); xml = head[:eot].decode('utf-8')
gh = xml[xml.index('UFS_IMAGE_GENERAL_HEADERS'):]
gr = [tuple(map(int, r.split())) for r in re.findall(r'Name="UFS_IMAGE_DIMENSION_RANGE"[^>]*>([^<]*)<', gh)][:5]
Wp = (gr[0][2]-gr[0][0]+gr[0][1])//gr[0][1]; Hp = (gr[1][2]-gr[1][0]+gr[1][1])//gr[1][1]
ox, oy = gr[0][0], gr[1][0]; L = (gr[3][2]-gr[3][0]+gr[3][1])//gr[3][1]; top = L-1; pad = (3 << L) - 3
P = np.arange(16384); area = P//16; within = P % 16
perm = ((area//32)*4 + within//4)*128 + (area % 32)*4 + within % 4
def unpack(masks, planes, ncoef):
    arr = np.frombuffer(planes, dtype=np.uint8); res = np.zeros((ncoef, 16384), dtype=np.uint16); p = 0
    for k in range(ncoef):
        for i in range(16):
            if (masks[k] >> i) & 1:
                bits = np.unpackbits(arr[p*2048:(p+1)*2048], bitorder='little').astype(np.uint16); p += 1
                res[k] |= bits << (15 if i == 0 else i-1)
    out = np.zeros_like(res); out[:, perm] = res
    mag = (out & 0x7FFF).astype(np.int32); sgn = (out >> 15).astype(bool)
    return np.where(sgn, -mag, mag).astype(np.int16).reshape(ncoef, 128, 128)
need = lambda n, s: -(-n // (1 << (s+1))) + 2
rows = [dict(x=int(r['x']), y=int(r['y']), s=int(r['scale']), k=int(r['coef']), o=int(r['offset']), z=int(r['size']))
        for r in csv.DictReader(open(csvp)) if int(r['color']) == 0]
H = {s: np.zeros((3, need(Hp, s), need(Wp, s)), dtype=np.int16) for s in range(L)}
LLtop = np.zeros((need(Hp, top), need(Wp, top)), dtype=np.int16)
f = open(path, 'rb')
for r in rows:
    if r['z'] <= 8: continue
    s, k = r['s'], r['k']; off = (3 << s) - 2 + ((1 << s) if k == 0 else 0)
    tx = (r['x']-ox-off)//(256 << s); ty = (r['y']-oy-off)//(256 << s)
    f.seek(r['o']); masks, planes, info = decode_block_v1(f.read(r['z']), 1 if k == 0 else 3)
    v = unpack(masks, planes, 1 if k == 0 else 3)
    tgt = LLtop[None] if k == 0 else H[s]
    y0, x0 = ty*128, tx*128; hh = min(128, tgt.shape[1]-y0); ww = min(128, tgt.shape[2]-x0)
    if hh > 0 and ww > 0: tgt[:, y0:y0+hh, x0:x0+ww] = v[:, :hh, :ww]
print('decoded Y blocks %.0fs' % (time.time()-t0), flush=True)
def lift_axis(Lc, Hc, axis):                       # copied from reconstruct_np.py (iSyntax inverse)
    Lc = np.moveaxis(Lc, axis, -1).astype(np.int32); Hc = np.moveaxis(Hc, axis, -1).astype(np.int32); N = Lc.shape[-1]
    Hn = np.concatenate([Hc[..., 1:], Hc[..., -1:]], axis=-1); xL = Lc - ((Hc+Hn+2) >> 2)
    xLp = np.concatenate([xL[..., :1], xL[..., :-1]], axis=-1); xH = Hc + ((xLp+xL) >> 1)
    out = np.empty(Lc.shape[:-1] + (2*N,), dtype=np.int32); out[..., 0::2] = xH; out[..., 1::2] = xL
    return np.moveaxis(out, -1, axis)
cur = LLtop.astype(np.int32)
for s in range(top, -1, -1):
    HL, LH, HH = H[s][0], H[s][1], H[s][2]; n_r, n_c = need(Hp, s), need(Wp, s); cur = cur[:n_r, :n_c]
    out = lift_axis(lift_axis(cur, HL[:n_r, :n_c], 1), lift_axis(LH[:n_r, :n_c], HH[:n_r, :n_c], 1), 0)
    nxt = np.zeros_like(out); nxt[1:, 1:] = out[:-1, :-1]; cur = nxt
pix = np.zeros((Hp, Wp), dtype=np.int32); pix[1:, 1:] = out[:Hp-1, :Wp-1]
Y = pix[pad:Hp-pad, pad:Wp-pad]; print('Y plane', Y.shape, 'Y min/max', Y.min(), Y.max(), '%.0fs' % (time.time()-t0), flush=True)
# ---- standard J2K reversible 5/3, one level, one axis (even = low), whole-sample symmetric extension
def j2k_axis(x, axis):
    x = np.moveaxis(x, axis, 0).astype(np.int32); ev, od = x[0::2], x[1::2]
    evn = np.concatenate([ev[1:], ev[-1:]], axis=0)                 # x[2n+2], mirrored: x[N] = x[N-2] = ev[-1]
    d = od - ((ev + evn) >> 1); dp = np.concatenate([d[:1], d[:-1]], axis=0)   # d[-1] = d[0]
    s = ev + ((dp + d + 2) >> 2)
    return np.moveaxis(s, 0, axis), np.moveaxis(d, 0, axis)
def j2k_level(x):                                                   # T.800 2D_SD: vertical first, then horizontal
    lo, hi = j2k_axis(x, 0)
    LL, HL = j2k_axis(lo, 1); LH, HH = j2k_axis(hi, 1)
    return LL, HL, LH, HH
def test(x0, y0, F, D, label):
    fr = Y[y0:y0+F, x0:x0+F].astype(np.int32); q_x = x0 + pad; q_y = y0 + pad; cur = fr; res = []
    for j in range(D):
        s = j; LL, HL, LH, HH = j2k_level(cur); cur = LL; n = F >> (j+1)
        ix = (q_x + 2) // (1 << (s+1)) - 1 + np.arange(n); iy = (q_y + 2) // (1 << (s+1)) - 1 + np.arange(n)
        tot = 0; bad = 0; widest = 0
        for b, loc in enumerate((HL, LH, HH)):
            ry = iy - 1 if b == 0 else iy          # HL: vertical-low rows use the low-type index (one less)
            cx = ix - 1 if b == 1 else ix          # LH: horizontal-low columns use the low-type index
            g = H[s][b][np.ix_(ry, cx)].astype(np.int32); m = (g != loc); tot += m.size; bad += int(m.sum())
            if m.any():
                ys, xs = np.nonzero(m); e = np.minimum.reduce([ys, n-1-ys, xs, n-1-xs]); widest = max(widest, int(e.max())+1)
        res.append('L%d: %d/%d differ (%.1f%%, edge depth<=%d)' % (j+1, bad, tot, 100*bad/tot, widest))
    print(label, 'F=%d D=%d | ' % (F, D) + ' | '.join(res))
a, b = 237, 157
x0, y0 = 1 + 64*a, 1 + 64*b
print('frame origin (valid coords) x0=%d y0=%d  (x0 mod 64 = %d)' % (x0, y0, x0 % 64))
for F in (256, 512): test(x0, y0, F, 5, 'ALIGNED (x0,y0 = 1 mod 64):')
test(x0+1, y0, 256, 5, 'x shifted +1:        ')
test(x0, y0+1, 256, 5, 'y shifted +1:        ')
test(x0-1, y0-1, 256, 5, 'both shifted -1:     ')
