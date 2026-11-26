# Strip-based halo recomputation: iSyntax (global 5/3) coefficients -> per-frame JPEG2000-convention 5/3 coefficients.
# Read-only on the input. Per channel: decode, then for each test frame:
#   copy interior coefficients (stored H bands + global LL arrays), recompute ONLY a ring of depth h at every frame edge
#   and level from strips (mirror extension), then (a) compare with a full per-frame forward transform,
#   (b) run a standard J2K inverse on the assembled coefficients and compare with the iSyntax-decoded pixels.
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
def lift_axis(Lc, Hc, axis):                       # iSyntax inverse (copied from reconstruct_np.py)
    Lc = np.moveaxis(Lc, axis, -1).astype(np.int32); Hc = np.moveaxis(Hc, axis, -1).astype(np.int32); N = Lc.shape[-1]
    Hn = np.concatenate([Hc[..., 1:], Hc[..., -1:]], axis=-1); xL = Lc - ((Hc+Hn+2) >> 2)
    xLp = np.concatenate([xL[..., :1], xL[..., :-1]], axis=-1); xH = Hc + ((xLp+xL) >> 1)
    out = np.empty(Lc.shape[:-1] + (2*N,), dtype=np.int32); out[..., 0::2] = xH; out[..., 1::2] = xL
    return np.moveaxis(out, -1, axis)
# ---- standard JPEG 2000 reversible 5/3 (T.800): even = low, whole-sample mirror at the ends
def fwd_axis(x, axis):
    x = np.moveaxis(x, axis, 0).astype(np.int32); ev, od = x[0::2], x[1::2]
    evn = np.concatenate([ev[1:], ev[-1:]], axis=0); d = od - ((ev + evn) >> 1)
    dp = np.concatenate([d[:1], d[:-1]], axis=0); s = ev + ((dp + d + 2) >> 2)
    return np.moveaxis(s, 0, axis), np.moveaxis(d, 0, axis)
def inv_axis(s, d, axis):
    s = np.moveaxis(s, axis, 0).astype(np.int32); d = np.moveaxis(d, axis, 0).astype(np.int32)
    dp = np.concatenate([d[:1], d[:-1]], axis=0); ev = s - ((dp + d + 2) >> 2)
    evn = np.concatenate([ev[1:], ev[-1:]], axis=0); od = d + ((ev + evn) >> 1)
    out = np.empty((2*s.shape[0],) + s.shape[1:], dtype=np.int32); out[0::2] = ev; out[1::2] = od
    return np.moveaxis(out, 0, axis)
def fwd2(x):                                       # vertical then horizontal
    lo, hi = fwd_axis(x, 0); LL, HL = fwd_axis(lo, 1); LH, HH = fwd_axis(hi, 1); return LL, HL, LH, HH
def inv2(LL, HL, LH, HH):                          # horizontal then vertical
    return inv_axis(inv_axis(LL, HL, 1), inv_axis(LH, HH, 1), 0)
def full_pyramid(fr, D):
    cur = fr; lv = []
    for _ in range(D):
        LL, HL, LH, HH = fwd2(cur); lv.append((HL, LH, HH)); cur = LL
    return cur, lv
def inv_pyramid(LLD, lv):
    cur = LLD
    for HL, LH, HH in reversed(lv): cur = inv2(cur, HL, LH, HH)
    return cur
# ---- strip halo recomputation
def build(fr, baseH, baseLL, D, h, stats):
    curf = fr.astype(np.int32); lv = []
    for j in range(D):
        n = curf.shape[0]; nn = n//2
        HLc, LHc, HHc = [a.copy() for a in baseH[j]]; LLc = baseLL[j].copy()
        if nn <= 2*h:                                              # ring covers everything: full forward
            LL, HL, LH, HH = fwd2(curf); LLc, HLc, LHc, HHc = LL.copy(), HL.copy(), LH.copy(), HH.copy(); stats['work'] += n*n
        else:
            W = 2*h + 4
            def strip(win_rows=None, win_cols=None):
                w = curf
                if win_cols is not None: w = w[:, win_cols]
                if win_rows is not None: w = w[win_rows, :]
                stats['work'] += w.size; return fwd2(w)
            outs = [HLc, LHc, HHc, LLc]
            def put(res, rs, cs, t_r, t_c):                        # res=(LL,HL,LH,HH) strip results; copy valid part
                for dst, src in zip((LLc, HLc, LHc, HHc), res): dst[rs, cs] = src[t_r, t_c]
            r = strip(win_cols=slice(0, W));            put(r, slice(0, nn), slice(0, h), slice(0, nn), slice(0, h))                    # left
            r = strip(win_cols=slice(n-W, n));          put(r, slice(0, nn), slice(nn-h, nn), slice(0, nn), slice(2, h+2))              # right
            r = strip(win_rows=slice(0, W));            put(r, slice(0, h), slice(0, nn), slice(0, h), slice(0, nn))                    # top
            r = strip(win_rows=slice(n-W, n));          put(r, slice(nn-h, nn), slice(0, nn), slice(2, h+2), slice(0, nn))              # bottom
        lv.append((HLc, LHc, HHc)); curf = LLc
        stats['total'] += n*n
    return curf, lv
# ---- main
FR = [('tissue 256', 237, 157, 256), ('tissue 512', 237, 157, 512), ('tissue edge 256', 137, 137, 256),
      ('blank glass 256', 37, 37, 256), ('scan-area edge 256', 385, 197, 256)]
D = 5; HS = [1, 2, 4]
allok = True
for c, cname in enumerate(('Y', 'Co', 'Cg')):
    rows = [dict(x=int(r['x']), y=int(r['y']), s=int(r['scale']), k=int(r['coef']), o=int(r['offset']), z=int(r['size']))
            for r in csv.DictReader(open(csvp)) if int(r['color']) == c]
    H = {s: np.zeros((3, need(Hp, s), need(Wp, s)), dtype=np.int16) for s in range(L)}
    LLtop = np.zeros((need(Hp, top), need(Wp, top)), dtype=np.int16); f = open(path, 'rb')
    for r in rows:
        if r['z'] <= 8: continue
        s, k = r['s'], r['k']; off = (3 << s) - 2 + ((1 << s) if k == 0 else 0)
        tx = (r['x']-ox-off)//(256 << s); ty = (r['y']-oy-off)//(256 << s)
        f.seek(r['o']); masks, planes, info = decode_block_v1(f.read(r['z']), 1 if k == 0 else 3)
        v = unpack(masks, planes, 1 if k == 0 else 3); tgt = LLtop[None] if k == 0 else H[s]
        y0, x0 = ty*128, tx*128; hh = min(128, tgt.shape[1]-y0); ww = min(128, tgt.shape[2]-x0)
        if hh > 0 and ww > 0: tgt[:, y0:y0+hh, x0:x0+ww] = v[:, :hh, :ww]
    cur = LLtop.astype(np.int32); curs = {}
    for s in range(top, -1, -1):
        n_r, n_c = need(Hp, s), need(Wp, s); cur = cur[:n_r, :n_c]; curs[s] = cur     # curs[s] = global LL_{s+1} (input of step s)
        out = lift_axis(lift_axis(cur, H[s][0][:n_r, :n_c], 1), lift_axis(H[s][1][:n_r, :n_c], H[s][2][:n_r, :n_c], 1), 0)
        nxt = np.zeros_like(out); nxt[1:, 1:] = out[:-1, :-1]; cur = nxt
    print('\n=== channel %s: decoded in %.0fs' % (cname, time.time()-t0), flush=True)
    for name, a, b, F in FR:
        x0, y0 = 1 + 64*a, 1 + 64*b; qx, qy = x0 + pad, y0 + pad
        fr = out[qy-1:qy-1+F, qx-1:qx-1+F]                          # padded pixel p = out[p-1]
        baseH = []; baseLL = []
        for j in range(D):
            s = j; nn = F >> (j+1); hx = ((qx+2) >> (s+1)) - 1 + np.arange(nn); hy = ((qy+2) >> (s+1)) - 1 + np.arange(nn)
            baseH.append((H[s][0][np.ix_(hy-1, hx)], H[s][1][np.ix_(hy, hx-1)], H[s][2][np.ix_(hy, hx)]))
            baseLL.append(curs[s][np.ix_(hy-1, hx-1)])
        refLL, refLV = full_pyramid(fr, D)
        line = []
        for h in HS:
            st = dict(work=0, total=0); LLb, lvb = build(fr, baseH, baseLL, D, h, st)
            nbad = int((LLb != refLL).sum()) + sum(int((a1 != a2).sum()) for l1, l2 in zip(lvb, refLV) for a1, a2 in zip(l1, l2))
            rec = inv_pyramid(LLb, lvb); pix_ok = bool((rec == fr).all())
            line.append('h=%d: coeff mismatches %d, pixels exact %s, strip work %.0f%% of full forward' % (h, nbad, pix_ok, 100*st['work']/st['total']))
            if h == 4 and (nbad or not pix_ok): allok = False
        # how many coefficients the strips actually changed vs plain copy (the true halo), at h=4
        LLb, lvb = build(fr, baseH, baseLL, D, 4, dict(work=0, total=0))
        ch = sum(int((lb[i] != baseH[j][i]).sum()) for j, lb in enumerate(lvb) for i in range(3)); tot = sum(a.size for l in baseH for a in l)
        print(' %-19s | copy-only would be wrong at %d/%d detail coeffs (%.2f%%) | %s' % (name, ch, tot, 100*ch/tot, ' | '.join(line)), flush=True)
print('\nALL OK (h=4):', allok)
