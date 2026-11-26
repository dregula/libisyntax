# Windowed top-down IDWT + strip halo recomputation: iSyntax v1 (global 5/3) -> per-frame JPEG 2000 5/3 coefficients.
# Read-only on the input file; writes nothing. Needs numpy (repo .venv). usage:
#   python3 -I windowed_idwt_halo.py <file.isyntax> <blocks.csv>        (blocks.csv from parse_v1.py)
#
# What it shows (see iSyntax_to_DICOM_WSI_Transcoding_Claude.md):
#  * DICOM-WSI layer k (a separate SOP Instance) is the pixel plane the iSyntax decoder produces at scale k ("out_k").
#    Each of its frames is an independent JPEG 2000 codestream (D decomposition levels, mirror extension at the frame
#    edges), so the frame must carry its own LL_D band. Nothing refers to another layer.
#  * Frame coefficients = interior copied from the stored H blocks + a ring of depth h recomputed from strips.
#    The strips (pixel ring of out_k, rings of the global LL arrays, and the small LL_D) are produced by a WINDOWED
#    top-down IDWT that decodes only a rectangle around the frame edges, never the whole image.
#  * Everything the builder does not read is filled with random garbage, so an exact result proves it was not read.
#  * Beyond the top stored scale the layer's extra J2K levels are an ordinary local forward transform of the frame's LL.
import sys, os, re, csv, time, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hulsken_v1 import decode_block_v1

# ---------------------------------------------------------------- file / block decoding (as reconstruct_np.py)
def open_file(path, csvp):
    f_ = open(path, 'rb'); head = f_.read(8 * 1024 * 1024)
    while head.find(b'\r\n\x04') < 0:                 # the XML header of the larger files exceeds 8 MB
        more = f_.read(8 * 1024 * 1024)
        if not more: raise SystemExit('no EOT marker found')
        head += more
    eot = head.find(b'\r\n\x04'); xml = head[:eot].decode('utf-8')
    gh = xml[xml.index('UFS_IMAGE_GENERAL_HEADERS'):]
    gr = [tuple(map(int, r.split())) for r in re.findall(r'Name="UFS_IMAGE_DIMENSION_RANGE"[^>]*>([^<]*)<', gh)][:5]
    g = dict(Wp=(gr[0][2] - gr[0][0] + gr[0][1]) // gr[0][1], Hp=(gr[1][2] - gr[1][0] + gr[1][1]) // gr[1][1],
             ox=gr[0][0], oy=gr[1][0], L=(gr[3][2] - gr[3][0] + gr[3][1]) // gr[3][1])
    g['top'] = g['L'] - 1; g['pad'] = (3 << g['L']) - 3
    rows = [dict(x=int(r['x']), y=int(r['y']), c=int(r['color']), s=int(r['scale']), k=int(r['coef']),
                 o=int(r['offset']), z=int(r['size'])) for r in csv.DictReader(open(csvp))]
    return g, rows
P = np.arange(16384); _area = P // 16; _within = P % 16
PERM = ((_area // 32) * 4 + _within // 4) * 128 + (_area % 32) * 4 + _within % 4
def unpack(masks, planes, ncoef):
    arr = np.frombuffer(planes, dtype=np.uint8); res = np.zeros((ncoef, 16384), dtype=np.uint16); p = 0
    for k in range(ncoef):
        for i in range(16):
            if (masks[k] >> i) & 1:
                bits = np.unpackbits(arr[p * 2048:(p + 1) * 2048], bitorder='little').astype(np.uint16); p += 1
                res[k] |= bits << (15 if i == 0 else i - 1)
    out = np.zeros_like(res); out[:, PERM] = res
    mag = (out & 0x7FFF).astype(np.int32); sgn = (out >> 15).astype(bool)
    return np.where(sgn, -mag, mag).astype(np.int16).reshape(ncoef, 128, 128)
need = lambda n, s: -(-n // (1 << (s + 1))) + 2
def lift_axis(Lc, Hc, axis):                       # iSyntax inverse lifting, copied from reconstruct_np.py
    Lc = np.moveaxis(Lc, axis, -1).astype(np.int32); Hc = np.moveaxis(Hc, axis, -1).astype(np.int32); N = Lc.shape[-1]
    Hn = np.concatenate([Hc[..., 1:], Hc[..., -1:]], axis=-1); xL = Lc - ((Hc + Hn + 2) >> 2)
    xLp = np.concatenate([xL[..., :1], xL[..., :-1]], axis=-1); xH = Hc + ((xLp + xL) >> 1)
    out = np.empty(Lc.shape[:-1] + (2 * N,), dtype=np.int32); out[..., 0::2] = xH; out[..., 1::2] = xL
    return np.moveaxis(out, -1, axis)

# ---------------------------------------------------------------- standard JPEG 2000 reversible 5/3 (T.800 convention)
def fwd_axis(x, axis):
    x = np.moveaxis(x, axis, 0).astype(np.int32); ev, od = x[0::2], x[1::2]
    evn = np.concatenate([ev[1:], ev[-1:]], axis=0); d = od - ((ev + evn) >> 1)
    dp = np.concatenate([d[:1], d[:-1]], axis=0); s = ev + ((dp + d + 2) >> 2)
    return np.moveaxis(s, 0, axis), np.moveaxis(d, 0, axis)
def inv_axis(s, d, axis):
    s = np.moveaxis(s, axis, 0).astype(np.int32); d = np.moveaxis(d, axis, 0).astype(np.int32)
    dp = np.concatenate([d[:1], d[:-1]], axis=0); ev = s - ((dp + d + 2) >> 2)
    evn = np.concatenate([ev[1:], ev[-1:]], axis=0); od = d + ((ev + evn) >> 1)
    out = np.empty((2 * s.shape[0],) + s.shape[1:], dtype=np.int32); out[0::2] = ev; out[1::2] = od
    return np.moveaxis(out, 0, axis)
def fwd2(x):                                       # vertical, then horizontal
    lo, hi = fwd_axis(x, 0); LL, HL = fwd_axis(lo, 1); LH, HH = fwd_axis(hi, 1); return LL, HL, LH, HH
def inv2(LL, HL, LH, HH):                          # horizontal, then vertical
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

# ---------------------------------------------------------------- windowed decoder
class Windowed:
    """Decodes rectangles of out_s (the iSyntax IDWT output at scale s) from rectangles of the stored coefficients.
    out_s[i] needs L[n], H[n] for n in [(i0>>1)-1, ((i1-1)>>1)+3) (one extra on the left, two on the right absorb the
    window-edge duplication of the lifting); L = cur_s[n] = out_{s+1}[n-1] is requested recursively, up to LLtop.
    At the true array ends the window is clamped to the array, so the lifting duplicates exactly as the full decoder
    does (that duplication is the real boundary of the global transform, not an artefact); cur_s[-1] is the zero row/column."""
    def __init__(self, H, LLtop, top, Hp, Wp):
        self.H, self.LLtop, self.top = H, LLtop, top; self.touched = set()
        self.lim = {s: (need(Hp, s), need(Wp, s)) for s in range(top + 1)}
    def shape(self, s):                            # shape of out_s
        return 2 * self.lim[s][0], 2 * self.lim[s][1]
    def touch(self, s, kind, r0, r1, c0, c1):
        for br in range(r0 // 128, (r1 - 1) // 128 + 1):
            for bc in range(c0 // 128, (c1 - 1) // 128 + 1): self.touched.add((s, kind, bc, br))
    def cur_region(self, s, r0, r1, c0, c1):       # global LL input of step s (= LL_{s+1}), indices may start at 0
        assert 0 <= r0 < r1 <= self.lim[s][0] and 0 <= c0 < c1 <= self.lim[s][1], 'cur_region outside the array'
        if s == self.top:
            self.touch(s, 'LL', r0, r1, c0, c1); a = self.LLtop[r0:r1, c0:c1].astype(np.int32)
            assert a.shape == (r1 - r0, c1 - c0); return a
        pr = 1 if r0 == 0 else 0; pc = 1 if c0 == 0 else 0                         # cur_s[0] = 0 (the shift)
        if r1 - pr <= 0 or c1 - pc <= 0: return np.zeros((r1 - r0, c1 - c0), dtype=np.int32)
        sub = self.out_region(s + 1, r0 - 1 + pr, r1 - 1, c0 - 1 + pc, c1 - 1)
        return np.pad(sub, ((pr, 0), (pc, 0)))
    def out_region(self, s, r0, r1, c0, c1):
        R, C = self.lim[s]; assert 0 <= r0 < r1 <= 2 * R and 0 <= c0 < c1 <= 2 * C, 'out_region outside the array'
        nr0 = max(0, (r0 >> 1) - 1); nr1 = min(R, ((r1 - 1) >> 1) + 3); nc0 = max(0, (c0 >> 1) - 1); nc1 = min(C, ((c1 - 1) >> 1) + 3)
        L = self.cur_region(s, nr0, nr1, nc0, nc1)
        self.touch(s, 'H', nr0, nr1, nc0, nc1)
        HL, LH, HH = [self.H[s][b][nr0:nr1, nc0:nc1] for b in range(3)]
        assert HL.shape == L.shape
        out = lift_axis(lift_axis(L, HL, 1), lift_axis(LH, HH, 1), 0)
        return out[r0 - 2 * nr0:r1 - 2 * nr0, c0 - 2 * nc0:c1 - 2 * nc0]

def rect_fill(arr, getter, W):
    """Fill only the four rings of depth W of arr (all of it if the rings cover it); getter(r0,r1,c0,c1) in local coords."""
    n = arr.shape[0]
    if n <= 2 * W: arr[:, :] = getter(0, n, 0, n); return
    arr[:, :W] = getter(0, n, 0, W); arr[:, n - W:] = getter(0, n, n - W, n)
    arr[:W, :] = getter(0, W, 0, n); arr[n - W:, :] = getter(n - W, n, 0, n)

def build_windowed(dec, k, i0y, i0x, F, D, h, rng):
    """JPEG 2000 coefficients of the F x F frame of DICOM layer k whose first pixel is out_k[i0y, i0x].
    Alignment: (i0 + 3) must be a multiple of 2^D (frame origin must sit on the JPEG 2000 grid; see the .md)."""
    top = dec.top; W = 2 * h + 4; qy = i0y + 1; qx = i0x + 1
    assert (qy + 2) % (1 << D) == 0 and (qx + 2) % (1 << D) == 0, 'frame origin not aligned'
    curf = rng.integers(-1000, 1000, (F, F)).astype(np.int32)                  # garbage everywhere ...
    rect_fill(curf, lambda r0, r1, c0, c1: dec.out_region(k, i0y + r0, i0y + r1, i0x + c0, i0x + c1), W)   # ... except the rings
    lv = []; work = total = 0
    for j in range(D):
        s = k + j; n = curf.shape[0]; nn = n // 2; total += n * n
        if s > top:                                                              # beyond the stored scales: local transform
            LL, HL, LH, HH = fwd2(curf); lv.append((HL, LH, HH)); curf = LL; work += n * n; continue
        hx = ((qx + 2) >> (j + 1)) - 1 + np.arange(nn); hy = ((qy + 2) >> (j + 1)) - 1 + np.arange(nn)
        dec.touch(s, 'H', hy[0] - 1, hy[-1], hx[0], hx[-1] + 1); dec.touch(s, 'H', hy[0], hy[-1] + 1, hx[0] - 1, hx[-1])
        dec.touch(s, 'H', hy[0], hy[-1] + 1, hx[0], hx[-1] + 1)
        HLc = dec.H[s][0][np.ix_(hy - 1, hx)].astype(np.int32); LHc = dec.H[s][1][np.ix_(hy, hx - 1)].astype(np.int32)
        HHc = dec.H[s][2][np.ix_(hy, hx)].astype(np.int32)                      # interior copy straight from the blocks
        by = ((qy + 2) >> (j + 1)) - 2; bx = ((qx + 2) >> (j + 1)) - 2
        LLc = rng.integers(-1000, 1000, (nn, nn)).astype(np.int32)
        if j == D - 1 or s == top or nn <= 2 * W:                               # complete LL needed (output / next level local)
            LLc[:, :] = dec.cur_region(s, by, by + nn, bx, bx + nn)
        else:
            rect_fill(LLc, lambda r0, r1, c0, c1: dec.cur_region(s, by + r0, by + r1, bx + c0, bx + c1), W)
        if nn <= 2 * h:                                                          # ring covers the whole level
            LL, HL, LH, HH = fwd2(curf); LLc, HLc, LHc, HHc = LL.copy(), HL.copy(), LH.copy(), HH.copy(); work += n * n
        else:
            def strip(rows=None, cols=None):
                w = curf
                if cols is not None: w = w[:, cols]
                if rows is not None: w = w[rows, :]
                nonlocal work; work += w.size; return fwd2(w)
            def put(res, rs, cs, tr, tc):
                for dst, src in zip((LLc, HLc, LHc, HHc), res): dst[rs, cs] = src[tr, tc]
            r = strip(cols=slice(0, W));        put(r, slice(0, nn), slice(0, h), slice(0, nn), slice(0, h))
            r = strip(cols=slice(n - W, n));    put(r, slice(0, nn), slice(nn - h, nn), slice(0, nn), slice(2, h + 2))
            r = strip(rows=slice(0, W));        put(r, slice(0, h), slice(0, nn), slice(0, h), slice(0, nn))
            r = strip(rows=slice(n - W, n));    put(r, slice(nn - h, nn), slice(0, nn), slice(2, h + 2), slice(0, nn))
        lv.append((HLc, LHc, HHc)); curf = LLc
    return curf, lv, work / total

# ---------------------------------------------------------------- main
def chunk_of(s, kind, bc, br):                       # chunk = 3 scales: (0,1,2) and (3,4,5); LL lives in the top chunk
    g = s // 3; sh = 3 * g + 2 - s; return (g, bc >> sh, br >> sh)
if __name__ == '__main__':
    path, csvp = sys.argv[1], sys.argv[2]; t0 = time.time()
    g, rows = open_file(path, csvp); Hp, Wp, top, pad = g['Hp'], g['Wp'], g['top'], g['pad']
    size = {}
    for r in rows:
        off = (3 << r['s']) - 2 + ((1 << r['s']) if r['k'] == 0 else 0)
        size[(r['c'], r['s'], r['k'], (r['x'] - g['ox'] - off) // (256 << r['s']), (r['y'] - g['oy'] - off) // (256 << r['s']))] = r['z']
    D, F, HS = 5, 256, (2, 4)
    # tissue centre in level-0 padded coordinates (from tile (60,40) of the earlier tests)
    cx0, cy0 = 15169 + pad + 128, 10049 + pad + 128
    frames = []
    for k in range(top + 1):                         # one tissue frame per DICOM layer
        ix = ((cx0 + 2) >> k) - 3; iy = ((cy0 + 2) >> k) - 3
        frames.append(('layer %d tissue' % k, k, 32 * round((iy - F // 2 + 3) / 32) - 3, 32 * round((ix - F // 2 + 3) / 32) - 3, F))
    frames += [('layer 0 tissue 512', 0, 64 * 157 + pad, 64 * 237 + pad, 512), ('layer 0 tissue edge', 0, 64 * 137 + pad, 64 * 137 + pad, F),
               ('layer 0 blank glass', 0, 64 * 37 + pad, 64 * 37 + pad, F), ('layer 0 scan edge', 0, 64 * 197 + pad, 64 * 385 + pad, F)]
    allok = True; foot = {}
    for c, cname in enumerate(('Y', 'Co', 'Cg')):
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
        cur = LLtop.astype(np.int32); outs = {}; curs = {}
        for s in range(top, -1, -1):                 # FULL decode: reference only (ground truth for the checks)
            n_r, n_c = need(Hp, s), need(Wp, s); cur = cur[:n_r, :n_c]; curs[s] = cur
            out = lift_axis(lift_axis(cur, H[s][0][:n_r, :n_c], 1), lift_axis(H[s][1][:n_r, :n_c], H[s][2][:n_r, :n_c], 1), 0)
            outs[s] = out; nxt = np.zeros_like(out); nxt[1:, 1:] = out[:-1, :-1]; cur = nxt
        print('\n=== channel %s (full decode for reference: %.0fs)' % (cname, time.time() - t0), flush=True)
        dec = Windowed(H, LLtop, top, Hp, Wp); rng = np.random.default_rng(5)
        # (1) window decoder alone vs the full decode, random rectangles at every scale
        bad = 0
        for s in range(top + 1):
            for _ in range(20):
                lo = 16 << (top - s)                  # keep the recursion away from the array origin (index halves per scale)
                if outs[s].shape[0] - 400 <= lo or outs[s].shape[1] - 400 <= lo: continue
                h_, w_ = rng.integers(3, 300, 2); r0 = int(rng.integers(lo, outs[s].shape[0] - 400)); c0 = int(rng.integers(lo, outs[s].shape[1] - 400))
                if r0 + h_ > outs[s].shape[0] - 8 or c0 + w_ > outs[s].shape[1] - 8: continue
                bad += int((dec.out_region(s, r0, r0 + h_, c0, c0 + w_) != outs[s][r0:r0 + h_, c0:c0 + w_]).sum())
        print(' windowed IDWT vs full decode on random rectangles (scales 0..%d): %d mismatching samples' % (top, bad))
        if bad: allok = False
        for name, k, i0y, i0x, Fr in frames:
            if i0y < 8 or i0x < 8 or i0y + Fr > outs[k].shape[0] - 8 or i0x + Fr > outs[k].shape[1] - 8:
                print(' %-20s skipped (frame outside layer: origin (%d,%d), layer %s)' % (name, i0y, i0x, outs[k].shape)); continue
            fr = outs[k][i0y:i0y + Fr, i0x:i0x + Fr]; refLL, refLV = full_pyramid(fr, D); line = []
            for h in HS:
                dec.touched = set(); LLb, lvb, wk = build_windowed(dec, k, i0y, i0x, Fr, D, h, rng)
                nbad = int((LLb != refLL).sum()) + sum(int((a1 != a2).sum()) for l1, l2 in zip(lvb, refLV) for a1, a2 in zip(l1, l2))
                pix_ok = bool((inv_pyramid(LLb, lvb) == fr).all()); line.append('h=%d: %d mismatches, pixels exact %s, strip work %.0f%%' % (h, nbad, pix_ok, 100 * wk))
                if h == 4:
                    if nbad or not pix_ok: allok = False
                    foot[(c, name)] = set(dec.touched)
            nloc = sum(1 for j in range(D) if k + j > top)
            print(' %-20s k=%d F=%d (stored scales %d..%d, %d pure-local level(s)) | %s' % (name, k, Fr, k, min(top, k + D - 1), nloc, ' | '.join(line)), flush=True)
    print('\n--- footprint of the windowed approach (blocks and chunks read to build one frame; h=4) ---')
    for name, k, i0y, i0x, Fr in frames:
        if (0, name) not in foot: continue
        same = all(foot.get((c, name)) == foot[(0, name)] for c in (1, 2)); T = foot[(0, name)]
        chunks = {chunk_of(*t) for t in T}; nbytes = sum(size.get((c, t[0], 1 if t[1] == 'H' else 0, t[2], t[3]), 0) for c in range(3) for t in T)
        print(' %-20s %3d codeblocks per colour, %2d chunks, %.2f MB compressed (3 colours); same blocks for all channels: %s' % (name, len(T), len(chunks), nbytes / 1e6, same))
    print('\nALL OK (h=4):', allok, ' total %.0fs' % (time.time() - t0))
