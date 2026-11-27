# Strip-wise (bounded-memory) independent numpy reconstruction of a v1 .isyntax WSI. Read-only; writes only to <outdir>.
# Replaces the whole-image approach of reconstruct_np.py (which keeps every coefficient array and three full int32/int16
# planes in memory: about 30 bytes per padded pixel, 15 GB for the 0.55 Gpx file 1.isyntax, impossible for 7 Gpx files).
# usage: reconstruct_strips.py <file.isyntax> <blocks.csv> <outdir> [downsample=32]
# env:   JOBS=<worker processes, default 12>   BAND_MB=<target MB of one int32 band plane, default 64>
#        VERIFY_PLANES=<dir with Y.npy Co.npy Cg.npy from reconstruct_np.py>   (bit-exact comparison, per band)
# Method (all verified against the whole-image reconstruction in reconstruct_np.py, see VERIFY_PLANES):
#   The 5/3 inverse at one scale is, along one axis,   xL[n] = L[n] - ((H[n]+H[n+1]+2)>>2),  xH[n] = H[n] + ((xL[n-1]+xL[n])>>1),
#   out[2n] = xH[n], out[2n+1] = xL[n]   and the next finer LL is  L'[n] = out[n-1]  (reconstruct_np.py).
#   Applied to a window of coefficient rows [n0,n1) with mirroring at the window ends, the output rows [2n0+1, 2n1-2) are exact
#   (the first and last 1-2 rows see mirrored instead of real neighbours). Going top-down, to get exact out_0 rows [A,B) one needs
#   at scale s the window n0 = floor((A-1)/2), n1 = ceil((B+2)/2) and then out_{s+1} rows [n0-1, n1-1) for the next coarser scale.
#   Windows are clipped to [0, N): at the true image border the mirroring is the same as in the whole-image run.
#   Horizontal lifting always uses the full width (a strip is a full-width band of rows).
#   Decoded 128-row block rows are cached and evicted once the band has moved past them.
import sys, os, csv, re, struct, zlib, time, collections, threading, resource
import numpy as np
import multiprocessing as mp
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hulsken_v1 import decode_block_v1

path, csvp, outdir = sys.argv[1], sys.argv[2], sys.argv[3]
DS = int(sys.argv[4]) if len(sys.argv) > 4 else 32
JOBS = int(os.environ.get('JOBS', '12')); BAND_MB = float(os.environ.get('BAND_MB', '64'))
VERIFY = os.environ.get('VERIFY_PLANES')
os.makedirs(outdir, exist_ok=True)
t0 = time.time()
f = open(path, 'rb'); head = f.read(8 << 20)
while head.find(b'\r\n\x04') < 0:
    more = f.read(8 << 20)
    if not more: raise SystemExit('no EOT')
    head += more
eot = head.find(b'\r\n\x04'); xml = head[:eot].decode('utf-8'); del head
gh = xml[xml.index('UFS_IMAGE_GENERAL_HEADERS'):]
gr = [tuple(map(int, r.split())) for r in re.findall(r'Name="UFS_IMAGE_DIMENSION_RANGE"[^>]*>([^<]*)<', gh)][:5]
Wp = (gr[0][2] - gr[0][0] + gr[0][1]) // gr[0][1]; Hp = (gr[1][2] - gr[1][0] + gr[1][1]) // gr[1][1]
ox, oy = gr[0][0], gr[1][0]; L = (gr[3][2] - gr[3][0] + gr[3][1]) // gr[3][1]; top = L - 1
pad = (3 << L) - 3; W, H = Wp - 2 * pad, Hp - 2 * pad
print(f'padded {Wp}x{Hp} levels {L} origin ({ox},{oy}) pad {pad} valid {W}x{H} ({W*H/1e6:.0f} Mpx)', flush=True)
del xml, gh

def need(n, s): return -(-n // (1 << (s + 1))) + 2
# --- block index: (c, s, kind) -> {(tx, ty): (offset, size)}
IDX = collections.defaultdict(dict); nblk_total = 0
for r in csv.DictReader(open(csvp)):
    c, s, k, z = int(r['color']), int(r['scale']), int(r['coef']), int(r['size'])
    off = (3 << s) - 2 + ((1 << s) if k == 0 else 0)
    tx = (int(r['x']) - ox - off) // (256 << s); ty = (int(r['y']) - oy - off) // (256 << s)
    nblk_total += 1
    if z > 8: IDX[(c, s, 0 if k == 0 else 1)][(tx, ty)] = (int(r['offset']), z)
print(f'{nblk_total} block headers, {sum(len(v) for v in IDX.values())} non-empty', flush=True)

# --- block decode (worker side)
P = np.arange(16384); area = P // 16; within = P % 16
perm = ((area // 32) * 4 + within // 4) * 128 + (area % 32) * 4 + within % 4
def unpack(masks, planes, ncoef):
    arr = np.frombuffer(planes, dtype=np.uint8); res = np.zeros((ncoef, 16384), dtype=np.uint16); p = 0
    for k in range(ncoef):
        for i in range(16):
            if (masks[k] >> i) & 1:
                bits = np.unpackbits(arr[p * 2048:(p + 1) * 2048], bitorder='little').astype(np.uint16); p += 1
                res[k] |= bits << (15 if i == 0 else i - 1)
    out = np.zeros_like(res); out[:, perm] = res
    mag = (out & 0x7FFF).astype(np.int32); sgn = (out >> 15).astype(bool)
    return np.where(sgn, -mag, mag).astype(np.int16).reshape(ncoef, 128, 128)
_fh = None
def work(task):
    """task = (key, [(tx, ty, offset, size), ...]) -> (key, [(tx, ty, int16[ncoef,128,128]), ...])"""
    global _fh
    if _fh is None: _fh = open(path, 'rb')
    key, blks = task; ncoef = 1 if key[2] == 0 else 3; res = []
    for tx, ty, o, z in blks:
        _fh.seek(o); data = _fh.read(z)
        masks, planes, info = decode_block_v1(data, ncoef)
        res.append((tx, ty, unpack(masks, planes, ncoef)))
    return key, res

# --- memory monitor: sum of Pss (proportional set size) over this process and its children
PEAK = dict(pss=0, main_rss=0); STOP = False
def _pss(pid):
    try:
        for line in open(f'/proc/{pid}/smaps_rollup'):
            if line.startswith('Pss:'): return int(line.split()[1])
    except Exception: pass
    return 0
def monitor():
    me = os.getpid()
    while not STOP:
        tot = 0
        for d in os.listdir('/proc'):
            if not d.isdigit(): continue
            try:
                st = open(f'/proc/{d}/stat').read().rsplit(')', 1)[1].split()
                if int(d) == me or int(st[1]) == me: tot += _pss(int(d))
            except Exception: pass
        PEAK['pss'] = max(PEAK['pss'], tot)
        for line in open('/proc/self/status'):
            if line.startswith('VmRSS:'): PEAK['main_rss'] = max(PEAK['main_rss'], int(line.split()[1]))
        time.sleep(0.5)


pool = mp.get_context('fork').Pool(JOBS)   # fork: the sandbox forbids the AF_UNIX socket that the default forkserver needs
threading.Thread(target=monitor, daemon=True).start()
CACHE = {}            # (c, s, kind, ty) -> int16 [ncoef, 128, ncols]
ndecoded = 0
def ensure(rowsneeded):
    """rowsneeded: set of (c, s, kind, ty). decode the missing ones through the pool."""
    global ndecoded
    miss = [k for k in rowsneeded if k not in CACHE]
    tasks = []
    for (c, s, kind, ty) in miss:
        ncoef = 1 if kind == 0 else 3
        CACHE[(c, s, kind, ty)] = np.zeros((ncoef, 128, need(Wp, s) + 128), dtype=np.int16)
        blks = [(tx, ty, o, z) for (tx, ty_), (o, z) in IDX[(c, s, kind)].items() if ty_ == ty] if (c, s, kind) in IDX else []
        for i in range(0, len(blks), 16): tasks.append(((c, s, kind, ty), blks[i:i + 16]))
    for key, res in pool.imap_unordered(work, tasks):
        c, s, kind, ty = key
        for tx, ty_, v in res:
            CACHE[key][:, :, tx * 128:tx * 128 + 128] = v; ndecoded += 1
def window(c, s, kind, n0, n1):
    """coefficient rows [n0,n1) (n1<=need(Hp,s)), cols [0,need(Wp,s)) -> int16 [ncoef, n1-n0, ncols]"""
    ncols = need(Wp, s); ncoef = 1 if kind == 0 else 3; out = np.zeros((ncoef, n1 - n0, ncols), dtype=np.int16)
    for ty in range(n0 // 128, (n1 - 1) // 128 + 1):
        a, b = max(n0, ty * 128), min(n1, ty * 128 + 128)
        out[:, a - n0:b - n0] = CACHE[(c, s, kind, ty)][:, a - ty * 128:b - ty * 128, :ncols]
    return out

stats = dict(maxabs=0)
def lift(Lc, Hc, axis):
    Lc = np.moveaxis(Lc, axis, -1).astype(np.int32); Hc = np.moveaxis(Hc, axis, -1).astype(np.int32)
    N = Lc.shape[-1]
    Hn = np.concatenate([Hc[..., 1:], Hc[..., -1:]], axis=-1)
    xL = Lc - ((Hc + Hn + 2) >> 2)
    xLp = np.concatenate([xL[..., :1], xL[..., :-1]], axis=-1)
    xH = Hc + ((xLp + xL) >> 1)
    out = np.empty(Lc.shape[:-1] + (2 * N,), dtype=np.int32); out[..., 0::2] = xH; out[..., 1::2] = xL
    stats['maxabs'] = max(stats['maxabs'], int(np.abs(out).max()))
    return np.moveaxis(out, -1, axis)

def plan(A, B):
    """windows[s] = (n0, n1) coefficient-row window at scale s, for exact out_0 rows [A,B)"""
    wins = {}
    for s in range(0, L):
        Nr = need(Hp, s)
        n0 = max(0, (A - 1) // 2); n1 = min(Nr, -(-(B + 2) // 2))
        wins[s] = (n0, n1); A, B = max(0, n0 - 1), n1 - 1
    return wins

def band_planes(A, B):
    """Y, Co, Cg int32 arrays of padded rows [A,B) (all padded columns): pixel[p] = out_0[p-1]."""
    wins = plan(A - 1, B - 1)
    want = set()
    for c in range(3):
        for s in range(L):
            n0, n1 = wins[s]
            for ty in range(n0 // 128, (n1 - 1) // 128 + 1): want.add((c, s, 1, ty))
        n0, n1 = wins[top]
        for ty in range(n0 // 128, (n1 - 1) // 128 + 1): want.add((c, top, 0, ty))
    ensure(want)
    res = []
    for c in range(3):
        n0, n1 = wins[top]; cur = window(c, top, 0, n0, n1)[0].astype(np.int32)
        for s in range(top, -1, -1):
            n0, n1 = wins[s]; nc = need(Wp, s)
            Hw = window(c, s, 1, n0, n1); cur = cur[:, :nc]
            lowrows = lift(cur, Hw[0], axis=1); highrows = lift(Hw[1], Hw[2], axis=1)
            out = lift(lowrows, highrows, axis=0)                  # abs rows [2*n0, 2*n1)
            del lowrows, highrows, Hw
            if s == 0: break
            m0, m1 = wins[s - 1]; ncn = need(Wp, s - 1)
            nxt = np.zeros((m1 - m0, ncn), dtype=np.int32)
            a = max(m0, 1)                                           # L'[n] = out[n-1], L'[0] = 0
            if m1 > a:
                src = out[a - 1 - 2 * n0:m1 - 1 - 2 * n0, :ncn - 1]
                assert src.shape[0] == m1 - a, (src.shape, m1 - a)
                nxt[a - m0:, 1:1 + src.shape[1]] = src
            cur = nxt; del out
        # out = abs rows [2*n0, 2*n1) of out_0 ; pixel row p = out row p-1
        o0 = 2 * wins[0][0]
        pix = np.zeros((B - A, Wp), dtype=np.int32)
        a_out = A - 1; assert a_out >= o0 and B - 1 <= 2 * wins[0][1], (A, B, o0, wins[0])
        pix[:, 1:] = out[a_out - o0:a_out - o0 + (B - A), :Wp - 1]
        res.append(pix); del out, cur
    return res

def png(p, arr):
    hh, ww, _ = arr.shape; raw = b''.join(b'\x00' + arr[i].tobytes() for i in range(hh))
    def chunk(t, d): return struct.pack('>I', len(d)) + t + d + struct.pack('>I', zlib.crc32(t + d) & 0xffffffff)
    open(p, 'wb').write(b'\x89PNG\r\n\x1a\n' + chunk(b'IHDR', struct.pack('>IIBBBBB', ww, hh, 8, 2, 0, 0, 0)) + chunk(b'IDAT', zlib.compress(raw, 6)) + chunk(b'IEND', b''))
def rgb(Yv, Cov, Cgv, mode):
    sh = (lambda a: a >> 1) if mode == 'floor' else (lambda a: np.trunc(a / 2).astype(np.int32))
    tmp = Yv - sh(Cgv); G = tmp + Cgv; B = tmp - sh(Cov); R = B + Cov
    return R, G, B

if VERIFY: VP = [np.load(os.path.join(VERIFY, n + '.npy'), mmap_mode='r') for n in ('Y', 'Co', 'Cg')]
band = max(64, min(H, int(BAND_MB * 1e6 / 4 / Wp) // 16 * 16))
print(f'band {band} rows (target {BAND_MB} MB per int32 plane), {JOBS} workers', flush=True)
acc = collections.Counter(); sums = np.zeros(3); dmin = np.array([9, 9, 9]); dmax = np.array([-9, -9, -9])
Ymin = Cmin = 10**9; Ymax = Cmax = -10**9; Comin = Cgmin = 10**9; Comax = Cgmax = -10**9
prev_rows = []; ph = (H // DS) * DS; pw = (W // DS) * DS; pend = 0; accum = None
seam_c_sum = np.zeros(256); seam_c_cnt = np.zeros(256); seam_r_sum = np.zeros(256); seam_r_cnt = np.zeros(256)
prev_last = None; verify_bad = 0; verify_n = 0; tb = time.time()
for y0 in range(0, H, band):
    y1 = min(H, y0 + band)
    Yp, Cop, Cgp = band_planes(y0 + pad, y1 + pad)
    Yb, Cob, Cgb = Yp[:, pad:Wp - pad], Cop[:, pad:Wp - pad], Cgp[:, pad:Wp - pad]
    assert -32768 <= min(int(Yb.min()), int(Cob.min()), int(Cgb.min())) and max(int(Yb.max()), int(Cob.max()), int(Cgb.max())) <= 32767
    if VERIFY:
        for pl, ref in zip((Yb, Cob, Cgb), VP):
            verify_n += pl.size; verify_bad += int((pl != ref[y0:y1].astype(np.int32)).sum())
    acc['pixels'] += Yb.size; acc['Y_negative'] += int((Yb < 0).sum())
    Ymin = min(Ymin, int(Yb.min())); Ymax = max(Ymax, int(Yb.max()))
    Comin = min(Comin, int(Cob.min())); Comax = max(Comax, int(Cob.max())); Cgmin = min(Cgmin, int(Cgb.min())); Cgmax = max(Cgmax, int(Cgb.max()))
    Ya = np.abs(Yb)
    Rf, Gf, Bf = rgb(Ya, Cob, Cgb, 'floor'); Rt, Gt, Bt = rgb(Ya, Cob, Cgb, 'trunc'); Rn, Gn, Bn = rgb(Yb, Cob, Cgb, 'floor')
    dd = [Rf - Rt, Gf - Gt, Bf - Bt]
    acc['floor_vs_trunc_differ'] += int(((dd[0] != 0) | (dd[1] != 0) | (dd[2] != 0)).sum())
    for k in range(3): sums[k] += np.abs(dd[k]).sum(); dmin[k] = min(dmin[k], dd[k].min()); dmax[k] = max(dmax[k], dd[k].max())
    acc['rgb_below0'] += int(((Rf < 0) | (Gf < 0) | (Bf < 0)).sum()); acc['rgb_above255'] += int(((Rf > 255) | (Gf > 255) | (Bf > 255)).sum())
    acc['abs_Y_changes_rgb'] += int(((Rn != Rf) | (Gn != Gf) | (Bn != Bf)).sum())
    acc['neg_odd_Cg'] += int(((Cgb < 0) & ((Cgb & 1) == 1)).sum()); acc['neg_odd_Co'] += int(((Cob < 0) & ((Cob & 1) == 1)).sum())
    acc['Cg_neg'] += int((Cgb < 0).sum()); acc['Co_neg'] += int((Cob < 0).sum())
    # seam profile (as reconstruct_np.py): columns: |dY| along x on every 16th row ; rows: |dY| along y on every 16th column
    rsel = [r for r in range(y0, y1) if r % 16 == 0]
    if rsel:
        d = np.abs(np.diff(Yb[[r - y0 for r in rsel]], axis=1)); idx = (np.arange(W - 1) + 1 + pad) % 256
        np.add.at(seam_c_sum, idx, d.sum(axis=0)); np.add.at(seam_c_cnt, idx, d.shape[0])
    Yc = Yb[:, ::16]
    if prev_last is not None:
        dr = np.abs(Yc[0] - prev_last); idr = (y0 - 1 + 1 + pad) % 256; seam_r_sum[idr] += dr.sum(); seam_r_cnt[idr] += dr.size
    if Yc.shape[0] > 1:
        dr = np.abs(np.diff(Yc, axis=0)); idr = (np.arange(y0, y1 - 1) + 1 + pad) % 256
        np.add.at(seam_r_sum, idr, dr.sum(axis=1)); np.add.at(seam_r_cnt, idr, dr.shape[1])
    prev_last = Yc[-1].copy()
    # preview: box average over DS x DS cells (accumulate rows, emit complete DS-row groups)
    rg = np.stack([np.clip(X, 0, 255) for X in (Rf, Gf, Bf)], axis=-1).astype(np.float32)[:, :pw]
    rows_here = rg.shape[0]; i = 0
    if accum is None: accum = np.zeros((pw // DS, 3), dtype=np.float64); acc_rows = 0
    while i < rows_here:
        gy = y0 + i; take = min(rows_here - i, DS - gy % DS)
        accum += rg[i:i + take].reshape(take, pw // DS, DS, 3).sum(axis=(0, 2)); acc_rows += take; i += take
        if (gy + take) % DS == 0:
            if gy + take <= ph: prev_rows.append((accum / (DS * DS)).astype(np.float32))
            accum = np.zeros((pw // DS, 3), dtype=np.float64)
    el = time.time() - tb
    print(f'  band rows {y0}-{y1} of {H}: {el:.0f}s, decoded {ndecoded} blocks, pss peak {PEAK["pss"]/1e6:.2f} GB', flush=True)
    # evict cached block rows that lie entirely above the next band's window
    nxt0 = y1 + pad
    if nxt0 < Hp:
        nw = plan(nxt0 - 1, min(Hp, nxt0 + band) - 1)
        for k in [k for k in CACHE if (k[3] + 1) * 128 <= nw[k[1]][0]]: del CACHE[k]
    del Yp, Cop, Cgp, Yb, Cob, Cgb, Ya, Rf, Gf, Bf, Rt, Gt, Bt, Rn, Gn, Bn, dd, rg
STOP = True
print('FINDINGS (valid region):')
print(f"  pixels {acc['pixels']}; Y min {Ymin} max {Ymax}; Co range [{Comin},{Comax}]; Cg range [{Cgmin},{Cgmax}]")
print(f"  negative Y after IDWT: {acc['Y_negative']}  (abs(Y) changes RGB of {acc['abs_Y_changes_rgb']} pixels)")
print(f"  RGB before clamping: below 0: {acc['rgb_below0']}, above 255: {acc['rgb_above255']} pixels")
print(f"  floor (>>1, SIMD path) vs truncate (/2, scalar path): {acc['floor_vs_trunc_differ']} pixels differ = {100*acc['floor_vs_trunc_differ']/acc['pixels']:.1f}%")
print('   difference floor-minus-trunc range per channel R/G/B:', list(zip(dmin.tolist(), dmax.tolist())), '| mean |diff|:', [round(v / acc['pixels'], 3) for v in sums])
print(f"  fraction with negative Cg {acc['Cg_neg']/acc['pixels']:.3f} (odd {acc['neg_odd_Cg']/acc['pixels']:.3f}), negative Co {acc['Co_neg']/acc['pixels']:.3f} (odd {acc['neg_odd_Co']/acc['pixels']:.3f})")
for name, sm, ct in (('columns', seam_c_sum, seam_c_cnt), ('rows', seam_r_sum, seam_r_cnt)):
    ok = ct > 0; prof = sm[ok] / ct[ok]; resid = np.nonzero(ok)[0]
    print(f'  seam profile ({name}): mean |dY| by position mod 256: min {prof.min():.4f} median {np.median(prof):.4f} max {prof.max():.4f}; max/median {prof.max()/np.median(prof):.3f} (argmax residue {int(resid[prof.argmax()])})')
if VERIFY: print(f'  VERIFY against {VERIFY}: {verify_bad} of {verify_n} samples differ')
print('  max |intermediate| over all levels:', stats['maxabs'], '(int16 limit 32767)')
if prev_rows:
    prev = np.clip(np.stack(prev_rows), 0, 255).round().astype(np.uint8); png(os.path.join(outdir, f'recon_ds{DS}.png'), prev); print('wrote preview', prev.shape)
tt = time.time() - t0
print(f'MEMORY peak: whole process tree Pss {PEAK["pss"]/1e6:.2f} GB; main process RSS {PEAK["main_rss"]/1e6:.2f} GB; per-pixel (padded) {PEAK["pss"]*1024/(Wp*Hp):.2f} B/px (units: kB from /proc, so GB = kB/1e6 approx)')
print(f'blocks decoded {ndecoded}; total {tt:.0f}s')
