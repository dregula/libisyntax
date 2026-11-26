# Compare libisyntax tile PNGs (a run folder of tools/run_c_decoder.sh) with the Philips SDK's LZW TIFF pyramid levels (read-only).
#   python3 -B -I compare_sdk_tiff_tiles.py <run_folder> <test_materials_dir> [stem=PAS]
# The SDK TIFFs are expected as <dir>/<stem>_level<L>_LZW_cache1G_.tiff (levels 2, 3, 5; the level-1 files are truncated, IFD offset 0).
# For every tiles/L<L>_x<x>_y<y>.png: find the offset (dx, dy) of the tile inside the TIFF level (search +-12 px around the padding
# estimate), report the share of identical RGB pixels and the per-channel difference histogram. The two builds of the C decoder
# (scalar truncating '/2' vs SIMD floor '>>1', see CLAUDE.md) differ in about 64 % of the pixels of a busy tile, so the SDK output matches
# at most one of them in a busy tile; that identifies the SDK's YCoCg -> RGB inverse.
import sys, os, re, glob, numpy as np
from PIL import Image
Image.MAX_IMAGE_PIXELS = None
run, mat = sys.argv[1], sys.argv[2]; stem = sys.argv[3] if len(sys.argv) > 3 else 'PAS'
est = {5: 23, 3: 96, 2: 191, 1: 382}          # padding at level L ~ 765 / 2^L (measured at level 5: 23)
cache = {}
def level(L):
    if L not in cache:
        f = os.path.join(mat, '%s_level%d_LZW_cache1G_.tiff' % (stem, L))
        cache[L] = np.array(Image.open(f).convert('RGB')).astype(np.int16) if os.path.exists(f) else None
    return cache[L]
print('run folder:', run)
for p in sorted(glob.glob(os.path.join(run, 'tiles', 'L*_x*_y*.png'))):
    m = re.match(r'L(\d+)_x(\d+)_y(\d+)\.png', os.path.basename(p)); L, tx, ty = map(int, m.groups())
    T = level(L)
    if T is None: continue
    P = np.array(Image.open(p).convert('RGB')).astype(np.int16); best = None
    for dy in range(-12, 13):
        for dx in range(-12, 13):
            y0 = ty * 256 - est[L] + dy; x0 = tx * 256 - est[L] + dx
            if y0 < 0 or x0 < 0 or y0 + 256 > T.shape[0] or x0 + 256 > T.shape[1]: continue
            w = T[y0:y0 + 256, x0:x0 + 256]; s = float(np.mean((w == P).all(-1)))
            if best is None or s > best[0]: best = (s, dy, dx, w)
    if best is None: print('%-22s outside the TIFF level' % os.path.basename(p)); continue
    s, dy, dx, w = best; d = w - P; hist = ' '.join('%s:%d' % (c, int((d[..., i] != 0).sum())) for i, c in enumerate('RGB'))
    busy = float(np.mean(np.abs(P[..., 0] - P[..., 2]) > 0))
    print('%-22s offset(dy,dx)=(%+d,%+d) identical %6.2f%% | pixels differing per channel %s | max|d| %d | chroma-active px %.0f%%' % (
        os.path.basename(p), dy, dx, 100 * s, hist, np.abs(d).max(), 100 * busy))
