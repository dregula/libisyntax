# Does OpenJPEG's reversible 5/3 forward transform use the same conventions as the NumPy transform used in all other tests here
# (low-pass at even local coordinates, whole-sample mirror at the frame edge, vertical-then-horizontal)?  Indirect test through
# Pillow + OpenJPEG: build an image as the inverse transform of a SPARSE coefficient set; if OpenJPEG's forward transform returns
# exactly that sparse set, the codestream is tiny; with a one-pixel phase error the coefficients become dense and the file large.
# Also checks lossless round trip and that mct=1 applies the reversible colour transform.  Needs Pillow (with OpenJPEG) and numpy.
#   python3 -B -I openjpeg_convention_check.py
import sys, os, io, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from PIL import Image, features
from windowed_idwt_halo import full_pyramid, inv_pyramid
print('Pillow jpg_2000:', features.check('jpg_2000'), 'OpenJPEG', features.version('jpg_2000'))
def enc(arr, **kw):
    im = Image.fromarray(arr.astype(np.uint8)); b = io.BytesIO(); im.save(b, format='JPEG2000', no_jp2=True, irreversible=False, **kw); return b.getvalue()
def dec(data): return np.array(Image.open(io.BytesIO(data)))
rng = np.random.default_rng(4)
for F, D in ((256, 5), (512, 5), (1024, 5)):
    for attempt in range(50):
        n = F >> D; LL = rng.integers(90, 160, (n, n)).astype(np.int32); lv = []; nz = 0
        for j in range(D):
            m = F >> (j + 1); bands = []
            for b in range(3):
                a = np.zeros((m, m), dtype=np.int32); k = max(3, m * m // 4000)
                idx = rng.integers(0, m, (k, 2)); a[idx[:, 0], idx[:, 1]] = rng.integers(-12, 13, k); nz += int((a != 0).sum()); bands.append(a)
            lv.append(tuple(bands))
        x = inv_pyramid(LL, lv)
        if x.min() >= 0 and x.max() <= 255: break
    else: print('no in-range image found'); continue
    # sanity: our own forward gives back exactly the sparse set
    LL2, lv2 = full_pyramid(x, D); ok_own = (LL2 == LL).all() and all((a == b).all() for l1, l2 in zip(lv, lv2) for a, b in zip(l1, l2))
    res = {}
    for name, img in (('matched', x), ('shift x+1', np.roll(x, 1, axis=1)), ('shift y+1', np.roll(x, 1, axis=0)), ('shift x+2', np.roll(x, 2, axis=1))):
        d = enc(img, num_resolutions=D + 1, codeblock_size=(64, 64)); res[name] = len(d)
        assert (dec(d) == img).all(), 'lossless round trip failed'
    print('F=%d D=%d: %d non-zero coefficients, own transform recovers them: %s | OpenJPEG bytes: %s | all round trips lossless' % (
        F, D, nz + LL.size, ok_own, ', '.join('%s %d' % kv for kv in res.items())))
# colour transform: a strongly correlated RGB image (G = R + small, B = R - small) codes much smaller with mct=1 if RCT is applied
r = rng.integers(0, 256, (256, 256)); rgb = np.stack([r, np.clip(r + rng.integers(-2, 3, r.shape), 0, 255), np.clip(r - rng.integers(0, 4, r.shape), 0, 255)], -1)
s0 = len(enc(rgb, num_resolutions=6, mct=0)); s1 = len(enc(rgb, num_resolutions=6, mct=1))
print('RGB correlated test image: mct=0 %d bytes, mct=1 %d bytes; round trip lossless: %s' % (s0, s1, (dec(enc(rgb, num_resolutions=6, mct=1)) == rgb).all()))
