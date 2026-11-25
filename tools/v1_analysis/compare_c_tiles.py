# Compare tiles written by the C decoder (isyntax_example ... <level> <tx> <ty> out.png) with the independent numpy
# reconstruction (reconstruct_np.py with SAVE_PLANES=1). Level-0 tiles only. Needs numpy.
#   <venv>/python3 -I compare_c_tiles.py <run_folder_with_tiles> <planes_dir_with_Y.npy_Co.npy_Cg.npy>
# For each level-0 tile it searches a small window of pixel offsets around the predicted position
# (tile pixel j <-> padded column 256*tx + j + 1, i.e. valid-region column 256*tx + j + 1 - pad) and reports, for the best offset,
# the fraction of pixels whose RGB is identical to the numpy result under the two YCoCg-R roundings:
#   floor (">>1", the SIMD path)  and  trunc ("/2", the scalar path).  The better match tells which path the C build used.
import sys, os, re, zlib, struct, glob
import numpy as np

def read_png(path):
    d = open(path, 'rb').read(); assert d[:8] == b'\x89PNG\r\n\x1a\n'
    pos = 8; idat = b''; w = h = bd = ct = None
    while pos < len(d):
        ln, = struct.unpack('>I', d[pos:pos + 4]); typ = d[pos + 4:pos + 8]; body = d[pos + 8:pos + 8 + ln]; pos += 12 + ln
        if typ == b'IHDR': w, h, bd, ct, _, _, il = struct.unpack('>IIBBBBB', body); assert bd == 8 and il == 0 and ct in (2, 6)
        elif typ == b'IDAT': idat += body
    ch = 4 if ct == 6 else 3; stride = w * ch
    raw = zlib.decompress(idat); out = np.zeros((h, stride), dtype=np.uint8); prev = np.zeros(stride, dtype=np.int32); p = 0
    for y in range(h):
        f = raw[p]; line = np.frombuffer(raw[p + 1:p + 1 + stride], dtype=np.uint8).astype(np.int32); p += 1 + stride
        if f == 0: cur = line
        elif f == 1:
            cur = line.copy()
            for i in range(ch, stride): cur[i] = (cur[i] + cur[i - ch]) & 255
        elif f == 2: cur = (line + prev) & 255
        elif f == 3:
            cur = line.copy()
            for i in range(stride): cur[i] = (cur[i] + (((cur[i - ch] if i >= ch else 0) + prev[i]) >> 1)) & 255
        elif f == 4:
            cur = line.copy()
            for i in range(stride):
                a = cur[i - ch] if i >= ch else 0; b = prev[i]; c = prev[i - ch] if i >= ch else 0
                pa, pb, pc = abs(b - c), abs(a - c), abs(a + b - 2 * c)
                pr = a if (pa <= pb and pa <= pc) else (b if pb <= pc else c)
                cur[i] = (cur[i] + pr) & 255
        else: raise ValueError('bad PNG filter')
        out[y] = cur; prev = cur
    return out.reshape(h, w, ch)

def rgb_variants(Y, Co, Cg):
    Y = np.abs(Y.astype(np.int32)); Co = Co.astype(np.int32); Cg = Cg.astype(np.int32)   # abs(Y) as in isyntax.c:2399
    res = {}
    for name, sh in (('floor', lambda a: a >> 1), ('trunc', lambda a: np.trunc(a / 2).astype(np.int32))):
        tmp = Y - sh(Cg); G = tmp + Cg; B = tmp - sh(Co); R = B + Co
        res[name] = np.stack([np.clip(X, 0, 255) for X in (R, G, B)], axis=-1).astype(np.uint8)
    return res

if __name__ == '__main__':
    run, planes = sys.argv[1], sys.argv[2]
    pad = 189
    Y, Co, Cg = (np.load(os.path.join(planes, n + '.npy'), mmap_mode='r') for n in ('Y', 'Co', 'Cg'))
    H, W = Y.shape
    for f in sorted(glob.glob(os.path.join(run, 'tiles', 'L0_x*_y*.png'))):
        m = re.search(r'L0_x(\d+)_y(\d+)\.png', f); tx, ty = int(m.group(1)), int(m.group(2))
        t = read_png(f); rgb_c = t[..., :3]
        best = None
        for dy in range(-6, 7):
            for dx in range(-6, 7):
                x0 = 256 * tx + 1 - pad + dx; y0 = 256 * ty + 1 - pad + dy
                if x0 < 0 or y0 < 0 or x0 + 256 > W or y0 + 256 > H: continue
                v = rgb_variants(Y[y0:y0 + 256, x0:x0 + 256], Co[y0:y0 + 256, x0:x0 + 256], Cg[y0:y0 + 256, x0:x0 + 256])
                for name, a in v.items():
                    frac = float((a == rgb_c).all(axis=-1).mean())
                    if best is None or frac > best[0]: best = (frac, name, dx, dy, float(np.abs(a.astype(int) - rgb_c.astype(int)).mean()))
        print(f'tile L0 x{tx} y{ty}:', 'no overlap with the valid region' if best is None else
              f'best offset (dx={best[2]},dy={best[3]}) variant={best[1]} identical pixels={100*best[0]:.2f}% mean|diff|={best[4]:.3f}')
