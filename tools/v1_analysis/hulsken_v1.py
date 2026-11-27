# Independent pure-Python port of the Hulsken v1 decoder (read-only analysis; written from isyntax.c:2551 and the
# Philips MATLAB hulskendecompress/decodeserializedcoefficients/bitdeserialize4x2).
import struct, csv, sys, random, collections

def win(data, bitpos, nbytes=8):
    b = bitpos >> 3
    return int.from_bytes(data[b:b+nbytes].ljust(nbytes, b'\0'), 'little') >> (bitpos & 7)

def decode_block_v1(data, ncoef, npx=128*128):
    """returns (masks[list], planes_bytes(bytes), info dict) or raises"""
    total_bits = len(data) * 8
    ser_len = struct.unpack_from('<I', data, 0)[0]
    zr_sym = data[4]; zr_bits = data[5]
    pos = 48
    # --- Huffman tree (pre-order: 0=internal(left,right) 1=leaf+8bit symbol); code bit order LSB-first, left=0,right=1
    codes = {}
    def rd(code, depth):
        nonlocal pos
        if pos >= total_bits: raise ValueError('tree out of bounds')
        bit = win(data, pos, 2) & 1; pos += 1
        if bit:
            sym = win(data, pos, 2) & 0xFF; pos += 8
            codes[sym] = (code, depth)
        else:
            rd(code, depth + 1)
            rd(code | (1 << depth), depth + 1)
    rd(0, 0)
    # decode table keyed by (size, code)
    maxlen = max(d for c, d in codes.values()) if codes else 0
    lut = {}
    for sym, (c, d) in codes.items():
        dd = max(d, 1)                      # single-leaf tree: code length 1
        for hi in range(1 << (maxlen - dd)) if maxlen >= dd else [0]:
            lut[c | (hi << dd)] = (sym, dd) if maxlen >= dd else (sym, dd)
    # build direct table over maxlen bits
    maxl = max(maxlen, 1)
    table = [None] * (1 << maxl)
    for sym, (c, d) in codes.items():
        dd = max(d, 1)
        for hi in range(1 << (maxl - dd)):
            table[c | (hi << dd)] = (sym, dd)
    zr_code, zr_size = codes.get(zr_sym, (0, 0))[0], max(codes.get(zr_sym, (0, 0))[1], 1)   # zero-run symbol may be absent from the tree (never used): isyntax.c leaves code 0 in that case
    zr_mask = (1 << zr_size) - 1
    cmask = (1 << zr_bits) - 1
    out = bytearray()
    mask_l = (1 << maxl) - 1
    while pos < total_bits and len(out) < ser_len:
        w = win(data, pos)
        sym, size = table[w & mask_l]
        pos += size
        if sym == zr_sym:
            n = (w >> size) & cmask; pos += zr_bits
            if n > 0:
                while True:
                    w2 = win(data, pos)
                    if (w2 & zr_mask) == zr_code and pos < total_bits:
                        extra = (w2 >> zr_size) & cmask
                        n = (n << zr_bits) | extra; pos += zr_size + zr_bits
                        if len(out) + n >= ser_len or pos >= total_bits: break
                    else:
                        break
                out.extend(bytes(min(n, ser_len - len(out))))
                if len(out) >= ser_len: break
            else:
                out.append(sym)
        else:
            out.append(sym)
    info = dict(codes=codes, serialized=bytes(out), ser_len=ser_len, decoded=len(out), zr_sym=zr_sym, zr_bits=zr_bits, tree_leaves=len(codes), bits_used=pos, bits_total=total_bits)
    bytes_per_plane = npx // 8
    full = 16 * ncoef * bytes_per_plane
    if len(out) < full:     # v1: masks stored at the end of the serialized data when any plane is empty
        mb = bytes(out[len(out) - 2 * ncoef:])
        masks = [struct.unpack_from('<H', mb, 2 * k)[0] for k in range(ncoef)]
        planes = bytes(out[:len(out) - 2 * ncoef])
    else:
        masks = [0xFFFF] * ncoef; planes = bytes(out)
    nplanes = sum(bin(m).count('1') for m in masks)
    info['nplanes'] = nplanes; info['len_ok'] = (len(planes) == nplanes * bytes_per_plane)
    return masks, planes, info

def unpack_values(masks, planes, ncoef, npx=128*128):
    """v1 order: per coefficient: mask bit0=sign (value bit 15), mask bit i>=1 -> value bit i-1."""
    bpp = npx // 8; vals = [[0] * npx for _ in range(ncoef)]; p = 0
    for k in range(ncoef):
        for i in range(16):
            if not (masks[k] >> i) & 1: continue
            shift = 15 if i == 0 else i - 1
            plane = planes[p * bpp:(p + 1) * bpp]; p += 1
            v = vals[k]
            for j, byte in enumerate(plane):
                if byte:
                    base = j * 8
                    for bit in range(8):
                        if (byte >> bit) & 1: v[base + bit] |= 1 << shift
    return vals

if __name__ == '__main__':
    path, csvp, per_group, seed = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
    rows = list(csv.DictReader(open(csvp)))
    f = open(path, 'rb')
    groups = collections.defaultdict(list)
    sizes = collections.defaultdict(list)
    for r in rows:
        key = (int(r['color']), int(r['scale']), 'LL' if int(r['coef']) == 0 else 'H')
        sizes[key].append(int(r['size']))
        if int(r['size']) > 8: groups[key].append(r)
    random.seed(seed)
    print('color scale band | blocks  empty(<=8B)%  sampled | union of MASK bits (v1: bit0=sign, bit k = magnitude bit k-1) per subband | min nonzero |c| | all multiples of')
    for key in sorted(groups):
        c, s, band = key
        pick = random.sample(groups[key], min(per_group, len(groups[key])))
        ncoef = 1 if band == 'LL' else 3
        union = [0] * ncoef; allmask = [0xFFFF] * ncoef; per_bit_use = [[0] * 16 for _ in range(ncoef)]
        okc = 0
        for r in pick:
            f.seek(int(r['offset'])); data = f.read(int(r['size']))
            masks, planes, info = decode_block_v1(data, ncoef)
            okc += info['len_ok']
            for k in range(ncoef):
                union[k] |= masks[k]; allmask[k] &= masks[k]
                for i in range(16):
                    if (masks[k] >> i) & 1: per_bit_use[k][i] += 1
        n = len(sizes[key]); emp = sum(1 for z in sizes[key] if z <= 8)
        def bits(m): return ''.join('1' if (m >> i) & 1 else '.' for i in range(16))
        names = ['LL'] if band == 'LL' else ['HL', 'LH', 'HH']
        print(f'c{c} s{s} {band:2s} | {n:6d} {100*emp/n:5.1f}% {len(pick):3d} (len ok {okc}) | ' + ' '.join(f'{names[k]}:{bits(union[k])}' for k in range(ncoef)))
