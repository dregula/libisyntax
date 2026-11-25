# Philips compression code samples (MATLAB/Octave): review

Source folder: `/mnt/d/Whole_Slide/iSyntax_Reference_Implementation/CompressionCodeSamples/` (35 files, 2,659 lines, Copyright 2020 Koninklijke Philips N.V.). Companion to [CLAUDE.md](CLAUDE.md), [Example_FIC_Files_Claude.md](Example_FIC_Files_Claude.md) and [Third_Party_Documentation.md](Third_Party_Documentation.md).

**Evidence labels:** **[M]** = read in the Philips `.m` files; **[code]** = checked in this repo's C (`file:line`); **[FIC]** = checked in the sample headers; **[doc]** = Philips PDFs; **[derived]** = I worked it out by hand or by a small script; **[inferred]** = my deduction, not tested.

## 0. Scope and limits

- **Nothing was executed.** No Octave or MATLAB is installed. Everything is from reading the code, plus one hand-derived check (section 5.3) that I ran as a short Python script.
- **This is the encoder as well as the decoder.** Unlike libisyntax, the folder contains the whole compression path. That is what lets you answer goal 4 (how Quantizer and Deadzone reach the bitplanes) from Philips' own code rather than inference.
- **Container caveat.** `compressimage.m` / `decompressimage.m` write a simple test container, **not** the iSyntax file format (section 7). The codec pieces (`hulsken*`, wavelet, quantization) are the real ones.
- **Sibling folders:** `ImageFormatCodeSamples/` is now reviewed in [Image_Format_Code_Samples_Claude.md](Image_Format_Code_Samples_Claude.md) (it uses `hulskendecompress`, `wavelet2dilift`, `liftwave` and `ycocg2rgb` from this folder). `ViewingCodeSamples/` is reviewed in [Viewing_Code_Samples_Claude.md](Viewing_Code_Samples_Claude.md) (the two colour files there are the ones discussed in section 6).
- **Licence.** Philips permits redistribution of the source with its copyright notice and disclaimer. I quote only short lines here.

## 1. Main findings

1. **Quantizer and Deadzone are fully specified, and they are encoder-side in-place operations.** `quantanddeadzone.m` zeroes small coefficients and truncates the rest to a multiple of `2^Q`. Nothing is scaled, so the decoder needs no dequantization. This **confirms the prediction** in the earlier documents (section 3).
2. **The spec's deadzone wording is inverted.** The code zeroes coefficients with `abs(c) <= deadzone`. The spec says "greater or equal". The spec is wrong, and the "likely typo" flag in the earlier documents is now resolved.
3. **`Quality=2` in the derivation string is the "Q2" parameter set**, and the (Quantizer, Deadzone) values in every sample FIC are exactly that set (section 3.2).
4. **The reference code does not use symmetric extension.** It zero-pads the image, keeps only "valid" wavelet coefficients, and drops the rest. Our C constants `(3<<s) − 2` and `3·(2^L − 1)` fall out of that logic exactly (section 5.3).
5. **The reference integer 5/3 rounding is the JPEG 2000 reversible 5/3**, which is also what the C code implements (section 5.2).
6. **Philips ships two different colour transforms.** The compression samples use `fix()` (truncate toward zero). The viewing samples use `bitshift(·,-1)`. The two C paths in libisyntax each match one of them (section 6).
7. **The "snake" order is now fully defined** (section 4.1), and the C code's 4×4 reshuffle matches it.

## 2. File inventory

| Group | Files | Role |
|---|---|---|
| Top level | `compressimage.m` (157 lines), `decompressimage.m` (131) | Test container: header, block loop, wavelet, quantization, colour |
| **Hulsken codec** | `hulskencompress.m`, `hulskendecompress.m` | Compress or decompress one coefficient block |
| Bit-plane layer | `bitserialize4x2.m`, `bitdeserialize4x2.m`, `snakeindices4x2.m` | Coefficients ↔ bitplane bytes; the 4×2 snake order |
| Entropy layer | `determinecompressionparameters.m`, `encodeserializedcoefficients.m`, `decodeserializedcoefficients.m` | Zero-run symbol, run-length coding, Huffman coding |
| Huffman | `create8bithuffmantree.m`, `createlutfor8bithuffmantree.m`, `serialize8bithuffmantree.m`, `deserialize8bithuffmantree.m` | Tree build, code table, tree (de)serialization |
| Header | `serializecompressionparameters.m`, `deserializecompressionparameters.m` | Per-block header bits |
| Bit helpers | `dec2bin_lsbfirst.m`, `bin2dec_lsbfirst.m`, `parsevarargin.m` | LSB-first bit strings; argument parsing |
| **Quantization** | `quantanddeadzone.m` | Deadzone and in-place quantization |
| Wavelet | `liftwave.m`, `waveletlift.m`, `waveletilift.m`, `wavelet2dlift.m`, `wavelet2dilift.m`, `waveletzeropadding.m`, `validwaveletcoefficientposition.m`, `waveletfilterindices.m`, `waveletfilterintervals.m` | The lifting DWT and its geometry |
| Wavelet (MATLAB) | `matlablwt2.m`, `matlabilwt2.m` | Alternative using MATLAB's `lwt2`/`ilwt2` (option `matlabdwt=1`) |
| Colour | `rgb2ycocg.m`, `ycocg2rgb.m` | YCoCg-R |
| Index | `files.txt` | Just the list of the files above |

**Call graph (compress):** `compressimage` → `rgb2ycocg` → `wavelet2dlift` (→ `waveletlift`, `waveletzeropadding`) → `quantanddeadzone` → per block `hulskencompress` → `bitserialize4x2`, `determinecompressionparameters`, `encodeserializedcoefficients`, `serializecompressionparameters`. **Decompress** mirrors it: `decompressimage` → `hulskendecompress` → `deserializecompressionparameters`, `decodeserializedcoefficients`, `bitdeserialize4x2` → `wavelet2dilift` → `ycocg2rgb`.

## 3. Quantizer and deadzone (goal 4)

### 3.1 What `quantanddeadzone.m` does **[M]**

The file's own comment: "deadzone and quant exactly as current implementation". For each colour channel `c`, level `lvl` and each of the three high-pass subbands (HL, LH, HH):

1. **Deadzone first:** `abs(coef) <= deadZoneV(lvl, col)` → `0`.
2. **Then quantization:** `coef = fix(coef / 2^q) * 2^q`, with `q = quantV(lvl, col)`. `fix` is truncation toward zero.

where `col = min([2 colorInd])`, so:

- **Column 1 is luma (Y). Column 2 is both chroma channels (Co and Cg).**
- **Rows are levels, with row 1 the first (finest) DWT level.** `HL{1}` is the first decomposition of the full image (see `wavelet2dlift`).
- **LL is never touched** (`oLL = LL`).
- **Levels beyond the number of rows in the parameter matrices are copied unchanged**, so coarser subbands stay lossless.

The comment in the file claims "rounding half towards zero" ("0.5 → 1, −1.5 → −1"), which does not match `fix`. The code is authoritative: `fix` truncates, so |c| is rounded down to a multiple of `2^q`.

### 3.2 The parameter sets, and the match to the sample headers

The file header gives the sets for 8-bit UFS data:

```
Q1 -> quantV=[0 4 ; 0 2], deadZoneV=[0 16 ; 0 8]
Q2 -> quantV=[1 4 ; 0 2], deadZoneV=[2 16 ; 0 8]
```

Rows are levels, columns are [Y, chroma]. For **Q2**, read column by column:

| | Level 0 (finest) | Level 1 |
|---|---|---|
| **Y** (Quantizer, Deadzone) | (1, 2) | (0, 0) |
| **Co, Cg** | (4, 16) | (2, 8) |

The sample FICs carry exactly these values **[FIC]**: group A `(1,2) (0,0)` = Y, group B `(4,16) (2,8) (0,0)` = chroma. The derivation string says `Quality=2` **[FIC]**. So:

- **`Quality=2` = the Q2 set**, which is how the header's quantizer values were chosen.
- The header lists **one more level for chroma than the matrices have** (a third `(0,0)` entry). The matrices have two rows. I don't know why. It makes no difference to the data (zero Q and zero deadzone mean no change).
- Level numbering: header "level 0" = the first DWT level = **scale 0 in the C code**. This answers the earlier open question; the block-size evidence in Example_FIC_Files_Claude.md §6.4 already pointed the same way.

### 3.3 What it does to the bitplanes (the answer to goal 4)

For a subband coefficient `c` after quantization, in the sign-magnitude form the Hulsken coder uses (sign plane + 15 magnitude planes):

- **Magnitude bits `0 … Q−1` are always zero.** |c| is a multiple of `2^Q`, so the lowest `Q` magnitude planes are empty and are **not stored** (the corresponding mask bits are clear).
- **Every nonzero |c| exceeds the deadzone**, so small values never appear. All the high planes that would hold values ≤ deadzone are sparse or empty.
- **A block whose coefficients are all ≤ deadzone becomes all zero**: no bitplanes at all (the 8-byte empty blocks seen in the cluster table).
- The decoder just reads the planes back. Nothing needs rescaling because the encoder never divided anything.

Predicted mask content per setting. The numbering is the compressor-v1 mask: bit 0 = sign, bit `k` = magnitude bit `k−1`. All six samples are compressor v1.

| Subband group | (Q, DZ) | Mask bits forced clear | Smallest nonzero |c| |
|---|---|---|---|
| Y, level 0 | (1, 2) | bit 1 (magnitude bit 0) | 2 (from 3) |
| Y, level 1 and above | (0, 0) | none | 1 |
| Co/Cg, level 0 | (4, 16) | bits 1–4 | 16 (from 17) |
| Co/Cg, level 1 | (2, 8) | bits 1–2 | 8 (from 9) |
| Co/Cg, level 2 and above | (0, 0) | none | 1 |
| LL (all colours) | none | none | n/a |

Status: derived from the algorithm **[M]**, and **now verified on real codeblocks**: in `example_isyntax_files/1.isyntax` all 33,867 blocks decode and none violates the table (see [Example_iSyntax_File_Claude.md](Example_iSyntax_File_Claude.md) §3). The smallest nonzero coefficient per group equals the deadzone (2, 16, 8) and every value is a multiple of `2^Q`.

Two further observations:

- **The deadzone and quantizer are not a single formula.** Deadzone = `2^Q` for (1,2) and (4,16), but 8 = `2^3` for Q = 2. They are tuned values in a table, as the earlier "open question" suspected.
- **The deadzone is applied before the quantizer, on the original coefficient.** So with Q = 4 and DZ = 16, a coefficient of 17 becomes 16 (not 0), and 15 becomes 0.

### 3.4 Corrections to the earlier documents

| Earlier statement | Now |
|---|---|
| Spec: deadzone zeroes values "greater or equal" to the deadzone (flagged as a probable typo) | **Confirmed wrong.** `abs(c) <= dz → 0` |
| Quantizer = "rounded to 2^quantizer" (spec) | Truncated toward zero to a multiple of `2^Q`, not rounded to nearest |
| "Level 0 = scale 0 is only suggested by evidence" | **Confirmed** by the code (`HL{1}` = first DWT level) |
| First PerColor = Y, second = Co and Cg | **Confirmed** (`min([2 colorInd])`) |
| The AI note "Effect of iSyntax Q2-lossy" models `ĉ = sign(q)·⌊(|q|+r)Δ⌉` with a reconstruction offset `r` | **Wrong for iSyntax.** There are no quantization indices and no `r`. The coefficients keep their scale. |
| Predicted mask bits `1…Q` clear | Now **derived from Philips' encoder** (not yet checked on data) |

## 4. Hulsken codec, bit by bit

### 4.1 Bitplane serialization (`bitserialize4x2.m`)

- **Input** is an `[rows, cols, k]` array, where `k = 1` (LL) or `k = 3` (HL, LH, HH in that order **[M]**: the compressor stacks `HL, LH, HH` along dimension 4, then `permute`s them to the third dimension). Rows must be a multiple of 4, columns a multiple of 8. The compressor pads each level to those multiples (`zeropadto8x4multiple`).
- **Sign-magnitude:** the sign plane is `coef < 0`, and the magnitude planes are the bits of `abs(coef)`.
- **A bitplane that is entirely zero is not serialized.** Its bit in the per-coefficient `bitmask` stays clear. Mask bit 0 is the sign plane.
- **Each stored plane is `rows·cols/8` bytes**, one bit per coefficient.
- **Byte layout (the "snake")**, from the file's own comment **[M]**:
  ```
  byte 0 = [x0y0 : x3y0 ; x0y1 : x3y1]     (4 wide x 2 tall, LSB = first)
  byte 1 = [x0y2 : x3y2 ; x0y3 : x3y3]
  byte 2 = [x4y0 : x7y0 ; x4y1 : x7y1]
  byte 3 = [x4y2 : x7y2 ; x4y3 : x7y3]
  ... repeated column first, row second over the plane
  ```
  So two bytes make a 4×4 area (rows 0–1, then rows 2–3), and 4×4 areas run left to right, then down. The C decoder's 16-value reshuffle into 4×4 areas (`isyntax.c:3014-3036`) implements exactly this **[code]**. The "snake" is therefore not a reversing path. It is the order of the 4×2 cells.
- **Plane order differs by compressor version:**

  | | Order of stored planes | Mask bit ↔ plane |
  |---|---|---|
  | **v1** (method 16) | per coefficient: sign, then magnitude bit 0 (LSB) … bit 14 (MSB); coefficient 1, then 2, then 3 | bit 0 = sign, bit `k` = magnitude bit `k−1` |
  | **v2** (method 19) | per bit position, all coefficients together: sign (coef 1,2,3), then magnitude bit 14 (MSB) … 0 | bit 0 = sign, bit `j` = magnitude bit `15−j` |

  Both match the C code **[code]** (`isyntax.c:2928-2977`).
- **v1 only:** if any plane is empty, the bitmasks are **appended after the planes** as `ceil(bits/8)` bytes per coefficient (2 bytes for 16 bits), LSB first. If no plane is empty, no masks are stored and the decoder assumes all ones. The C code does the same **[code]** (`isyntax.c:2889-2906`).
- **Value range:** the coefficient width is a parameter (`bits`, default 16). Magnitudes use `bits−1` planes. iSyntax uses 16.

### 4.2 Zero-run and Huffman parameters (`determinecompressionparameters.m`)

- **Symbols** are the serialized bytes (0–255).
- **Zero-run symbol** = the **least frequent nonzero byte** value in the block (`min(frequencies(2:end))`). If that value actually occurs in the data, it must be escaped (`escapezerorunsymbol = true`).
- **What counts as a run:** two or more consecutive zero bytes. Runs are measured in bytes, and a single zero is coded as an ordinary symbol.
- **Counter width (`zerorunbits`)** is chosen to minimise total size, from a histogram of how many bits each run length needs. Long runs use several consecutive counter blocks.
- **`minimumzerorunlength`** is only used by v2: the run length at which a run symbol plus counter is shorter than repeated Huffman zeros (at least 2). For v1 it is 0, so every run of ≥ 2 zeros is a run.
- **Counters:** v1 stores the run length, **v2 stores length − 1**. The C code agrees **[code]** (`isyntax.c:2817`).
- **The bitplane seektable** (v2 only) holds the offsets of the 2nd, 3rd, … stored bit-position groups. Offsets are in the **encoded bit stream**, measured from the start of the offset table (`bitoffset` includes the table and the serialized tree). The bit width per offset is `floor(log2(serializedlength)) + 5`. The number of entries is (number of bit positions set in any mask) − 1. All of this agrees with the C reader **[code]** (`isyntax.c:2628-2641`), **but the C code never uses the offsets**.
- **Seektable caveat, verified by reading [M]:** the encoder tests for a plane boundary only at the top of its loop (`i > offset`). A zero run that straddles a boundary is written as one run, so the recorded offset points **after** the true start of that plane, with no information about how many zeros belong to the previous plane. The offsets are therefore not safe for plane-parallel decoding. (An earlier AI note made this claim from the encoder files; reading the code confirms it.)

### 4.3 Encoding (`encodeserializedcoefficients.m`)

For each position in the serialized stream:

- **Two or more zeros in a row** start a run. The run is the maximal stretch of zeros.
- **If the next non-zero byte is the zero-run symbol** (and escaping is on), the run is shortened by one and the last zero is coded separately. This makes sure a run can never be followed directly by an unescaped symbol, so the decoder never mistakes it.
- **A run shorter than `minimumzerorunlength`** is written as repeated Huffman zeros.
- **Otherwise the run is coded as** `[Huffman(zerorun symbol)] [counter]`. If the counter needs more than `zerorunbits` bits, it is split into several `[symbol][counter block]` pairs, **most significant block first** (so no block after the first can be zero, which is the escape).
- **A literal zero-run symbol** is escaped as `[Huffman(symbol)] [zerorunbits zeros]`, i.e. a counter of 0.
- **The final symbol** is handled after the loop.

The C decoder implements the inverse: a counter of 0 means a literal symbol, further zero-run symbols extend the counter MSB-first, and the run length is `count` (v1) or `count + 1` (v2) **[code]** (`isyntax.c:2811-2857`).

### 4.4 Huffman tree

- **Construction** (`create8bithuffmantree.m`): symbols sorted by **descending frequency using MATLAB's stable sort**, then repeatedly merging the two least frequent. The file notes that a different sort gives different output with the same ratio. Ties therefore do not matter to a decoder, because the tree is stored.
- **Serialization** (`serialize8bithuffmantree.m`): pre-order traversal. An internal node is the bit `0` followed by its left then right subtree. A leaf is the bit `1` followed by the 8-bit symbol (LSB first). Maximum size 511 + 256·8 bits.
- **Code assignment:** left child = bit 0, right child = bit 1. The **first** decision from the root is the **least significant** bit of the code, and codes are read from the stream LSB first.
- **A tree with one leaf** gets a 1-bit code `0`. The C decoder handles this ("code_size == 0 → 1") **[code]**.
- **Bit order everywhere:** the bit string is LSB-first per byte (`dec2bin_lsbfirst`), and the final stream is zero-padded to whole bytes.

### 4.5 Per-block header

```
v1:  [serialized length : 32 bits]
v2:  [bitmask per coefficient : 16 bits each]                    (1 or 3 masks)
     -- if serialized length > 0 --
     [zero-run symbol : 8] [zero-run counter bits : 8]
v2:  [bitplane offsets : (N-1) x (floor(log2(len))+5) bits]
     [serialized Huffman tree]
     [Huffman / zero-run coded data]  ... padded to a whole byte
```

The fields are **bit-packed, not byte-aligned**, after the first two bytes. The C reader matches **[code]** (`isyntax.c:2590-2731`). If the serialized length is 0 (an all-zero block), nothing follows the length.

### 4.6 M versus C (what agrees, what differs)

| Item | M reference | libisyntax C | Status |
|---|---|---|---|
| v1/v2 header fields and widths | as above | `isyntax.c:2590-2643` | agree |
| Plane order and mask meaning, v1 and v2 | `bitserialize4x2.m` | `isyntax.c:2928-2977` | agree |
| 4×2 / 4×4 cell order | `snakeindices4x2.m` | `isyntax.c:3014-3036` | agree |
| Huffman tree format, code bit order | `*8bithuffmantree.m` | `isyntax.c:2650-2730` | agree |
| Zero-run, escape, v2 count − 1, multi-block counter | `encode…`, `decode…` | `isyntax.c:2811-2857` | agree |
| Bitplane seektable | written, bit offsets | read, **unused** | C ignores it |
| Sign/magnitude to two's complement | none (works in sign-magnitude, signed values in MATLAB) | `isyntax.c:3039` | C-side detail |
| Empty block | `serializedlength = 0`, so only the 32-bit length (4 bytes) | treated as empty when ≤ 8 bytes (`isyntax.c:2577`) | **Real files use ~8 bytes** **[FIC]**, so the production encoder differs from this reference for empty blocks. Confirmed by Philips' `decompress_block.m`, which says its reference decoder "cannot interpret" the 8-byte blocks and substitutes zeros (see [Image_Format_Code_Samples_Claude.md](Image_Format_Code_Samples_Claude.md) §4.3) |
| End-of-stream handling | decoder stops when `bitpos >= length` | stops at `serialized_length` | equivalent in practice |

## 5. Wavelet transform

### 5.1 Structure **[M]**

- **Lifting scheme `rbio2.2` ("LeGall 5/3") with the `Int2Int` option** (`liftwave.m`): predict step `d` with taps `[-0.5 -0.5]`, update step `p` with taps `[0.25 0.25]`. The scaling factors `√2` and `√0.5` are used only for the non-integer variant.
- **`rbio4.4`** (a 9/7) is also defined, selected by the container's wavelet field (2). iSyntax uses wavelet 1.
- **Forward 2-D (`wavelet2dlift`):** zero-pad the image on all sides; then per level, per colour, **vertical first** (on the transpose), **then horizontal** on both the low and high outputs. LL, HL, LH, HH are returned, with LL fed to the next level.
- **Inverse 2-D (`wavelet2dilift`):** **horizontal first** (LH+HH, and LL+HL), **then vertical**. This is the order the C code uses **[code]** (`isyntax.c:1843-1870`).
- **Subband naming:** `HL` comes from the horizontal high-pass of the vertical low-pass output.

### 5.2 Integer rounding is JPEG 2000's reversible 5/3 **[derived from M]**

The code adds `floor(y + 0.5)` to each lifting target (forward) and `ceil(y − 0.5)` (inverse, with negated taps). Working through the two steps:

- Predict (forward): `d = x_odd + floor(−(a+b)/2 + 1/2) = x_odd − floor((a+b)/2)`.
- Update (forward): `s = x_even + floor((d_prev + d_next)/4 + 1/2) = x_even + floor((d_prev + d_next + 2)/4)`.

Both are the JPEG 2000 reversible 5/3 formulas. The C inverse uses `even − ((d1+d2+2)>>2)` and `odd + ((s1+s2)>>1)` **[code]** (`isyntax_dwt.c:99-102`). So the C code is the exact inverse of the Philips forward transform, in the integer arithmetic (the 16-bit overflow caveats in Third_Party_Documentation.md §3.2 still apply). A `rounding='fix'` option also exists in the M code but is not the default.

### 5.3 Zero padding and "valid" coefficients (why the geometry looks the way it does)

- The M forward transform **never extends the signal**. It computes a coefficient only where its whole filter support lies inside the (zero-padded) image. These are the "valid" coefficients (`validwaveletcoefficientposition.m`); the rest are discarded. There is no mirroring at image edges.
- So the padding must be large enough that every pixel of the real image is covered by valid coefficients at all levels (`waveletzeropadding.m`), and **this is what the spec's "padding" and the 4-pixel overlap between neighbouring tiles come from**: tiles overlap because each tile's coefficients need real neighbour data, not an extension.
- **Level 1 first valid positions** (0-based, padded coordinates): **L at 2, H at 1** (filter taps: L spans −2…2, H spans −1…1; derived by hand from `waveletfilterindices.m`). So **in padded-image coordinates L sits on even positions and H on odd positions.**
- **Per level `lvl`:** `firstH = 3·2^(lvl−1) − 2`, `firstL = 2^(lvl+1) − 2`. With the C scale `s = lvl − 1`:
  `firstH = (3<<s) − 2` and `firstL = firstH + (1<<s)`.
  These are `get_first_valid_coef_pixel` and `get_first_valid_ll_pixel` in `isyntax.c:3047-3055` **[code]**. I confirmed the match for all levels 1–9 with a script **[derived]**, using the tap sets above.
- **Padding:** `waveletzeropadding` gives `firstInvalidL + max(fIL) + 1`, where `firstInvalidL = 2^L − 2` and the largest L tap at `L` levels is `2(2^L − 1)`, i.e. `3(2^L − 1)` = `(3<<L) − 3`. That agrees with the spec, with `isyntax.c:692`, and with the valid-data envelopes in every sample FIC **[FIC]**. The tap-extent step is my hand derivation; I did not execute `waveletzeropadding.m`.
- **Consequence for the `cas = 1` question** raised in Third_Party_Documentation.md §5.3: the first valid H is at an odd padded coordinate (1) and the first valid L at an even one (2). The C code's stitched tile buffer puts the low band at odd positions of its own local index. This is consistent with the M code if the buffer starts at an odd padded coordinate. I have not traced the buffer origin, so treat that as **[inferred]**. An independent numpy reconstruction of the real file that used exactly this geometry (H first, low band at odd positions, first valid positions `(3<<s)−2` and `+(1<<s)`) produced a clean, seam-free image ([Example_iSyntax_File_Claude.md](Example_iSyntax_File_Claude.md) §9), which supports the parity and the constants.
- **For JPEG 2000 transcoding:** the Philips transform extends by zeros (and keeps valid coefficients), while JPEG 2000 extends by whole-sample symmetric mirroring. The boundary differences the AI notes try to correct are therefore **not** "zero vs mirror at tile edges" but "continuous-image coefficients vs mirrored per-tile coefficients". The first valid coefficients also sit at odd offsets that the notes' formulas do not assume.

## 6. Colour transform

- **`CompressionCodeSamples/rgb2ycocg.m` and `ycocg2rgb.m` use `fix()`** (truncation toward zero) on `Co/2` and `Cg/2`:
  ```
  forward: Co = R-B; tmp = B + fix(Co/2); Cg = G - tmp; Y = tmp + fix(Cg/2)
  inverse: tmp = Y - fix(Cg/2); G = Cg + tmp; B = tmp - fix(Co/2); R = B + Co
  ```
- **`ViewingCodeSamples/rgb2ycocg.m` and `ycocg2rgb.m` use `bitshift(x, -1)`** on `int16`, with the comment "as seen in H.Malvar and G. Sullivan 2003". These are the same formulas as in the viewing PDF (Third_Party_Documentation.md §2.2).
- **The two Philips versions are different files that disagree on negative odd values** if `bitshift` on negative `int16` is an arithmetic shift (floor). I cannot confirm how Octave or MATLAB defines `bitshift` for negative signed integers, so the exact difference is **[unverified]**.
- **How libisyntax lines up:**
  - the scalar C path (`Cg/2` in C truncates toward zero) matches the **compression** samples' `fix()` (`isyntax.c:1664-1667`);
  - the SIMD and NEON paths (`srai`, floor) match the **viewing** samples' `bitshift` (`isyntax.c:1693-1696`, `:1720-1723`);
  - this settles the question left open in Third_Party_Documentation.md §5.1: **both rounding conventions are Philips'**, so neither C path is "wrong" by the vendor, but they are not the same as each other, and pixel results depend on which path runs.
- **The build system decides which one runs:** CMake (`-march=native`) gives the floor path, the Meson default on x86 the truncating scalar path (CLAUDE.md, "Build / run").
- **Measured on the real file ([Example_iSyntax_File_Claude.md](Example_iSyntax_File_Claude.md) §9.1):** the two roundings give different RGB for 64.5% of pixels (R and B by up to 2 levels, G by up to 1), so the choice is not a corner case.
- The compressor only converts when the image has more than one channel and `ycocg = 1` (default). Y has range 0–255, and Co, Cg have range −255…255.
- Wavelet-domain colour conversion is not addressed anywhere in these files. Quantization and the DWT are applied **after** the colour transform, so the quantized coefficients are YCoCg-R, which is why the Y, Co and Cg columns differ.

## 7. The test container (`compressimage.m` / `decompressimage.m`)

**It is not the iSyntax file format.** It is a minimal wrapper to exercise the codec.

- **Defaults:** `bits 16`, `version 1`, `ycocg 1`, `wavelet 1` (5/3), `dwtlevels 3`, `blocksize [128 128]`, `quantV []`, `deadZoneV []`.
- **Header (little-endian):** version (u8), levels (u8), bit depth (u8), width (u16), height (u16), channels (u8, written as 1 byte even though the comment says u16), `ycocg + 16·wavelet` (u8), block size (2 × u16), LL width and height (2 × u16), then HL, LH, HH width and height (6 × u16) per level.
- **Blocks:** each level's HL, LH, HH are padded to multiples of 4 rows and 8 columns and cut into 128×128 blocks. Each block is 3 coefficients deep. LL is cut and coded only at the top level (one coefficient). Each block goes through `hulskencompress`.
- **Block order in the container:** column fastest, then row, colour, level (finest first), and finally the coefficient group (H blocks, then LL). This is **different** from the real file's cluster-based order (Example_FIC_Files_Claude.md §5.2). There is no XML header, seektable, cluster structure, EOT marker or DICOM item framing.
- **What is the same as the real file:** the whole image is zero-padded once for all levels (`wavelet2dlift`), the codec per block is identical, and LL is stored only at the top level.
- **What differs:** real files use the 5-D block addressing and the block/cluster templates from the FIC header, and group blocks into clusters rather than one raster per level.

## 8. Quality notes on the reference code

These do not change the findings, but they matter if you port or run the files.

- **Debug output left in:** `decodeserializedcoefficients.m` prints in almost every step (statements like `bitpos=bitpos+1` have no semicolon, and `fprintf` calls report progress). `hulskendecompress.m` prints too. Expect very noisy, very slow runs.
- **"TRY 2020-07-09" and "2020-08-09" edits** (preallocation) are scattered through the code. In `bitserialize4x2.m`, `serialized` is created as numeric (`zeros('like', bits)`) but assigned with cell syntax (`serialized{...} = ...`), which MATLAB rejects for a numeric array; it may work only in Octave. I did not run it.
- **Trailing-zero truncation** (`bitserialize4x2.m`, `serialized(1:indLastNonZero)`) removes trailing zero bytes. That makes `numel(input)` smaller than `rows·cols/8 × planes`. The decoder recomputes the length from the masks, and `bitdeserialize4x2.m` raises an error if the lengths differ. The offset bit width (`floor(log2(len)) + 5`) is also computed from the truncated length by the encoder but from the full length by the decoder. This is a latent mismatch for rare blocks **[inferred from reading; not run]**.
- **Single-symbol edge cases:** if there are no zero runs and no escapes, `runlengthbitshistogram` is empty and `zerorunbits` would be empty (an AI note raised this; it is plausible from the code, not run).
- **Performance:** the encoder builds bit strings by repeated concatenation, so it is quadratic in the block size.

## 9. What this adds for your four goals

| Goal | What the Philips code gives |
|---|---|
| 1. File format and FIC header | Little. The container here is a toy. The real (v1) file layout is in [Image_Format_Code_Samples_Claude.md](Image_Format_Code_Samples_Claude.md); Philips' code does not cover the v2 layout at all. |
| 2. 5/3 wavelet | The exact lifting steps and integer rounding; forward order (vertical then horizontal); zero padding and valid-coefficient logic, from which the padding and first-valid-pixel constants follow. |
| 3. Hulsken compression | The complete encoder and decoder, bit by bit (section 4). The only gap is the real iSyntax block framing, which comes from the spec and the C code. |
| 4. Quantizer and Deadzone | Complete: exact algorithm, parameter sets, column/row meaning, the match with the FIC values, and the predicted bitplane masks (section 3). |

## 10. Updates made to the other documents, and next steps

- I updated [CLAUDE.md](CLAUDE.md) (key finding and open questions), [Example_FIC_Files_Claude.md](Example_FIC_Files_Claude.md) (section 6 and the open questions) and [Third_Party_Documentation.md](Third_Party_Documentation.md) (sections 3.4, 5.1 and the gaps table) to reflect section 3.4 and section 6 of this file.
- **Best next test:** take one real `.isyntax` file, decode a level-0 chroma codeblock, and confirm mask bits 1–4 are clear (and bit 1 for Y level 0). A small C tool using `isyntax_hulsken_decompress` would do it, in a new file as `CLAUDE.md` requires.
- **Alternative without real data:** port `quantanddeadzone.m` and `bitserialize4x2.m` to Python (a few dozen lines), push random Laplacian coefficients through them, and check the predicted masks. Python 3 is available here; Octave is not.
- **Folder review status:** all three Philips sample folders are now reviewed ([Viewing_Code_Samples_Claude.md](Viewing_Code_Samples_Claude.md) for CLAHE and sharpening). See [Image_Format_Code_Samples_Claude.md](Image_Format_Code_Samples_Claude.md) for the file-format side, including a confirmation from `decompress_block.m` that the reference decoder cannot read the 8-byte empty blocks that real files contain (this note's section 4.6).
