# Philips image-format code samples (Python + MATLAB): review

Source folder: `/mnt/d/Whole_Slide/iSyntax_Reference_Implementation/ImageFormatCodeSamples/` with `Python/` (4 files, 983 lines) and `Matlab/` (8 files, 516 lines). Philips sample code that reads the iSyntax *file* and reassembles an image from its codeblocks. Companion to [Philips_Compression_Code_Samples_Claude.md](Philips_Compression_Code_Samples_Claude.md) (the codec itself), [CLAUDE.md](CLAUDE.md), [Example_FIC_Files_Claude.md](Example_FIC_Files_Claude.md) and [Third_Party_Documentation.md](Third_Party_Documentation.md).

**Evidence labels:** **[PY]** = read in the Python files; **[M]** = read in the MATLAB files; **[code]** = this repo's C (`file:line`); **[FIC]** = sample headers; **[ran]** = I executed it; **[derived]** = worked out by hand; **[inferred]** = deduction, untested.

## 0. Scope and limits

- **The Python metadata code was run** (section 3.1) on copies of two sample headers (`1.fic`, `PAS.fic`) with the `\r\n\x04` delimiter appended, in a scratch directory. Nothing was written next to your files.
- **Nothing else was run.** There is no real `.isyntax` file (so the extract scripts have nothing to read), and no MATLAB/Octave (so none of the `.m` files ran).
- **These samples only understand the older (v1) file layout.** Every sample header you gave me is data model 100.4 (v2), which this code cannot process past the metadata step (section 3.1).

## 1. Inventory and data flow

| File | Lines | Role |
|---|---|---|
| `Python/isyntax_metadata.py` | 251 | Find the XML header; read dimensions, levels, codeblock size, block count |
| `Python/extract_block_header.py` | 473 | Decode the block header table and seektable; compute each codeblock's file offset and size; write a CSV |
| `Python/extract_codeblocks.py` | 178 | Write each codeblock to a text file; write `properties.csv` |
| `Python/extract_macro_label_image.py` | 81 | Decode the base64 JPEGs of the label and macro images |
| `Matlab/reconstruct_image.m` | 85 | Driver: reassemble the full-resolution image and write a TIFF |
| `Matlab/generate_struct.m` | 157 | Enumerate every possible codeblock file name and coordinate |
| `Matlab/decompress_block.m` | 53 | Call `hulskendecompress` on one codeblock file |
| `Matlab/stitch_LL.m`, `stitch_coef_level.m`, `stitch_color.m`, `stitch_row.m` | 44, 45, 40, 45 | Assemble blocks into LL and HL/LH/HH planes |
| `Matlab/idwt_at_level.m` | 47 | One inverse-DWT level with the padding the geometry needs |

```
iSyntax file --Python--> <name>_block_header.csv, <name>/codeblock_<x>_<y>_<c>_<s>_<k>.ssv, properties.csv
                                        |
codeblock .ssv files + properties.csv --MATLAB--> bigRgbImg.tiff (LZW)
   generate_struct -> decompress_block (hulskendecompress, from CompressionCodeSamples)
   -> stitch_* -> idwt_at_level per level, coarsest first -> crop -> ycocg2rgb
```

The MATLAB side calls `hulskendecompress`, `wavelet2dilift`, `liftwave` and `ycocg2rgb` from the sibling `CompressionCodeSamples` folder (reviewed in Philips_Compression_Code_Samples_Claude.md).

## 2. What the pair demonstrates

Together they are the spec's "reading WSI images" recipe (Third_Party_Documentation.md §2.1) in executable form: find the header, get the block table, extract the codeblocks, Hulsken-decode, stitch spatially, inverse-DWT level by level, convert YCoCg to RGB. They are **demonstration code**: one slide at a time, everything in memory, and text files as the exchange format.

## 3. Python

### 3.1 `isyntax_metadata.py`

- **Finding the header [PY]:** it seeks to every byte offset in the file in turn and reads 3 bytes until it sees `0D 0A 04`. That is correct but needs one seek and read per byte, which is very slow for a multi-gigabyte file (the C code scans a 1 MB buffer for `\x04`, `isyntax.c:3173`). If a payload contained `\r\n\x04` earlier it would be misread, as in the C code.
- **Dimensions [PY]:** each `UFS_IMAGE_DIMENSION_RANGE` is `start step end`, with `size = (end − start + step) / step`. Dimensions are paired with ranges **by index** (x, y, component, scale, waveletcoef).
- **Image size without padding [PY]:** `level_padding = (3 << levels) − 3`, `width = padded_x − 2·level_padding`. This matches `isyntax.c:692-694` **[code]**.
- **Codeblock size [PY]:** from the **first** block header template's first two ranges, `(end − start + step)/step`. For the sample headers that gives 128×128 from `0 2 254`.
- **Running it on the v2 sample headers [ran]:**

  | | `1.fic` | `PAS.fic` |
  |---|---|---|
  | levels | 6 | 8 |
  | image size (valid) | 25222 × 21126 | 101894 × 77318 |
  | codeblock | 128 × 128 | 128 × 128 |
  | blocks | 33,867 | 498,357 |
  | `block_header_table()` | **fails** | **fails** |

  The first four rows agree exactly with Example_FIC_Files_Claude.md §1, which independently confirms those derived sizes and the padding rule against Philips' own code. The last row fails with `AttributeError: 'NoneType' object has no attribute 'text'`, because the header has no `UFS_IMAGE_BLOCK_HEADER_TABLE` (0x2014). The code sets its "found" flag before reading `.text`, so the intended `ValueError` is never raised. Every later step (`extract_block_header`, `extract_codeblocks`, and therefore the MATLAB stage) is blocked for these files.

### 3.2 `extract_block_header.py`: the v1 block table and seektable

**Block header table (attribute 0x2014, base64) [PY]**, all little-endian:

```
u32  total size
repeat per block:
  FFFE E000  u32 item size        // 72 if offset/size are present, 40 if a seektable holds them
  301D 200E  u32 20               // block coordinate: u32 x, y, colour, scale, coefficient
  [ 301D 2010 u32 8  u64 data offset ]   // only when item size = 72
  [ 301D 2011 u32 8  u64 data size   ]   // only when item size = 72
  301D 2012  u32 4  u32 template id
```

`coefficient` is 0 for LL and otherwise 1. The 72-byte form is the "full" header and the 40-byte form the "partial" one, which is what `isyntax_full_block_header_t` / `isyntax_partial_block_header_t` implement **[code]** (`isyntax.h:125-155`).

**Seektable [PY]:** it follows the EOT. It starts with tag `301D 2015` and a u32 total size, then **40-byte entries** (an 8-byte item header, a 16-byte offset element and a 16-byte size element), so `entries = size / 40`. The C struct `isyntax_seektable_codeblock_header_t` is the same layout **[code]**.

**Joining the two (`fill_block_header`, `calculate_block_id`) [PY]**: each block header is mapped to a seektable index. In words:

1. `base_w = ceil(W / (cb_w << L)) << (L−1)` and `base_h` likewise (the tile grid at scale 0), using the padded `W`, `H` and the 128-wide codeblock.
2. Per colour, the seektable order is: all H tiles of scale 0 (raster), then scale 1, … scale L−1, then the LL tiles of the top scale. Tiles at scale `s` form a `(base_w >> s) × (base_h >> s)` grid.
3. `index = Σ_{s' < scale} (base_w>>s')(base_h>>s')` (+ the whole top-scale tile count if the block is LL) `+ tx + ty·(base_w >> scale) + colour · tiles_per_colour`, with `tx = (x − orig_x − offset) / (256 << scale)` and `offset = (3 << scale) − 2` for H or that plus `(1 << scale)` for LL.

This is the same computation as `isyntax.c:3300-3330`, which says it was "adapted from extract_block_header.py" **[code]**. The only difference is cosmetic (the C code sums tile counts up to `scale+1` for LL where Python adds the top-scale count separately). I checked the two by reading them; they are equivalent.

- **Zero-offset entries mean unscanned tissue.** Blocks outside the scanned regions have offset 0 and are filtered out (`extract_codeblocks.py`). That matches the spec ("not written").
- **Coordinates are shifted by the image origin** (`adjust_coordinates_with_orgin`) after the lookup, "to reduce the complication in matlab code".

### 3.3 `extract_codeblocks.py`

- Writes each selected codeblock to `<name>/codeblock_<x>_<y>_<colour>_<scale>_<coef>.ssv`, **one decimal byte value per line of text**. A 40 KB codeblock becomes roughly 120 KB of text. It is a teaching format, not a practical one.
- `x` is the **horizontal** coordinate in these file names. (The MATLAB side calls its row index `x`; see 4.2.)
- `-p` writes `properties.csv`: `num_level, wsi_image_width, wsi_image_height, codeblock_width, codeblock_height, requested_level`. The width and height are the **unpadded** size.
- A level argument of `-1` extracts everything.

### 3.4 `extract_macro_label_image.py`

Selects every `DPScannedImage` whose `PIM_DP_IMAGE_TYPE` is not `WSI` and base64-decodes `PIM_DP_IMAGE_DATA` into `<name>_<TYPE>.jpg`. This is what the stripped sample FICs no longer contain (their data fields are empty).

### 3.5 Limitations that matter to you

- **v1 layout only.** There is no code for cluster header templates (0x2016), the cluster header table (0x201f), valid-data envelopes, or the `UFS_IMAGE_DIMENSIONS_IN_CLUSTER` ordering. **Philips' reference code does not describe the v2 layout at all.** For v2 the only implementation available to us is this repo's C code (Example_FIC_Files_Claude.md §5 and §7).
- It reads **only** the first block-header template for the codeblock size, and assumes one template.
- It assumes a seektable follows the header when the 40-byte form is used.
- Error handling is minimal, and `__scanned_image` re-parses the entire XML on every call (slow for the multi-megabyte v2 headers).

## 4. MATLAB

### 4.1 `reconstruct_image.m`

1. Reads `properties.csv`, and calls `generate_struct` for the list of all possible codeblocks.
2. Stitches and decodes the **LL** blocks (coefficient 0) with `stitch_LL`.
3. For levels from the **coarsest to the finest**, stitches the H blocks of that level and runs one inverse-DWT step (`idwt_at_level`) on `LL` plus them, replacing `LL` with the larger result.
4. Crops the accumulated padding off the top-left and the image size off the bottom-right.
5. Converts YCoCg to RGB (`ycocg2rgb`, the `fix`-based compression version) and writes `bigRgbImg.tiff` with `bigimage` and LZW compression.

A comment in the file says the 2020-07-07 run "generates a one-level tiff with weird color (red background)". It is not resolved in the sample.

### 4.2 `generate_struct.m`: block coordinates, enumerated rather than read

It does not read the block table. It **generates** every codeblock file name from the geometry, then relies on `decompress_block` to turn missing or empty files into zeros. This confirms the geometry independently of the header:

- **H blocks at scale `s` start at coordinate `c_s`** with `c_0 = 1` and `c_{s+1} = 2·c_s + 2`, giving **1, 4, 10, 22, 46, 94, 190, …**, which is `(3 << s) − 2`.
- **LL blocks (top level only)** start at `((c_L − 2)/2) + 2^(L−1)`, which equals `first_valid_coef + (1 << s)` at `s = L − 1`.
- **The step between blocks** is `block_size · 2^(s+1)` (256·2^s pixels).

Both match `get_first_valid_coef_pixel` / `get_first_valid_ll_pixel` **[code]** (`isyntax.c:3047-3055`) and the cluster-range starts in the sample FICs (1, 22, 190 for scales 0, 3, 6) **[FIC]**. I verified the closed forms by hand **[derived]**; the MATLAB itself was not run.

Other points:

- **Naming:** file names are built as `codeblock_<coef_y>_<coef_x>_…`, where MATLAB's `coef_x` is the **row** index and `coef_y` the **column** index. So the name's first number is the horizontal coordinate, which agrees with the Python file naming. The swap of "x" and "y" inside the MATLAB code is internal.
- **Grid size** is computed from `wsi_image_width/height`, which `properties.csv` holds **without padding**, while the Python block-ID code uses the **padded** size. The two give the same tile count unless the unpadded and padded widths fall on opposite sides of a multiple of `128 << L`. For `1.fic` and `PAS.fic` they do not **[derived]**, but it is a latent mismatch **[inferred]**.
- The blocks are sorted by y, then x, then colour.

### 4.3 `decompress_block.m`

- `hulskendecompress(double(load(file)), 16, [w h(, 3)], 1)`: **16 coefficient bits and compressor version 1 are hard-coded.** LL blocks are decoded as 2-D, H blocks as 3 planes.
- **Any failure becomes a zero block.** The comment says: "The 8 bytes codeblocks from an iSyntax file refer to the non-tissue area … Reference code has a limitation that it cannot interpret these codeblocks." Two consequences:
  - It **confirms** that the production encoder writes an 8-byte block for an all-zero block, which this reference decoder cannot read (my note in Philips_Compression_Code_Samples_Claude.md §4.6 about the empty-block size). The C code treats ≤ 8 bytes as empty (`isyntax.c:2577`).
  - The comment says these are "background". In real files they are **also** the tissue blocks whose coefficients were all zeroed by the deadzone (Example_FIC_Files_Claude.md §6.4: up to 96% of chroma blocks at scale 0). Decoding them as zeros is still correct, but "non-tissue" understates what they are.

### 4.4 Stitching

- `stitch_row` lays blocks side by side (a strip), `stitch_color` stacks strips in order of the row index, `stitch_coef_level` splits the three subband planes (HL, LH, HH) by colour, and `stitch_LL` does the same for LL. Plane order in an H block is HL, LH, HH, as in the C code.
- **Known unresolved bug [M]:** a `TODO: 2020-07-09` in `stitch_row.m` records the error "Unable to perform assignment because the size of the left side is 128-by-128-by-3 and the size of the right side is 128-by-128." `row_x` always has three planes, but an LL block is 2-D, so stitching LL as written fails. `stitch_LL` runs first in `reconstruct_image`, so the shipped driver would fail at that step **[inferred from reading; not run]**. (An earlier AI note reached the same conclusion and suggested a fix; I did not test it.)
- `row_x` is allocated `(codeblock_width, codeblock_height·n, 3)` but filled with a column stride of `codeblock_width`. It works only because the blocks are square.

### 4.5 `idwt_at_level.m` and the padding bookkeeping

- Each level builds `HL{1}, LH{1}, HH{1}` for one level and calls `wavelet2dilift` (integer 5/3, `Int2Int`).
- It pads LL with 3 zero rows/columns at the top-left and 1 at the bottom-right, and pads HL (right), LH (bottom) and HH (both) with one zero row/column. This is how the stored valid-coefficient grids (L from coordinate 2, H from coordinate 1 at level 1) are made the same size for the inverse transform.
- The LL padding is applied only when **both** width and height mismatch (`&&`), so a non-square case where only one dimension differs would skip it **[inferred]**.
- **Final crop:** `top_left_pad = Σ_{levels ≥ 1} 3·2^level = 3(2^L − 2)`, and the last (level 0) inverse step crops 3 itself inside `wavelet2dilift`. In total `3(2^L − 2) + 3 = 3(2^L − 1) = (3 << L) − 3`, which equals the padding in the spec, in `isyntax.c:692` and in every sample FIC **[derived]**.

## 5. What this adds, and cross-checks against libisyntax

| Item | Philips sample | libisyntax C | Result |
|---|---|---|---|
| Header delimiter `0D 0A 04` | Python, byte scan | `isyntax.c:3173` (`\x04`) | agree |
| Valid image size = padded − 2·`(3<<L) − 3` | `wsi_image_dimension` | `isyntax.c:692-694` | agree [ran for the numbers] |
| Block header item sizes 72 / 40, tags `2010/2011/2012/200E` | `extract_block_header.py` | `isyntax.h:125-155` | agree |
| Seektable 40-byte entries, tag `0x2015` | Python | `isyntax.c:3349` | agree |
| Block ID / seektable index | `calculate_block_id` | `isyntax.c:3300-3330` | agree (C is a port) |
| First valid coefficient `(3<<s)−2`, LL `+ (1<<s)` | Python and `generate_struct.m` | `isyntax.c:3047-3055` | agree |
| Colour transform | `ycocg2rgb` (`fix`) | scalar path `/2`, SIMD `>>1` | scalar matches, SIMD differs (see Philips_Compression_Code_Samples_Claude.md §6) |
| Empty blocks | 8 bytes, reference decoder fails, treated as zeros | ≤ 8 bytes, zeros | agree in effect |
| v2 clusters, envelopes, dims-in-cluster | **not covered** | implemented | C only |
| Compressor version | hard-coded 1 | read from header (16→v1, 19→v2) | C is more general |

**Seektable order versus the spec.** The Python index order is per colour: scale 0 first, scale L−1 last, LL after. The spec's Figure 9 says "per level: NumLevels to 0 (increasing resolution)" for the *packaging* order and the text says IDs count "from the base level to the top level". The Python shows the **ID order** is base → top. The **file order** is a separate thing, determined by the offsets in the seektable (or, in v2, by the clusters). So the two statements in the spec are about different orders, as I suggested in Third_Party_Documentation.md §2.1; this code now supports that reading.

## 6. Open questions

1. ~~Does a v1 sample file exist that these scripts can read?~~ **Resolved:** `example_isyntax_files/1.isyntax` is a real v1 file (5.0), and the Philips Python reads it correctly (see [Example_iSyntax_File_Claude.md](Example_iSyntax_File_Claude.md) §2). `1.fic` is its v2 (100.4) counterpart. Original text: The earlier chat described `1.txt` as "the smallest, version-1 reference file", but `1.fic` here is interface version 100.4 and has no block header table, so the Philips Python cannot read it past the metadata step. Either that SDK file is a different file, or "version-1" meant something else (compressor v1 is true of all six samples). Worth checking which header a real v1-data-model file carries.
2. **How does Philips' own reference code read the v2 layout** (cluster table, envelopes)? Nothing in either sample folder does. The C code is the only description.
3. **Does the shipped MATLAB driver run at all?** (The geometry it implements is confirmed: an independent numpy reconstruction of the real v1 file using the same coordinates gives a clean image, see [Example_iSyntax_File_Claude.md](Example_iSyntax_File_Claude.md) §9. The LL stitching bug is separate.) Original question: The LL stitching TODO suggests not. Running it needs Octave or MATLAB and a v1 file.
4. **Grid-size basis** (padded vs unpadded width) in `generate_struct.m` vs the Python: harmless on the samples, latent in general.
5. **The "red background" note** in `reconstruct_image.m` is unexplained. A candidate cause is the zero fill for missing LL (Y = 0 and chroma 0, versus libisyntax's white Y = 255 dummy for missing top-scale LL), but I did not trace it.

## 7. Next steps

- A real **v1-layout** `.isyntax` file would let the Python scripts and the MATLAB driver run end to end, and would give codeblocks for the quantizer/deadzone mask check (Philips_Compression_Code_Samples_Claude.md §3.3).
- A small Python port of the block-ID computation could be checked against the C code on any file that has a seektable. I can write it as a new tool file if you want.
- `ViewingCodeSamples/` (CLAHE, sharpening) is now reviewed in [Viewing_Code_Samples_Claude.md](Viewing_Code_Samples_Claude.md).
