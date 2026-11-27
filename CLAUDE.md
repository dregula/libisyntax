# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

Fork of `libisyntax` (https://github.com/dregula/libisyntax, upstream by Pieter Valkema), the core decoder behind the SlideScape whole-slide-image viewer. It reads Philips iSyntax files. In this fork the library is used as a **reference implementation for studying the iSyntax format**, not only as a decoder.

## Goals of this fork (read this first)

The owner wants to understand, at the bit level:

1. The iSyntax file format and its XML header (called the **FIC header** in this project).
2. The reversible 5/3 (LeGall) wavelet and the IDWT convolution.
3. Hulsken compression of wavelet coefficients (Huffman + zero-run RLE over bitplanes).
4. How the FIC header's **Quantizer** and **Deadzone** parameters show up in the coefficient bitplanes, before Hulsken compression.

Prefer explaining and annotating over refactoring. Do not change decoder behavior unless asked: the regression tests compare decoded pixels against known output.

### Path convention

The parent folder of this repository (the one that holds `libisyntax/`, `example_isyntax_files/`, `c_decoder_output/` and `third_party_documentation/`) is written as `$SOURCEROOT` in this repository. Each developer sets that environment variable outside the repository (it is deliberately not committed). Relative paths such as `../example_isyntax_files/` also work from the repo root.

### Companion material (outside this repo)

| Path | Contents |
|---|---|
| `../third_party_documentation/` | Specs and papers: Philips iSyntax format PDF (+ markdown conversion), Hulsken 2016 paper, Hulsken patent US20180255321A1, webinar slides, and earlier `*_Claude.md` notes (YCoCg-R→RCT, codeblock ordering, JPEG2000/DICOM/OME transcoding comparisons). Check these before re-deriving anything. |
| `../example_FIC_files/` | Extracted XML headers: `1.fic 2.fic 4.fic 9.fic 10.fic PAS.fic`. |
| `../example_isyntax_files/` | `1.isyntax`: one complete real file (88.9 MB, data model v1, same scan as `1.fic`). `3` to `13.isyntax` (no 2; 0.5 to 1.7 GB, 7 to 13 Gpx, 8 or 9 levels): more v1 files (interface 5.0, Quality=2). Block tables of all are parsed; six (1, 3, 4, 5, 11, 12) are fully analysed. These larger files need the strip-wise script (`reconstruct_strips.py`); `reconstruct_np.py` would need 200 GB or more. |

More specs and FIC samples are expected later. The `.fic` files are plain XML text, so use `grep -a`.

## Key finding: the library never uses Quantizer/Deadzone

- `isyntax.c:629-632` handle `DP_WAVELET_QUANTIZER_SETTINGS_PER_COLOR` (0x1019), `..._PER_LEVEL` (0x101A), `DP_WAVELET_QUANTIZER` (0x101B) and `DP_WAVELET_DEADZONE` (0x101C) as **empty `{} break;` stubs**. The values are not stored anywhere.
- The decode path has **no dequantization step**. Coefficients go: bitplanes → `i16` → signed-magnitude to two's-complement → IDWT, with no scale factor or deadzone reconstruction offset.
- **Resolved by Philips' own encoder** (`CompressionCodeSamples/quantanddeadzone.m`, see [Philips_Compression_Code_Samples_Claude.md](Philips_Compression_Code_Samples_Claude.md) section 3): the encoder zeroes coefficients with `abs(c) <= deadzone`, then truncates the rest toward zero to a multiple of `2^Quantizer`. Nothing is divided, so the decoder needs no dequantization. Effect on the bitplanes: the lowest `Quantizer` magnitude planes are always empty and never stored, and small values never occur. LL is never quantized, and neither are levels beyond the table.
- The spec's "greater or equal" deadzone wording is **wrong** (the code uses `<=`).
- `Quality=2` in the derivation string is Philips' "Q2" parameter set, and the sample FIC values are exactly that set. Column 1 of the Philips matrices is Y, column 2 is Co and Cg. Header level 0 is the first (finest) DWT level, i.e. scale 0 in this code.
- **Verified on real data** ([Example_iSyntax_File_Claude.md](Example_iSyntax_File_Claude.md)): in `example_isyntax_files/1.isyntax` (88.9 MB, interface 5.0 = the v1 export of the same scan as `1.fic`), all 33,867 codeblocks decode, and **none** has a forced-clear mask bit set (v1 numbering: Y level 0 bit 1 clear; chroma level 0 bits 1-4; chroma level 1 bits 1-2). The smallest nonzero coefficient equals the deadzone (2, 16, 8). **Repeated on files 3, 4, 5, 11 and 12** (7 to 13 Gpx, all v1, Quality=2): every codeblock decodes and no forced-clear bit is set ([Multi_File_Strip_Analysis_Claude.md](Multi_File_Strip_Analysis_Claude.md)).

### What the sample headers contain (all 6 files are identical here)

`DERIVATION_DESCRIPTION = "PHILIPS UFS V1.8.6614 | Quality=2 | DWT=1 | Compressor=16"`, `LOSSY_IMAGE_COMPRESSION_RATIO=15`, method `PHILIPS_DP_1_0`. Block compression method is 16, which means Hulsken **compressor v1**.

Nesting (inside the WSI `DPScannedImage`, after `DP_COLOR_MANAGEMENT`; there is no `DPImagePostProcessing` wrapper in these samples): `DP_WAVELET_QUANTIZER_SETTINGS_PER_COLOR[]` → `DP_WAVELET_QUANTIZER_SETTINGS_PER_LEVEL[]` → `{QUANTIZER, DEADZONE}`. Observed (Quantizer, Deadzone) pairs:

```
group A: (1,2) (0,0)
group B: (4,16) (2,8) (0,0)
```

Resolved by the spec (Tables 11-12): group A is **Y**, group B is **Co and Cg**; entries are per level, level 0 first; quantizer means "rounded to 2^Q" and deadzone means small coefficients are zeroed. Details, and block-size evidence that supports this, are in [Example_FIC_Files_Claude.md](Example_FIC_Files_Claude.md) section 6.

Still open:
- Group B (chroma) has 3 header entries but Philips' Q2 matrices have only 2 rows. The third `(0,0)` entry changes nothing in the data; why it is listed is unknown.
- `grep -a -A1 'Element="0x101[bc]"' *.fic` extracts the values.

## Build / run

Two build systems (keep both working): CMake and Meson. Output binaries land in the repo root (CMake) or `builddir/` (Meson).

**Baseline build (decided 2026-11-26): Linux + Meson with `-mssse3`, directory `build-ssse3/`.** Use it for running and testing, and say so when reporting output; it selects the SIMD floor (`>>1`) colour path and matches the cmake-native output byte for byte (see the notes below). Create and build it with `meson setup build-ssse3 -Dc_args=-mssse3 && ninja -C build-ssse3`. `builddir/` (plain `-O3`, scalar `/2` colour path) is only for comparing the two rounding paths.

```bash
# Meson (a configured ./builddir exists)
meson setup builddir         # first time; add -Dtests=true for tests
ninja -C builddir
meson test -C builddir          # needs a test slide, see test/fetch.py
meson test -C builddir smoke_example_runs_no_args   # single test by name (names are in meson.build: `smoke_example_runs_with_test_slide_showing_levels` and `regression_example_tile_3_5_10_pixel_check` need the test slide; `smoke_thread_test` is skipped on macOS)

# CMake
cmake -B cmake-build -DCMAKE_BUILD_TYPE=Release && cmake --build cmake-build
# Debug + sanitizers: -DCMAKE_BUILD_TYPE=Debug -DENABLE_MEMORY_CHECK=ON
```

Targets: `libisyntax` (shared lib), `isyntax_example` (tile/region dump, PNG), `isyntax-to-tiff` (needs libtiff), `isyntax-dirwalk`, `thread_test`.

Notes:
- **The two build systems select different code paths (measured, see Example_iSyntax_File_Claude.md §9.1 and `builddir/compile_commands.json`).** CMake adds `-march=native`, so on this machine (ssse3 + avx2) `__SSSE3__` is defined and the YCoCg->RGB conversion runs the **SIMD floor (`>>1`)** path. The Meson build adds **no `-march` flag on x86** (only on Apple Silicon): the existing `builddir/` was compiled with plain `-O3`, only `__SSE2__` is defined, and the conversion runs the **scalar truncating (`/2`)** path. The wavelet lifting is integer and bit-exact in every variant, but the two colour paths give different RGB for about 64% of the pixels of the real file (R and B by up to 2 levels, G by up to 1). **Always record which build produced a decoder output.** To get the SIMD path from Meson use `meson setup build-ssse3 -Dc_args=-mssse3`.
- CMake also changes `VREG_INT_COUNT` (AVX2/SSE lifting paths) and uses NEON on Apple Silicon.
- `.gitignore` excludes `*.json *.xml *.isyntax build/ builddir/ .claude/ reconfigure/`. A `.gitignore` edit is currently uncommitted.
- There is no real `.isyntax` test data in the repo. `test/fetch.py` downloads it and `test/expected_output/` holds the references. One real file is available outside the repo: `$SOURCEROOT/example_isyntax_files/1.isyntax` (v1 data model, same scan as `1.fic`). Read-only Python analysis tools for such files are in `tools/v1_analysis/` (use `reconstruct_strips.py`, bounded memory about 2 to 3.5 GB, for the large files) (plus `tools/run_c_decoder.sh`, see "Running the C decoder" below; the transcoding experiments are listed in that folder's README.txt); they include an independent numpy reconstruction of the full image (needs the repo's `.venv`, Python 3.14 + numpy).

## Running the C decoder and storing its output for Claude

**Who runs what.** Earlier Claude sessions were blocked from executing compiled binaries (the permission system blocked it twice); a session on 2026-11-26 did run `isyntax_coefficients` from `build-ssse3/` without trouble. If a run is blocked, do not try to work around it: **you run the C decoder in a normal terminal, store the output in the folder below, and tell Claude the run folder**. Claude then reads the files (PNG, JPEG, logs, checksums) and compares them with the independent numpy reconstruction. Claude *can* run Python (`.venv/bin/python3 -I`, numpy 2.5) and the read-only tools in `tools/v1_analysis/`.

**Where the output goes** (outside the repository, next to the other example folders; override with `C_DECODER_OUTPUT`):

```
$SOURCEROOT/c_decoder_output/
  reference_numpy/                     Y.npy Co.npy Cg.npy  (int16 planes of the numpy reconstruction, ~3.2 GB; see step 1)
  <YYYYMMDD_HHMMSS>_<variant>_<file>/  one folder per run, never overwritten
      RUN.txt         date, variant, build flags, compiler, git commit, CPU flags, input name/size/sha256, every command, file list
      SHA256SUMS.txt  checksums of every file in the run folder
      levels.log      isyntax_example <file>  (metadata, level sizes, tile grid), with timing
      label.jpg macro.jpg  and their logs
      tile_L<level>_<x>_<y>.log  and  tiles/L<level>_x<x>_y<y>.png   (256x256 RGBA tiles)
      full_lzw.tif    only if FULL_TIFF=1 (isyntax-to-tiff, lossless, large)
```

**Step 1: one-time reference for comparison (about 3 minutes, about 15 GB RAM).** From the repository root:

```bash
mkdir -p ../c_decoder_output/reference_numpy
.venv/bin/python3 -I tools/v1_analysis/parse_v1.py ../example_isyntax_files/1.isyntax ../c_decoder_output/reference_numpy/blocks.csv
SAVE_PLANES=1 .venv/bin/python3 -I tools/v1_analysis/reconstruct_np.py ../example_isyntax_files/1.isyntax \
    ../c_decoder_output/reference_numpy/blocks.csv ../c_decoder_output/reference_numpy 32
```

(Claude can do this step itself.) The planes cover the valid region (25222 x 21126), i.e. padded coordinates minus the 189-pixel padding.

**Step 2: run the C decoder.** `tools/run_c_decoder.sh <variant> [isyntax_file]` builds the chosen variant, runs the example programs, and writes the run folder:

| variant | build | YCoCg->RGB path expected |
|---|---|---|
| `meson-default` | existing `./builddir` (`ninja -C builddir`) | scalar, truncating `/2` |
| `meson-ssse3` | `./build-ssse3` (`meson setup build-ssse3 -Dc_args=-mssse3`) | SIMD, floor `>>1` |
| `cmake-native` | `./cmake-build` (CMake, `-march=native`) | SIMD, floor `>>1` |

Recommended: run `meson-default` first (it is the build that already exists) and then `cmake-native`; the two should differ in about 64% of pixels, which tests the rounding finding.

```bash
cd $SOURCEROOT/libisyntax
tools/run_c_decoder.sh meson-default                 # default input: ../example_isyntax_files/1.isyntax
tools/run_c_decoder.sh cmake-native
# optional: more tiles, or the full lossless TIFF (large and slow)
EXTRA_TILES="0 40 40;1 20 20" tools/run_c_decoder.sh meson-default
FULL_TIFF=1 tools/run_c_decoder.sh cmake-native
```

The script writes level 0 tiles (35,35) (36,35) (35,36) (36,36) at a tissue edge, (60,40) tissue interior, (10,10) blank glass, (97,50) the edge of the scanned area, level 2 tile (15,9), level 3 tile (7,4), and **all 12 tiles of level 5** (the whole slide at 1/32 resolution, 4 x 3 tiles). Tile coordinates are `<level> <tile_x> <tile_y>` in that level's tile grid; level 0 is full resolution.

**Manual equivalent** (one tile; the example program is `builddir/isyntax_example`, or `./isyntax_example` for CMake):

```bash
builddir/isyntax_example FILE                              # metadata + levels
builddir/isyntax_example FILE <level> <tile_x> <tile_y> out.png
builddir/isyntax_example FILE label label.jpg              # and: macro macro.jpg
builddir/isyntax-to-tiff FILE out.tif --compression LZW    # lossless; JPEG is lossy, do not use it for comparisons
```

**Step 3: tell Claude the run folder**, e.g. "C decoder output is in $SOURCEROOT/c_decoder_output/20261125_101500_meson-default_1". Claude will run

```bash
.venv/bin/python3 -I tools/v1_analysis/compare_c_tiles.py <run_folder> ../c_decoder_output/reference_numpy
```

which, for every level 0 tile, finds the pixel offset between the C tile and the numpy planes, and reports the percentage of identical RGB pixels under the floor and truncating roundings. Expected: about 100% for the variant that matches the build's colour path (and about 35% for the other), with a best offset of (dx=-1, dy=-1) relative to the script's prediction. **Measured (run 20261125_141312_cmake-native, 7 level 0 tiles, all 100.00% identical under the floor path):** the real mapping is `valid column = 256*tx + j - 189` (padded column `256*tx + j`, no `+1`); the script's prediction `+ 1` was off by one and `compare_c_tiles.py` searches around it. The level 2, 3 and 5 tiles and the label/macro images are kept for later checks (level 5 shows the slide at 1/32 resolution; the macro image is the scanner thumbnail).

**Measured results (2026-11-25, both runs in `../c_decoder_output/`):** `20261125_141312_cmake-native_1` matches the numpy reconstruction 100.00% under the floor (`>>1`) path and `20261125_142012_meson-default_1` matches 100.00% under the truncating (`/2`) path, on all 7 level 0 tiles (best offset (-1,-1) in both). Comparing the two C runs directly, 62-68% of pixels differ in busy level 0 tiles (RGB difference 0 to 2), 29% in the tile at the edge of the scanned area, 2% in blank-glass level 5 tiles; this reproduces the 64.5% predicted from the numpy planes. So the two build systems really do give different pixels, and the finding in the Build notes is confirmed on C output. `20261125_142126_meson-ssse3_1` (Meson with `-mssse3`) also matches 100.00% under the floor path on all 7 level 0 tiles, and all its 21 tile PNGs plus label.jpg and macro.jpg are byte-identical to the cmake-native run: so the colour path is decided by `__SSSE3__` alone, and the other CMake/Meson differences (`-march=native`, AVX2 lifting) change no output byte. Only level 0 tiles were checked against numpy; level 2, 3, 5 and label/macro are kept but unchecked.

**Rules.** Keep every run folder (they are small without the TIFF). Never edit the input file or the decoder to produce output; if you change the code, say so in the message, because `RUN.txt` records the git commit and the number of dirty `src/` files. `RUN.txt` must stay with the images: the build variant is the key to interpreting any pixel difference. If a build step fails (for example `meson` or `libtiff` missing; the TIFF tool is only built when libtiff is found), the log is in the terminal; `isyntax_example` alone is enough for everything above.

## Source map (`src/`)

**Line numbers.** `src/isyntax/*.c` and `isyntax.h` now carry detailed explanatory comments (381 lines added, comments only; the code is token-for-token identical to commit `92471db` once comments are stripped, checked with `gcc -fpreprocessed -dD -E -P`). All `isyntax.c:NNNN` / `isyntax_dwt.c:NNN` references in this file and in the companion documents were written **before** that and refer to commit `92471db`. To see the code at those line numbers use `git show 92471db:src/isyntax/isyntax.c | sed -n 2551,2560p`; in the working tree search for the function name instead (the long description of the Hulsken decoder sits directly above `isyntax_hulsken_decompress`, the format and geometry notes at the top of `isyntax.h`).

| File | Role |
|---|---|
| `libisyntax.h/.c` | Public C API: open, level/tile queries, `libisyntax_tile_read`, `read_region`, label/macro images, metadata getters, cache. **Added in this fork:** `libisyntax_tile_read_coefficients` (stored H coefficients at any level, LL at the coarsest level; no IDWT, no colour conversion, no cache, no lock; thread-safe, each caller thread needs `libisyntax_thread_init`), implemented by `isyntax_read_tile_coefficients` in `isyntax_reader.c`; example `src/examples/isyntax_coefficients.c` (**run 2026-11-26 on `build-ssse3`**: output identical to the independent Python decoder on 3 tiles of `1.isyntax` at levels 0, 2 and 4; serial and 8 to 16 thread reads identical on 12,288 tiles of `1.isyntax` and `5.isyntax`). Purely additive: existing code paths are unchanged, so the comments-only identity check no longer holds for `isyntax_reader.c`. |
| `isyntax/isyntax.h` | All core structs and the header enums. **Read first.** |
| `isyntax/isyntax.c` (3.7k lines) | XML header parsing, Hulsken decoder, IDWT tile driver, YCoCg→RGB, `isyntax_open`. |
| `isyntax/isyntax_dwt.c` | 5/3 inverse lifting, adapted from OpenJPEG (BSD-2). **`#include`d into `isyntax.c`**, not compiled separately. |
| `isyntax/isyntax_reader.c` | Cache-backed tile reading, used by `libisyntax.c`. |
| `isyntax/isyntax_streamer.c` | Legacy SlideScape threaded streaming. **Not in either build**. Useful as documentation of chunk and prefetch logic. |
| `platform/`, `utils/`, `third_party/` | Threads/work queue, arenas, SIMD intrinsics, block allocator; `yxml` (XML), `stb_*`, `ltalloc`. |

## File format as implemented

```
XML header (FIC) | EOT "\r\n\x04" | [Seektable] | Codeblocks
```

- `isyntax_open` (`isyntax.c:3110`) reads 1 MB chunks and scans for `\x04` to find the end of the XML. The XML is parsed incrementally with `yxml`.
- Data model version: `<100` is iSyntax v1 (block header table, `UFS_IMAGE_BLOCK_HEADER_TABLE` 0x2014, base64 DICOM-tagged records; this is the only layout Philips' Python/MATLAB samples understand). `>=100` is v2 (cluster header templates, cluster header table 0x201f, valid-data envelopes), which is documented only by this repo's C code and the sample FICs. All six sample FICs are 100.4. The one real file, `example_isyntax_files/1.isyntax`, is v1 (5.0); its seektable and the v2 table of `1.fic` agree on all 33,867 blocks.
- Image types in one file: `MACROIMAGE`, `LABELIMAGE` (base64 JPEG) and `WSI` (wavelet pyramid).
- Dimensions are 5-D: `x, y, color(Y/Co/Cg), scale(0-8), waveletcoeff(LL/LH/HL/HH)`. A block header template gives a codeblock's `block_width/height` (128), `color`, `scale`, and `waveletcoeff` (1 = LL; 3 = HL+LH+HH).
- `PER_LEVEL_PADDING = 3`. Level-0 padding is `(3 << level_count) - 3` per side. Image width is padded width minus 2×padding (`isyntax.c:690`).
- **Codeblock** = one color × scale × tile position. A `coefficient==1` block carries 3 subbands (H); `coefficient==0` is the LL block, which exists only at the top scales. Codeblocks are grouped in **chunks** spanning up to 3 scales (e.g. levels 6,7,8). See the example table at `isyntax.c:2435`. Blocks of size ≤ 8 bytes are empty and decode to zeros.
- Color is **YCoCg** (Y, Co, Cg). Reconstruction is at `isyntax.c:1663`.
- `UFS_IMAGE_BLOCK_COMPRESSION_METHOD` (0x200F): 16 → compressor v1, 19 → v2 (`isyntax.c:727`).

## Hulsken decoder: `isyntax_hulsken_decompress` (`isyntax.c:2551`)

Per-codeblock pipeline, all little-endian and LSB-first:

1. **Header**
   - v1: `u32 serialized_length`, `u8 zerorun_symbol`, `u8 zero_counter_bits`.
   - v2: 1 or 3 `u16` **bitmasks** (one per coefficient), then the same two bytes, then a **bitplane seektable** (`log2(len)+5` bits per pointer). The serialized length is derived as `popcount(masks) × w·h/8`.
2. **Huffman table**: serialized tree, where bit 0 = descend and bit 1 = leaf followed by an 8-bit symbol. Codes are read LSB-first. Decoding uses an 11-bit fast table plus a slow path for longer codes.
3. **Huffman + zero-run RLE** decode to a byte buffer that is the concatenation of **bitplanes** (`w·h/8` bytes each, 1 bit per coefficient).
   - The `zerorun_symbol` is followed by a `zero_counter_bits` counter.
   - Counter 0 means an escaped literal symbol.
   - Consecutive zero-run symbols concatenate their counters into a wider count.
   - v2 stores count−1.
4. **Bitmask**: bit `i` says whether bitplane `i` of that coefficient is present at all. Absent planes are all-zero and are not stored. In v1 the masks are at the end of the decompressed data when planes are missing. In v2 they are in the header.
   - v1 plane order per coefficient: `bit0 = sign`, then LSB→MSB (`shift = i-1`).
   - v2 plane order: sign, then MSB→LSB, interleaved across the coefficients (`shift = 15-i`).
5. **Unpack** bitplanes into `u16` words (SIMD bit-spread at `isyntax.c:2994`), reorder 16 consecutive values into 4×4 areas then raster rows (`isyntax.c:3014`). Each bitplane byte is one 4×2 area (the webinar's "symbol"), so 16 consecutive values (two bytes) form a 4×4 area. The "snake" order is not visible in the code, which copies rows without reversing them. See Third_Party_Documentation.md section 2.4, then convert **sign-magnitude → two's complement** (bit 15 = sign).

Output is `coeff_count × w × h` `i16` values. In `coeff_h`, order is HL, LH, HH (offsets 0, +stride, +2·stride).

**This stage (4, plus the sign-magnitude form) is where quantization and deadzone must be visible**: which bitplanes are present, and how sparse the high planes are.

## Wavelet / IDWT path

- `icoeff_t` is `i16` (`DWT_COEFF_BITS 16`). Lifting is integer reversible 5/3, as in JPEG 2000:
  `s[n] = even − ((d[n-1]+d[n]+2) >> 2)`, then `odd' = d[n] + ((s[n]+s[n+1]) >> 1)`.
  Boundary cases use symmetric extension (see `isyntax_dwt.c` `cas0/cas1`).
- iSyntax always uses `cas = 1`: the low band is at odd positions, which differs from the usual JPEG 2000 tile layout. SIMD processes 8 columns (SSE2) or 16 (AVX2) at once for 16-bit.
- `isyntax_idwt` (`isyntax.c:1821`) takes a `2Q×2Q` buffer laid out as `LL | HL / LH | HH`. It runs rows first, then columns.
- `isyntax_idwt_tile_for_color_channel` (`isyntax.c:1971`) builds that buffer for one tile and one color:
  - Quadrant size is `block + 4 + 4` (`ISYNTAX_IDWT_PAD_L/R`). Margins are copied from the **8 neighboring tiles**, which is why neighbor tiles must be loaded.
  - Missing H coefficients are filled with a zero "black" block. Missing LL is filled with a constant.
  - The valid output starts at pixel `ISYNTAX_IDWT_FIRST_VALID_PIXEL = 7`. The result is a `2·block` tile (256 px for 128 blocks).
  - Each IDWT output is split into 4 child tiles (`isyntax.c:2330`). These become the **LL input of the next finer scale**. Scale 0 is full resolution and the pyramid is built top-down from `max_scale`.
- Quirks to keep in mind:
  - After the IDWT, the **Y channel is converted to an absolute value** (`isyntax.c:2399`). Co and Cg stay signed. Measured on the real file ([Example_iSyntax_File_Claude.md](Example_iSyntax_File_Claude.md) §9): Y never goes negative (range 48 to 257), so this step changes nothing on that slide.
  - The scalar `ycocg_to_rgb` uses `/2` (truncation) and clamps only at the top (`ATMOST(255, ...)`); the SIMD paths use `>>1` (floor) and saturate both ends. **Measured on the real file: 64.5% of pixels differ between the two roundings (R and B by 0-2 levels, G by 0-1)**, so which path runs decides most pixels. No pixel went below 0 on that slide, so the missing lower clamp did not matter there. See [Example_iSyntax_File_Claude.md](Example_iSyntax_File_Claude.md) §9.1. **On the larger slides the share differs (17 to 27 % on glass-heavy slides) and files 12 and 4 have 42 and 19 pixels below 0 before clamping**, where the scalar path (no lower clamp) would wrap ([Multi_File_Strip_Analysis_Claude.md](Multi_File_Strip_Analysis_Claude.md) §2). Y stays non-negative on every slide.

## DICOM WSI transcoding exploration (summary; details in [iSyntax_to_DICOM_WSI_Transcoding_Claude.md](iSyntax_to_DICOM_WSI_Transcoding_Claude.md); colour in [Color_Conversion_Claude.md](Color_Conversion_Claude.md))

Goal: transcode v1 iSyntax to DICOM WSI (YBR_RCT, JPEG 2000 lossless 5/3, sparse tiles, one SOP Instance per resolution layer, each layer self-contained) in the wavelet domain. Findings, all on `1.isyntax` and read-only (tools in `tools/v1_analysis/`):

- **Neither format stores overlap.** iSyntax is one seamless transform; a JPEG 2000 frame mirrors at its own edges. Interior coefficients are identical; only about 2 coefficients per edge, per level, per band differ (about 2 % of detail coefficients for a 256-px frame, 0.9 % for 512 px). **Verified.**
- **Alignment is required.** iSyntax `out_s[i]` sits at padded position 2^s (i+3) - 2, i.e. the J2K grid shifted by +2. The frame origin must satisfy `(i0+3) % 2^D == 0` (layer 0: valid x, y = 1 mod 2^D); a one-pixel shift breaks 38-99 % of coefficients. **Derived, verified.**
- **DICOM layer k = the decoder's pixel plane `out_k`;** its frames use iSyntax scales k and coarser as H bands. Levels beyond the top stored scale (5) are a plain local forward transform, so no layer reads another layer's output. **Verified for layers 0-5.**
- **Strip halo recomputation is exact:** copy interior H coefficients, recompute a ring of depth 2 (use 4 for margin) from mirror-extension strips, and a standard J2K inverse returns the iSyntax pixels bit for bit, for Y, Co and Cg, ten frames (`strip_halo.py`, `windowed_idwt_halo.py`).
- **A windowed top-down IDWT is enough:** it never reconstructs the whole image, only rectangles around the frame edges. (The scripts also run a full decode, but only as ground truth.) One 256-px layer-0 frame needs about 22 codeblocks per colour from 5 chunks, about 1 MB; the earlier conservative closure estimate (8-37 chunks, 3.2 MB median) was an upper bound.
- **Colour (measured, Color_Conversion_Claude.md C2; RGB target in C3):** Y' = Y, Cb = -Co/2 - Cg, Cr = Co/2 - Cg. The Philips SDK's RGB output (PAS TIFF levels 2, 3, 5) equals the truncating `/2` build in 100 % of 22 tiles and the floor build in 25-37 %, so the product's inverse is truncating (and has no lower clamp); the scanner's forward is then taken to be `fix` (inferred, not observable). Y' = Y is wrong by -1 in ~64 % of tissue pixels, Cb and Cr are exact. **Wavelet-domain conversion is biased, not ±1:** the error grows about -0.6 (Cb) and -0.4 (Cr) per decomposition level (D = 5: Cb mean -4.9, RGB off by up to 6-7 in ~100 % of pixels); a per-tile LL constant removes the mean but leaves std 1.3 (near-lossless at best). **Exact pixel-domain chroma destroys sparsity:** Cb level-1 detail coefficients go from 0.12 % to 63 % non-zero (0.02 to 1.4 bits/coefficient). Fidelity target is an open decision.
- **Chunks:** "chunk" is libisyntax's word; Philips says "cluster" (Hulsken 2016: for scanner streaming, tiles for viewer random access). A chunk is a storage unit, not a decode unit.
- **Tiles at 256, 512 and 1024 (`dicom_tile_tests.py`):** block-header coordinates are absolute level-0 positions (unit 0.25 um, origin x 77715, y 20614 units), sample i of layer k is at `ox + 2^k (i+3) - 2`, layer spacing 0.25 um x 2^k, tile origin `i0 = i_base + nF` with `(i_base+3) % 2^D == 0`. Per-tile exactness, standalone decode and 2x2 stitching pass for all three sizes and layers 0-5; the cross-layer offset check is exact. Blank (omittable) tiles exist only at layer 0 with F = 256 (457 of 8217) and 512 (49 of 2100), all pure white (255,0,0) and all in the right (about 421 px) and bottom (about 677 px) margin bands, i.e. synthetic filler; the 'blank glass' test regions are scanned glass (Y about 230) inside the envelope, where the exact colour route costs 3-4x. Edge tiles use black filler. PS3.3 was read (see the companion file, section 8b): regular tiling with a constant offset, overhanging frames with padding, absent tiles and the Recommended Absent Pixel CIELab Value, one pyramid per Series sharing a Frame of Reference, and YBR_RCT equations are all compatible with this plan; PS3.5 was read too (section 8c): `.90` lossless-only, raw codestream (no JP2 header), SGcod = 1 requires YBR_RCT, each frame encoded separately, decomposition levels and internal tiling unconstrained; **still open:** whether a frame may span several fragments, T.800 conventions checked only against my own NumPy transforms (no real codec yet), corner versus centre for the offsets, and the mapping of the iSyntax axes and origin to the slide system (Image Orientation (Slide)).
- **Quantizer/deadzone vs EBCOT (companion file 8e):** they empty out the finest scale (Y 1.55 % non-zero, chroma 0.01 %; chroma 96-99 % all-zero 64x64 blocks) and scale 1 of chroma; Y scales 1-5 and chroma scales 2-5 are unquantized and dense (51-69 % non-zero, 2.4-4 bits/coefficient) and hold nearly all bits. Zero-order total about 83 MB versus the 88.9 MB file; Cb+Cr cost about 10 % more than Co+Cg. T.800 Annex E (read): for truncated planes the reconstruction is q +- r 2^(Mb-Nb) with r chosen by the decoder (r = 0 is legal, r = 1/2 common, OpenJPEG effectively fixed at 1/2), so empty LSB planes can be dropped exactly only with an r = 0 decoder, only in code-blocks that do not touch a frame border (halo coefficients are arbitrary), and the gain is at most a few percent (companion file 8f); a stock decoder would be off by up to 4 RGB levels on 1024-px frames. 
- **Real EBCOT sizes (OpenJPEG 2.5.4 via Pillow 13.0, companion file 8g):** OpenJPEG's forward 5/3 matches our NumPy conventions (sparse-coefficient size test). On 2048x2048 real regions the wavelet-domain route codes to 0.87-0.99x of the iSyntax bytes (larger frames better: tissue 0.99x at 256 px, 0.93x at 1024, 0.92x at 2048); the exact pixel-domain route is 2.1-3.2x (floor) or 2.6-4.0x (fix); 32x32 code-blocks +3-4 %. The empty LSB planes cost at most 0.6-1.2 % of the stream, so dropping them is not worth it. A full six-layer DICOM pyramid is estimated at about 2.6x the iSyntax size (each layer carries the bits of scales k and coarser).
- **Omitting tiles (companion file 8h; every tile of every layer encoded with OpenJPEG; PNGs and per-tile CSVs in `$SOURCEROOT/transcoding_review/`):** the whole six-layer pyramid is about 204 / 199 / 196 MB for 256 / 512 / 1024 px frames (layer 0 alone 84.0 MB vs 88.9 MB iSyntax). Omitting only the exact-filler tiles saves almost nothing (0.08 MB at 256 px; a filler tile is about 170 bytes). Omitting near-constant scanned glass saves 23 % (robust range <= 4) to 34 % (<= 8) of the pyramid at 256 px frames (16-26 % at 512, 11-17 % at 1024), but is lossy and the robust rule also deletes dust specks (worst pixel up to 167 levels off); a strict no-speck rule (max-min <= 16) saves 13 %. One fill colour about RGB (228, 230, 233) fits the glass, not the white filler bands.
- **Not done:** real tile-by-tile pipeline output (DICOM file writing), the iSyntax-to-slide axis mapping, further sparse-frame rules, a real JPEG 2000 codestream (OpenJPEG), per-layer geometry metadata. See the companion file, section 8.

## Companion documents

- [iSyntax_to_DICOM_WSI_Transcoding_Claude.md](iSyntax_to_DICOM_WSI_Transcoding_Claude.md): wavelet-domain transcoding to DICOM WSI: halo geometry, strip algorithm, windowed IDWT, read footprint, spec reading, EBCOT sizes, tile omission, open items.
- [Color_Conversion_Claude.md](Color_Conversion_Claude.md): colour conversion for the DICOM transcoding (split out of the transcoding document): YCoCg-R to YBR_RCT (C1, C2: wavelet-domain bias versus exact but dense pixel domain) and to RGB (C3).
- [Pipeline_Plan_v2_Claude.md](Pipeline_Plan_v2_Claude.md): plan for the complete C++ iSyntax-to-DICOM pipeline (v2 supersedes [v1](Pipeline_Plan_v1_Claude.md): decisions 256 px / S = 16 / exact route / white absent frames / lossy tag always set, and the built OpenJPEG coefficient-input patch). Versioned: each iteration is a new `Pipeline_Plan_v<N>_Claude.md`; use the highest N.
- OpenJPEG coefficient-input patch: `patches/openjpeg-2.5.4-coefficient-input.patch`, built and tested with `tools/transcode_tests/build_patched_openjpeg.sh` and `ext/test-bin/coeff_input_test` (`ext/` is git-ignored).

- [Multi_File_Strip_Analysis_Claude.md](Multi_File_Strip_Analysis_Claude.md): the strip-wise (bounded-memory) reconstruction, verified bit-exact on `1.isyntax`, and its results on six v1 files up to 12.8 Gpx (memory trend 1.8 to 3.4 GB versus 16.9 GB for the whole-image script on the smallest file), the quantizer/deadzone check on all six, per-slide colour-rounding findings, a Python-port fix (zero-run symbol absent from the Huffman tree) and the single-core OpenJPEG Tier-1 timing.
- [Example_iSyntax_File_Claude.md](Example_iSyntax_File_Claude.md): analysis of the one real file: byte layout and chunk order, bitplane statistics, the Quantizer/Deadzone check on every codeblock, and an independent full-image reconstruction (clean, no tile seams).
- [Example_FIC_Files_Claude.md](Example_FIC_Files_Claude.md): analysis of the six sample headers, including quantizer/deadzone evidence.
- [Philips_Compression_Code_Samples_Claude.md](Philips_Compression_Code_Samples_Claude.md): review of Philips' MATLAB encoder/decoder (`/mnt/d/Whole_Slide/iSyntax_Reference_Implementation/CompressionCodeSamples`, outside the repo; all three Philips sample folders are now reviewed). It is the bit-level reference for Hulsken, quantization and the wavelet geometry.
- [Image_Format_Code_Samples_Claude.md](Image_Format_Code_Samples_Claude.md): review of Philips' Python and MATLAB for reading the file (`/mnt/d/Whole_Slide/iSyntax_Reference_Implementation/ImageFormatCodeSamples`). It covers the v1 layout only and cannot read the v2 sample headers past the metadata step. The C code's block-ID logic is a port of it.
- [Viewing_Code_Samples_Claude.md](Viewing_Code_Samples_Claude.md): review of Philips' CLAHE and sharpening reference (display post-processing; not decoding). libisyntax does none of it. All its constants are hard-coded: 40 px regions, 256 bins, clip limit 7 per bin, sharpening gain 2.0 halved per level.
- [Third_Party_Documentation.md](Third_Party_Documentation.md): review of the Philips specs, papers and earlier AI notes, with discrepancies and gaps. Read its section 5 for known issues in this code (swapped arguments in `isyntax_is_parent_tile_missing`, scalar YCoCg path rounding and clamping).

## Conventions

- C11, 4-space/tab mix as in the existing files. Keep upstream code style and the original copyright headers (Valkema, OpenJPEG BSD-2, stb, yxml).
- When adding analysis code (dump bitplanes, histogram bitmasks, print FIC quantizer values), put it in a **new example/tool file** and leave the decoder untouched. `isyntax_hulsken_decompress` has commented-out `dump_block` hooks you can mimic.
- Debug aids already present: `ISYNTAX_WANT_DEBUG_OUTPUT_PNG` (writes IDWT steps as PNG) and `DEBUG_OUTPUT_IDWT_STEPS_AS_PNG`.
- Mark claims as **verified in code**, **from spec**, or **hypothesis**. Cite `file:line` and the spec section.
- Source comments in `src/isyntax/` follow the same rule: facts measured on the real file are tagged as such, and nothing in them changes behaviour. When you add comments, re-run the identity check: `git show HEAD:src/isyntax/isyntax.c | gcc -fpreprocessed -dD -E -P -x c - | sed '/^\s*$/d' | cmp - <(gcc -fpreprocessed -dD -E -P -x c src/isyntax/isyntax.c | sed '/^\s*$/d')` (comments-only changes print nothing).
