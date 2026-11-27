# Multi-file analysis with the strip-wise reconstruction (v1 files, 2026-11-26)

Companion to [Example_iSyntax_File_Claude.md](Example_iSyntax_File_Claude.md) (the one-file analysis of `1.isyntax`). Tags: **[verified]** (measured here), **[hypothesis]**, **[open]**.

## 1. Why a strip-wise reconstruction

`tools/v1_analysis/reconstruct_np.py` keeps every coefficient array, several int32 intermediates and three full int16 planes in memory: about 31 bytes per padded pixel. **[verified]** 16.85 GB max RSS (195 s) for `1.isyntax` (0.55 Gpx). The eleven additional v1 files (all data model 5.0, Quality=2, compressor 16) reach 7 to 13 Gpx, which would need 200 to 400 GB.

`tools/v1_analysis/reconstruct_strips.py` produces the same findings in bounded memory:

- It works on full-width bands of rows (default target 64 MB per int32 plane), top-down from the coarsest scale.
- For exact output rows `[A,B)` at scale 0 it needs, at scale `s`, the coefficient-row window `n0 = floor((A-1)/2)`, `n1 = ceil((B+2)/2)`, and the next coarser scale must then deliver `out` rows `[n0-1, n1-1)`. The 5/3 inverse on a window with mirrored ends is exact for output rows `[2 n0 + 1, 2 n1 - 2)`; windows clipped at the true image border behave like the whole-image run.
- Codeblock rows (128 coefficient rows) are decoded on demand by a fork pool, cached, and evicted once the band has moved past them.
- The Y/Co/Cg to RGB statistics, the seam profile and the preview are accumulated band by band.

**[verified]** Bit-exact: on `1.isyntax`, with `VERIFY_PLANES` pointing at the planes saved by `reconstruct_np.py`, 0 of 1 598 519 916 Y/Co/Cg samples differ. Findings are identical to the old script (Y 48 to 257, 64.5 % floor-versus-truncate difference, no negative Y).

**[verified]** Memory (summed Pss over the process tree, sampled every 0.5 s): 1.78 GB for `1.isyntax`, against 16.85 GB for the old script. (A run with `VERIFY_PLANES` shows 4.68 GB because the memory-mapped reference planes are counted; use the clean figure.)

Implementation notes: Python 3.14 defaults to the forkserver start method, which needs an AF_UNIX socket that the sandbox forbids, so the script uses `get_context('fork')`. The pool is created before the monitor thread starts.

## 2. Files and results

All twelve files are v1 (interface version 5.0) with derivation string `PHILIPS UFS V1.8.6614 | Quality=2 | DWT=1 | Compressor=16`; the v1 header has no wavelet-quantizer elements, so the Philips Q2 table of `exhaustive.py` applies. The block tables of all twelve were parsed (`parse_v1.py`; the XML header of the larger files exceeds 8 MB, so the parser now reads until the EOT marker). Six files were analysed completely, in order of increasing size, until the memory trend repeated (stop rule agreed with the owner).

| File | Valid size (px) | Mpx | Levels | Non-empty blocks | Compressed codeblocks | Peak memory | Wall time |
|---|---|---|---|---|---|---|---|
| 1 | 25222 x 21126 | 533 | 6 | 17 626 | 84 MB | 1.78 GB | 88 s |
| 3 | 112134 x 65030 | 7 292 | 8 | 94 408 | 475 MB | 2.90 GB | 18.4 min |
| 11 | 117254 x 60934 | 7 145 | 8 | 109 264 | 612 MB | 3.22 GB | 18.6 min |
| 5 | 121350 x 66054 | 8 016 | 8 | 77 747 | 439 MB | 3.05 GB | 20.4 min |
| 12 | 207878 x 32774 | 6 813 | 9 | 107 593 | 510 MB | 3.32 GB | 19.4 min |
| 4 | 162822 x 78854 | 12 839 | 9 | 119 572 | 599 MB | 3.43 GB | 32.9 min |

Not analysed (block tables parsed only): 6, 7, 8, 9, 10, 13. There is no `2.isyntax`.

**[verified] Memory trend.** Peak memory is bounded: 24 times more area than file 1 costs 1.9 times more memory, and the five 6.8 to 12.8 Gpx slides lie between 2.90 and 3.43 GB. It does not follow image width alone (the widest slide, file 12, is not the highest). A working model **[hypothesis]** is about 2 GB fixed plus about 1 GB for the cached block rows, which grows slowly with width and level count. Wall time grows with area (about 2.5 min per Gpx) because the main thread does all IDWT and statistics; the decoder pool is mostly idle.

**[verified] Quantizer/Deadzone check (`exhaustive.py`) on files 1, 3, 4, 5, 11, 12:** every codeblock decodes, zero forced-clear mask-bit violations, zero length mismatches (after the decoder fix below). The finding of `Example_iSyntax_File_Claude.md` therefore holds on six slides. Stream statistics (`stream_stats.py`): compression ratio raw/compressed 17.4 to 20.7 to 1 on all of them.

**[verified] Colour-conversion findings per slide.**

| File | Y range | Floor vs truncating differ | RGB < 0 before clamp | RGB > 255 before clamp |
|---|---|---|---|---|
| 1 | 48 to 257 | 64.5 % | 0 | 1.89 M |
| 3 | 44 to 258 | 23.5 % | 0 | 23.5 M |
| 11 | 46 to 258 | 26.7 % | 0 | 3.1 M |
| 5 | 0 to 258 | 18.4 % | 0 | 4.7 M |
| 12 | 37 to 258 | 25.7 % | **42** | 11.8 M |
| 4 | 0 to 258 | 17.4 % | **19** | 12.3 M |

- Y is never negative on any slide, so the absolute-value step after the IDWT changes nothing. The differences between the roundings stay within R/B 0 to 2 and G 0 to 1 levels.
- Files 12 and 4 are the first slides with pixels below 0 before clamping (42 and 19 pixels). The scalar `ycocg_to_rgb` has no lower clamp, so on the Meson default build those pixels would wrap; this was not observable on `1.isyntax`.
- The share of pixels affected by the rounding depends on the slide (tissue-rich 1.isyntax 64.5 %, glass-heavy slides 17 to 27 %).
- Seam profile: max / median of the mean |dY| by position modulo 256 is 1.08 to 1.53 on all slides (1.30 on `1.isyntax`, documented as clean), and the position of the maximum differs per file, so there is no consistent tile-seam position.

## 3. Decoder issue found: zero-run symbol absent from the Huffman tree

`exhaustive.py` failed on two top-level Y LL codeblocks of `3.isyntax` with `KeyError` in the Python port (`hulsken_v1.py`): the block's Huffman tree contains no leaf for the stored zero-run symbol, i.e. the block never uses zero runs. `isyntax.c` handles this implicitly (`huffman.code[zerorun_symbol]` stays 0, size forced to 1, and the zero-run branch is never taken). The port now does the same (`codes.get(zr_sym, (0, 0))`). **[verified]** After the fix all blocks of `3.isyntax` decode with no violations. The C decoder is unaffected.

## 4. Tier-1 timing of OpenJPEG on real frames (`tier1_profile.py`)

Stock OpenJPEG Tier-1 is the code the patched coefficient-input encoder runs (Pipeline_Plan_v2 section 2.1). Single core, `1.isyntax`, 256-px frames of the wavelet-domain colour route, best of five repeats; the cost of DWT + colour transform + headers is measured on constant frames (0.7 ms) and removed. The machine was not idle (the strip runs were using most cores), so the figures may be pessimistic.

| Layer | Frames | Median frame bytes | Median encode | Tier-1+2 speed |
|---|---|---|---|---|
| 0 | 8 217 | 6 975 | 2.25 ms | about 6.0 MB/s coded per core |
| 1 | 2 100 | 24 446 | 4.46 ms | about 6.6 |
| 2 | 525 | 53 654 | 9.40 ms | about 7.5 |
| 3 to 5 | 197 | 41 to 50 KB | 7 to 8 ms | about 7.2 to 7.5 |

**[verified]** Whole six-layer pyramid of `1.isyntax` (11 039 frames, 204.4 MB coded): about 34 core-seconds for Tier-1 + Tier-2, about 3.1 ms per frame; the exact (coefficient-input) route avoids about 18 % of the stock encode time. **[hypothesis]** Scaling by coded size (pyramid about 2.3 times the iSyntax file) gives about 8 core-minutes for a 1.4 GB file, about 30 to 40 s on 16 cores; Tier-1 is therefore not the bottleneck of a C++ pipeline, and I/O and Hulsken decoding are more likely to dominate. The fits for layers 3 to 5 rest on 12 to 26 frames each. The C Hulsken decode speed is **[open]**.

## 5. Related code added in this branch

- `libisyntax_tile_read_coefficients()` and helpers (`libisyntax.h`, `isyntax_reader.c`): stored H (any level) and LL (coarsest level only) coefficients of a tile, lock-free and thread-safe, no pixels. Example program `isyntax_coefficients` (serial-versus-parallel self-check). Compiled with Meson but **not run** (Claude sessions cannot execute compiled binaries); the owner should run `isyntax_coefficients FILE 0 all 8`.
- Per-file results: `$SOURCEROOT/c_decoder_output/strip_analysis/<file>/` (blocks.csv, exhaustive.log, stream_stats.log, strips/out.log, strips/time.log; file 1 also `old_reconstruct/`, `tier1_profile.log`).

## 6. Reproduce

```bash
.venv/bin/python3 -B -I tools/v1_analysis/parse_v1.py ../example_isyntax_files/N.isyntax BLOCKS.csv
.venv/bin/python3 -B -I tools/v1_analysis/reconstruct_strips.py ../example_isyntax_files/N.isyntax BLOCKS.csv OUTDIR 64   # env JOBS=12 BAND_MB=64
.venv/bin/python3 -B -I tools/v1_analysis/exhaustive.py ../example_isyntax_files/N.isyntax BLOCKS.csv
.venv/bin/python3 -B -I tools/v1_analysis/stream_stats.py ../example_isyntax_files/N.isyntax BLOCKS.csv
.venv/bin/python3 -B -I tools/v1_analysis/tier1_profile.py FILE BLOCKS.csv ../transcoding_review OUT.csv   # per-frame sizes come from tile_stats_F256_L*.csv
```
