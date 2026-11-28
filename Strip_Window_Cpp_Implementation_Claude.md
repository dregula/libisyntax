# C++ strip-window transcoder and DICOM layer writer (2026-11-27)

Written for the project owner. Tags: **[verified]** (measured here), **[hypothesis]**, **[open]**. Related: [Halo_Strategy_Comparison_Claude.md](Halo_Strategy_Comparison_Claude.md) (why the strip window), [Pipeline_Plan_v2_Claude.md](Pipeline_Plan_v2_Claude.md), [iSyntax_to_DICOM_WSI_Transcoding_Claude.md](iSyntax_to_DICOM_WSI_Transcoding_Claude.md), [Multi_File_Strip_Analysis_Claude.md](Multi_File_Strip_Analysis_Claude.md).

## 1. What was built

New files only; the decoder is untouched. All of it uses the public `libisyntax.h` API (which now has `extern "C"` guards) plus the lock-free `libisyntax_tile_read_coefficients`.

| File | Role |
|---|---|
| `src/transcode/strip_transcoder.cpp` | Strip-window transcoder (executable `isyntax_strip_transcode`). |
| `src/transcode/dicom_layer_writer.{h,cpp}` | One DICOM file per layer: DCMTK 3.6.9 metadata, hand-written encapsulated Pixel Data, Extended Offset Table patched in afterwards. |
| `src/transcode/srgb_lcms.icc` | sRGB profile (588 bytes, made with Pillow/lcms) for the optical path. |
| `tools/transcode_tests/build_strip_transcoder.sh` | Build against `build-ssse3` (baseline build), the patched OpenJPEG in `ext/`, and `libdcmtk-dev`. |
| `tools/v1_analysis/dump_frame_coeffs.py` | Python reference frames for cross-checking. |

**Algorithm (per layer k, 256-px frames, D = 5, ring depth h = 4):** frames are processed in raster rows over one sliding window of decoded codeblocks. A planner (indices only) lists the blocks each frame touches; before a row runs, the window keeps the blocks of the current and next row, drops the rest, and decodes the missing ones in parallel. Frame tasks of a row then run in parallel over the read-only window. Per frame (port of `windowed_idwt_halo.py`): interior H coefficients are copied from the stored blocks; the ring is recomputed from thin strips of a windowed top-down IDWT; frames at the array border (0 to 9 % per layer) use the local transform of the zero-padded decoded pixels. Then wavelet-domain colour conversion (Y' = Y, Cb = -q(Co) - Cg, Cr = Cb + Co, q = truncating halving, -128 on LL of Y'), Mallat layout, and Tier-1 through the patched OpenJPEG (`opj_set_coefficient_input`, in-memory streams, one codec per frame, 64x64 code-blocks, 6 resolutions).

**Omission and constant frames (owner decisions 2026-11-27):** on layers 0 and 1 a frame is omitted when the strict range of its RGB (truncating inverse, clamped) is <= S = 16. **Omitted frames are not left absent: they are written as colour-constant frames** (166 bytes for white), so the files are dense (TILED_FULL) and a viewer needs no absent-frame convention. The colour is the frame's mean RGB (constant frames from the shortcut: their exact value). A coefficient prefilter (maximum |H| of the frame's own scales <= T = 3 S = 48) avoids the pixel decode of obvious tissue. **Constant-frame shortcut (exact):** if every block of the frame's footprint is empty, the top LL is constant over the footprint, and no window reached index 0 (the zero shift L'[0] = 0), the frame is constant; it is built without any lifting, and in the omission test it needs no pixel decode.

## 2. Verification **[verified]**

| Check | Result |
|---|---|
| Builder vs full per-frame forward transform of the decoded pixels (T.800 5/3), all colours, sampled interior frames | 0 mismatches: 4 422 frames on `1.isyntax`, 17 420 on `5`, 17 024 on `3` (study without omission) plus the samples of every production run |
| Constant-frame shortcut vs the full builder | 0 mismatches: 8 381 frames (5), 5 740 (3) |
| C++ vs the Python implementation (`dump_frame_coeffs.py`), floor colour, byte for byte | 10 of 10 identical frames of `1.isyntax` (layers 0 to 5, tissue, glass, scan edge, border frames) |
| Codestream bytes vs the Python/Pillow route (`tile_stats_F256_L*.csv`) | 4 466 of 4 487 frames without out-of-range pixels are identical; **all** 2 687 differing frames have pixels outside 0..255 after decoding our coefficients (the Python route clips to 8 bit, ours keeps the coefficients); no differing frame is free of clipping |
| Tier-1 round trip: stock decoder vs the inverse of our own coefficients (catches range overflow in Tier-1) | 0 bad of 286 (file 1, floor), 17 (1, trunc), 102 (5), 125 (3) sampled frames |
| Omission audit: exact strict-range test on every layer-0/1 interior frame vs the prefilter | 3 865 omittable frames on `1`, 123 226 on `5`, 99 787 on `3`; **0 rejected by the prefilter at T = 48**. Largest |H| seen in an omittable frame: 17/16/15 (file 1), **32/32/16 (file 5)**, 17/20/16 (file 3) for Y/Co/Cg. The first default T = S + 4 = 20 wrongly rejected 9 frames on file 5; the bound is therefore about 2 S, not S + 1 **[empirical, three files; no proof]** |
| `dciodvfy` (dicom3tools) on the DICOM files | after adding Total Pixel Matrix Focal Planes, empty Acquisition Context and the two Issuer sequences and removing Number/Distance of Focal Planes: no errors; one warning (`PHILIPS_DP_1_0` is not a defined term of Lossy Image Compression Method; kept on purpose to record the source). The OV Extended Offset Table is skipped by this older tool |
| wsidicom 0.37 / pydicom 3.0.2 open the directory of six files | **one pyramid of 6 levels** (25 413 x 21 317 ... 795 x 667, mpp 0.25), thumbnail and region reads work. Before the fix only 3 levels were found: wsidicom groups instances into pyramids by image origin (threshold 0.01 mm) and my per-layer frame grids put the layer origins up to 0.25 mm apart |
| Stitched layers 4 and 5 from the DICOM files viewed | tissue, glass, white margins as expected (black padding in the first 31 pixels at top and left, see section 5) |

## 3. Results (16 threads, `build-ssse3`, trunc colour, Tier-1 on, omitted frames as constants)

| File | Valid size | Frames (6 layers) | Omitted (of which no pixel decode) | Fast-path frames when omission is off | Coded output | Wall | Peak RSS | Peak window |
|---|---|---|---|---|---|---|---|---|
| `1.isyntax` | 0.5 Gpx | 11 039 | 3 865 (0) | 0 | 177 MB | 9.2 s | 0.22 GB | 406 blocks, 120 MB |
| `5.isyntax` | 8.0 Gpx | 164 377 | 123 226 (54 934) | 58 544 (36 %) | 946 MB | 81 s | 0.71 GB | 1 910 blocks, 563 MB |
| `3.isyntax` | 7.3 Gpx | 149 526 | 99 787 (37 750) | 40 177 (27 %) | 1 049 MB | 92 s | 0.67 GB | 1 768 blocks, 521 MB |

- The shortcut cut the omission pixel-decode time from 770 to 422 thread-seconds and the wall time from 108 s to 81 s on `5` (92 vs 103 s on `3`). The omission pixel test is still the largest single CPU item for the glass-heavy files; Tier-1 is second. `1.isyntax` is tissue-rich and has no constant frames: its glass carries noise in every scale.
- Block decodes equal the per-layer unique count of the footprint simulation (each block once per layer, 1.6 x the file).
- Memory is bounded by the window, not by the image: 0.2 to 0.7 GB for 0.5 to 8 Gpx.

**Common origin (found with wsidicom).** The J2K alignment `(i0 + 3) % 32 == 0` puts the possible origin of layer k on a grid of 2^(k+5) level-0 units, so a free choice per layer scatters the layer origins by up to 0.25 mm. All layers now start on the common grid of the coarsest layer (`i_base_k = 32 * 2^(top-k) * M - 3`, same M for all layers); the remaining origin spread is the sample-centre offset of the wavelet pyramid (at most 2^(top-1) units = 4 um at top = 5), below wsidicom's threshold. Cost: the matrix starts at the padded origin, so layer 0 carries about 190 extra pixels of padding at top and left (100 x 84 frames for file 1). **Padding is now white** (Y = 255, Co = Cg = 0, as the iSyntax margin filler) instead of black, and frames wholly outside the decoded array are constant white frames.

## 4. DICOM files **[verified structurally; semantics not validated]**

`isyntax_strip_transcode FILE --dicom OUTDIR [--orientation a,b,c,d,e,f] [--icc file]` writes `<name>_L0.dcm` ... `<name>_L5.dcm` (one SOP instance per layer, shared Study, Series, Frame of Reference and Specimen UIDs; UIDs `2.25.<128 bit>` derived deterministically from the file name). `1.isyntax` gives 62.3, 64.6, 37.3, 9.9, 2.6 and 0.7 MB (177 MB), identical over three runs. Label and macro images are skipped (owner instruction).

- VL Whole Slide Microscopy Image Storage, transfer syntax JPEG 2000 Lossless Only, YBR_RCT, 8 bit, 256 x 256 frames, **TILED_FULL** with the two Dimension Index pointers, shared functional groups (pixel measures: 0.00025 mm x 2^k, frame type, optical path identification), one optical path with the sRGB profile and brightfield/full-spectrum codes, specimen module with dummy identifiers, Lossy Image Compression = 01 with method `PHILIPS_DP_1_0\ISO_15444_1` and ratios (source, ours), Total Pixel Matrix (layer 0: 25 253 x 21 157), matrix origin from the iSyntax origin, Imaged Volume.
- Image Type is `DERIVED\PRIMARY\VOLUME\NONE` (layer 0) and `...\RESAMPLED` (layers 1 to 5). Patient, study and equipment values are dummy or taken from the iSyntax header (manufacturer, model, serial, acquisition time).
- Pixel Data: empty Basic Offset Table, one fragment per frame (padded to even length), sequence delimiter; Extended Offset Table and Lengths patched in place. `dcmftest` and `dcmdump` read all files; the offsets and lengths match the real item positions; Pillow/OpenJPEG decodes frames of every layer.

## 5. Open items and cautions **[open]**

1. **Slide orientation and origin.** Image Orientation (Slide) defaults to (0,-1,0, -1,0,0). The identity (1,0,0, 0,1,0) is **not a valid slide orientation** (wsidicom warns that the column direction must be the row direction rotated by -90 degrees), so a right-handed one is used. The matrix origin is the iSyntax corner of the layer-0 matrix placed in the lower-left of the slide bounding box for the chosen orientation. This is a placeholder: the owner will establish how SlideScape orients the slide; use `--orientation` and, if needed, a change of the origin mapping.
2. **Validators:** dciodvfy and wsidicom now run clean (section 2). Still unchecked: the Extended Offset Table lengths convention (padded item length used, from memory of CP 2388; dciodvfy skips OV) and OmniHisto behaviour (multi-file pyramid, TILED_FULL, YBR_RCT, ICC).
3. ~~Black padding~~ is now white (see the common-origin paragraph).
4. **Colour:** wavelet-domain conversion is biased (Color_Conversion_Claude.md): expect RGB differences of several levels against SlideScape, and a generic sRGB profile that may differ from what SlideScape applies.
5. Prefilter bound: empirical (about 2 S); keep the audit (`--audit`) in any regression run. The omission pixel test remains the largest CPU item on glass-heavy slides.
6. The coefficient-range guard of Pipeline_Plan_v2 2.1 is only checked by sampled round trips (`--roundtrip`); it is not enforced per frame.
7. Not done: label/macro images (skipped by decision), per-frame functional groups (not needed for TILED_FULL), DCMTK-based reading test, C++ unit tests, packaging.

## 6. A bug worth remembering

The first thread pool reset its "active workers" counter at the start of every job; a worker that woke late for the previous job could then end a job early or deadlock the next one. Symptoms: a missing codestream on layer 5 (about one run in several) and two study runs hanging for hours. The pool now gives every job its own counters; three repeated runs of file 1 and the reruns of files 5 and 3 are clean. Results from before the fix were checked against the grid size (frames processed = frames in the grid) and were complete.

## 7. Run

```bash
ninja -C build-ssse3 && tools/transcode_tests/build_strip_transcoder.sh            # needs ext/ (build_patched_openjpeg.sh) and libdcmtk-dev
build-ssse3/isyntax_strip_transcode ../example_isyntax_files/1.isyntax --threads 16 --dicom ../c_decoder_output/dicom_out/1
#   verification flags: --verify-every K --roundtrip --audit --oob --csv out.csv --omit -1 (no omission, fast-path study) --dump k ty tx dir
```
Logs of all runs: `$SOURCEROOT/c_decoder_output/strip_analysis/<file>/cpp_*.log`.
