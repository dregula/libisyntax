# Halo correction for the iSyntax to DICOM transcoder: strip versus recruiting individual tiles (2026-11-27)

Written for the project owner, who must transcode about a million iSyntax files. Tags: **[verified]** (measured or simulated from the real block tables), **[hypothesis]**, **[open]**. Companion documents: [iSyntax_to_DICOM_WSI_Transcoding_Claude.md](iSyntax_to_DICOM_WSI_Transcoding_Claude.md) (halo geometry, windowed IDWT), [Pipeline_Plan_v2_Claude.md](Pipeline_Plan_v2_Claude.md) (pipeline), [Multi_File_Strip_Analysis_Claude.md](Multi_File_Strip_Analysis_Claude.md) (Tier-1 timing, large files).

## 1. What is compared

A DICOM frame (JPEG 2000, 256 px, D = 5) of layer k needs the J2K coefficients of the frame. Interior coefficients are copied from stored H blocks; a ring of depth h = 4 along the frame edges differs, because iSyntax is one seamless transform and J2K mirrors at each frame edge. The ring is recomputed from thin pixel strips, which come from a windowed top-down IDWT. Two ways to organise this work were compared (my reading of the request; **[open]** if "recruiting tiles" meant something else):

| Name | Meaning |
|---|---|
| **STRIP** | The planned pipeline: windowed IDWT of thin strips (width 2h + 4) around frame edges, shared between neighbouring frames, with a sliding window of decoded block rows (every block decoded once per layer, or once overall when the six layers advance together). |
| **THIN** | Same algorithm, but every frame is an independent task: it fetches (and decodes) its own footprint of blocks and computes its own four rings. With a shared LRU block cache it converges on STRIP. |
| **RECRUIT** | Halo correction by recruiting whole neighbouring iSyntax tiles the way `libisyntax_tile_read` does: per frame, at every scale from the frame's own up to the top, the tiles under the frame footprint (ancestors above the frame's scales) plus their 8 neighbours are loaded, and every footprint tile is inverse-transformed in full (2 x 136)^2 samples per colour and tile. |

## 2. Method

- `tools/v1_analysis/halo_strategy_footprint.py` replays the exact window arithmetic of `windowed_idwt_halo.py` (`build_windowed`, class `Windowed`) with indices only, over every interior frame of every layer (F = 256, D = 5, h = 4), and records the touched codeblocks (only non-empty ones, size > 8 bytes, cost anything), the compressed bytes, the lifted output samples and the strip forward-transform samples. **Check [verified]:** it reproduces the footprints measured by the real decoder for the tissue frames of `1.isyntax`: 22 / 12 / 8 blocks per colour at layers 0 / 3 / 5. Edge frames (partly outside the decoded array: 0 to 9 % of frames per layer, mostly at the coarse layers) are skipped.
- Cache effect: frames in raster order through an LRU cache of C blocks (a block = the three colours' codeblocks at one position).
- Transform sharing: with `UNION=1` the rectangles of all `out_region` calls of a layer are merged (each IDWT output cell computed once), a lower bound for shared strips. Done for `1.isyntax` only (the difference arrays are too large for the 7 to 13 Gpx files).
- Time constants, all measured in C on `1.isyntax` with `build-ssse3` (Linux, Meson, `-mssse3`), one core (`tools/transcode_tests/bench_decode.c`, `bench_idwt.c`; build with `gcc -O2 -Isrc X.c -Lbuild-ssse3 -lisyntax -Wl,-rpath,$PWD/build-ssse3`):
  - Hulsken decode (the new `libisyntax_tile_read_coefficients`, page-cache warm): all 17 590 non-empty H blocks of the file in 1.18 s. Least squares over the six levels: **t = 30.9 us per block + 7.7 ms per MB compressed (130 MB/s)**, within 5 % per level.
  - IDWT: 800 cold level-0 pixel tiles in 0.93 s; subtracting the decode share gives about **3 ns per lifted sample** (includes RGB conversion, so slightly pessimistic).
  - Tier-1 + Tier-2 of OpenJPEG: 6.0 to 7.5 MB/s coded per core (`tier1_profile.py`, noisy); coded pyramid taken as 2.3 x the file **[hypothesis; measured on `1.isyntax` only]**.
- `halo_strategy_costmodel.py` turns the footprint logs into core-seconds.

## 3. Results **[verified]** (simulation from real block tables; times use the constants above)

Totals over all interior frames and layers; "x shared" is relative to decoding every non-empty block once (= the file).

| | `1.isyntax` (89 MB, 0.5 Gpx, 10 827 frames) | `5.isyntax` (483 MB, 8.0 Gpx, 164 k frames) | `4.isyntax` (680 MB, 12.8 Gpx, 262 k frames) |
|---|---|---|---|
| STRIP, every block once | 17 626 decodes, 84 MB | 77 747, 439 MB | 119 572, 599 MB |
| STRIP, each layer its own window | 1.60 x, 215 MB | 1.62 x, 1.1 GB | 1.60 x, 1.5 GB |
| THIN, per-frame, no cache | 18.7 x, 4.5 GB read (53 x file) | 31.2 x, 34 GB (77 x) | 31.8 x, 49 GB (82 x) |
| THIN + LRU 256 blocks | 3.4 x | 4.6 x | 4.5 x |
| THIN + LRU 1024 blocks | **1.60 x** | 3.65 x | 4.0 x |
| THIN + LRU 4096 blocks | 1.60 x | **1.62 x** | **1.60 x** |
| RECRUIT, no cache | 103 x, 26 GB (315 x) | 184 x, 209 GB (477 x) | 172 x, 276 GB (461 x) |
| RECRUIT + LRU 4096 blocks | 1.60 x | 1.62 x | 7.2 x |
| Lifted samples per colour: THIN / RECRUIT | 0.84 G / 8.7 G | 14.2 G / 158 G | 23.8 G / 274 G |

Findings:

1. **Tile-independent frame tasks re-decode each block 19 to 32 times** (THIN) and 100 to 180 times (RECRUIT) without a shared cache; read volume is 53 to 82 times and 315 to 477 times the file.
2. **A shared cache removes it, but its size must follow the image width.** To reach the strip's 1.6 x the cache must hold the block rows the ring windows touch: 1 024 blocks suffice for the 25 k px wide file, 4 096 for the 120 to 166 k px wide ones (about 3 block rows at each scale, about 7 KB per pixel of width at 3 colours x 98 KB per block, i.e. about 1.2 GB at 166 k px). That cache **is** the strip's sliding window; a per-frame task pool with a big enough shared cache and raster-order scheduling is STRIP in another form. RECRUIT needs a still larger cache (7.2 x at 4 096 blocks on the widest file).
3. **The 1.6 x for "each layer its own window"** arises because a block of scale s serves layers s-4..s. Running the six layers in lockstep gives 1.0 x, but the saving is about 1 to 6 core-seconds per file, negligible next to Tier-1.
4. **Shared strips do much less lifting.** On `1.isyntax` the union of all requested IDWT windows is 293 M samples per colour against 840 M summed over frames: **65 % saved** (66 % at layer 0, falling to 41 % at layer 4). Adjacent frames request overlapping windows, especially at the coarse scales. RECRUIT lifts 10 to 12 times more than THIN, because whole tiles with 4-sample margins are transformed for every frame. Pixel-bounce (full IDWT then full forward DWT) is not modelled.

### Core-seconds per file **[mixed: decode and transform from the model above; Tier-1 and the 65 % sharing transferred from `1.isyntax`, so hypothesis for files 5 and 4]**

| Strategy | `1.isyntax` | `5.isyntax` | `4.isyntax` |
|---|---|---|---|
| Tier-1 + 2 (same for all) | 31 | 171 | 241 |
| STRIP (per-layer windows, shared cells): decode + transform | 2.5 + 4.9 | 12.5 + 79 | 17.6 + 130 |
| **STRIP total (model 0.44 to 0.57 core-s per MB of file; measured end to end later: 1.1 to 1.4, see Strip_Window_Cpp_Implementation_Claude.md and DpsBulk_Collection_Claude.md)** | **39 (1.00)** | **263 (1.00)** | **388 (1.00)** |
| THIN + sufficient cache | 44 (1.13 x) | 346 (1.32 x) | 528 (1.36 x) |
| THIN, no cache | 86 (2.2 x) | 669 (2.6 x) | 1006 (2.6 x) |
| RECRUIT, no cache | 379 (9.7 x) | 3781 (14 x) | 5674 (15 x) |

Where the CPU goes (STRIP): Tier-1 about 80 % for the tissue-rich file, about 62 to 65 % for the glass-heavy large ones, where the transforms (lifting plus ring forward transforms, proportional to image **area**) become a comparable second cost. Decode is only 3 to 7 % once it is done once and in C (130 MB/s per core).

## 4. Fleet view for a million files **[hypothesis; arithmetic on the above]**

Mean size of the 12 example files is 840 MB (a guess for the population). STRIP at about 0.5 core-s per MB (the model; the measured transcoder costs 1.1 to 1.4, so multiply by about 2.5) gives about 420 core-s per file, **about 4.2 x 10^8 core-seconds = 117 000 core-hours** (about 5 days on 1 000 cores, 1 day on 5 000 cores). Pyramid output is about 2.3 x the input (measured on one file): about 1.9 PB for 1 M files of 840 MB; storage and write bandwidth may bind before the CPUs do.

| Strategy | CPU total (x STRIP) | Bytes read per file | Notes |
|---|---|---|---|
| STRIP, one single-threaded worker per file | 1.0 | 1 x (sequential), 1.6 x per-layer windows | memory about 7 KB x image width (about 0.2 GB at 25 k px, 1.2 GB at 166 k px, C++, int16), no locks, no shared state |
| THIN + shared cache | 1.1 to 1.4 | 1.6 to 2.6 x | cache contention across threads of one file; equals STRIP if scheduled in raster order |
| THIN, no cache | 2.2 to 2.6 | 53 to 82 x | acceptable only if the file is staged on local NVMe or in RAM; from network or object storage it multiplies traffic 50 to 80 times |
| RECRUIT, no cache | 10 to 15 | 315 to 477 x | prototype only |

Parallelism: with a million files the parallelism comes from the files, so STRIP's sequential dependence inside a file costs nothing; the largest file here (1.7 GB) is about 15 core-minutes. Intra-file parallelism (frame tasks) matters only for tail latency, partial jobs, or fewer, larger inputs.

## 5. Contrast

| | STRIP (sliding window) | Independent tiles (THIN / RECRUIT) |
|---|---|---|
| CPU per file | lowest | +10 to 40 % with a right-sized cache, 2.2 to 2.6 x without (THIN), 10 to 15 x (RECRUIT) |
| I/O pattern | sequential, once | random, repeated; needs local staging |
| Shared IDWT work | yes (about 65 % less lifting) | no |
| Memory | bounded by width (about 7 KB per px of width) | small per task, but a shared cache of the same size is needed to be efficient |
| Parallelism | across files; within a file only per stage | trivially across frames of one file |
| Restart / retry | row granularity; re-decode one window | per frame |
| Random access / partial output (a region, a few frames, streaming from object storage) | poor (must prime the window) | good; the only sensible choice, THIN preferred |
| Code complexity | window management, edge handling, ordered output | simplest; RECRUIT is the simplest and can reuse `libisyntax_tile_read` |
| Verified so far | per-frame exactness (via the windowed builder), decode-once and exactness of the full-width band IDWT in `reconstruct_strips.py` (bit-exact); **shared strips for frame rings not implemented** | per-frame exactness verified; RECRUIT only simulated |

**Recommendation [hypothesis]:** implement the frame task as in THIN (self-contained: gets coefficients and rings from a block-row window), but schedule frames in raster bands over one shared window per file, i.e. STRIP; run one such worker per file and many files in parallel. This keeps the CPU at the STRIP minimum, bounds memory by image width, and still has restartable, independently testable frame tasks. Keep THIN with a block cache for partial or random-access jobs. Do not use RECRUIT beyond prototypes.

The bigger lever for a million files is not the halo strategy but Tier-1 (65 to 80 % of CPU), and for the glass-heavy files the area-proportional transforms: omitting near-constant glass frames (13 to 34 % of the pyramid, lossy, see section 8h of the transcoding document), a fast path for frames whose whole footprint is empty (no IDWT needed), and a faster Tier-1 implementation would each shift the total more than the choice between STRIP and THIN.

## 6. Caveats

- The footprint simulation is exact for the footprint, not for time. Times rest on `1.isyntax` timings (warm page cache, SSSE3 build; the AVX2/`-march=native` lifting may be faster) and on the Tier-1 fit from a loaded machine.
- Coded pyramid = 2.3 x the file and the 65 % sharing saving are carried over from one file; the glass-heavy files probably code smaller, which makes transforms a larger share than shown.
- Edge frames are excluded (0 to 9 % of frames per layer); they use the local-transform path of the plan.
- Not implemented: shared-strip ring computation across frames, a C++ strip window, a RECRUIT pipeline.

## 7. Reproduce

```bash
.venv/bin/python3 -B -I tools/v1_analysis/halo_strategy_footprint.py ../example_isyntax_files/N.isyntax BLOCKS.csv 256 5 4   # UNION=0 for the big files
python3 -B -I tools/v1_analysis/halo_strategy_costmodel.py FOOTPRINT.log FILE_BYTES
```
Logs: `$SOURCEROOT/c_decoder_output/strip_analysis/<file>/halo_strategy_footprint_F256.log`.
