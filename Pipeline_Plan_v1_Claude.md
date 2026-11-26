# iSyntax to DICOM WSI pipeline: plan v1

**Version:** v1 (2026-11-26). Each iteration of the plan is a new file `Pipeline_Plan_v<N>_Claude.md`; older versions are kept unchanged. A new version starts with a changelog against the previous one.
**Status:** plan only, no pipeline code exists yet. Companion documents: [iSyntax_to_DICOM_WSI_Transcoding_Claude.md](iSyntax_to_DICOM_WSI_Transcoding_Claude.md) (geometry, halo algorithm, tiling, spec reading, EBCOT sizes, tile omission), [Color_Conversion_Claude.md](Color_Conversion_Claude.md) (colour conversion, C1 to C3).
Tags: **[verified]** = measured in this repo; **[spec]** = from the standards; **[hypothesis]**; **[unverified]** = from web search or memory, not checked against the source; **open** = a decision or measurement still missing.

## Changelog

- v1: first plan, from the owner's four constraints (below), the measurements in the two companion documents, and a review of open-source DICOM tools.

## 1. Constraints given by the project owner (v1)

1. Computation is restricted to the **wavelet domain** (no full-image pixel decode).
2. Colour transform YCoCg to YBR (YBR_RCT) in the wavelet domain **for now**; an exhaustive manual review of colour strategies comes later, so the colour step must be swappable.
3. **Strict-range omission** of tiles that contain no tissue (glass), **only on layers 0 and 1**.
4. **OpenJPEG** in the pipeline, with the possibility of other JPEG 2000 encoders (Kakadu, DCMJP2K) later.
5. Goal: a complete iSyntax-to-DICOM pipeline, in **C++**, parallel on a large multi-core CPU.
6. For now: DICOM tag conversion is out of scope, required tags are filled with placeholder values; LABEL and MACRO JPEG files are deferred.
7. SOP Instance UID and Series Instance UID are taken from the FIC when present, otherwise placeholder UIDs are stored. The FIC's `Attribute` elements for standard DICOM tags carry `Group` and `Element` attributes that map directly to the DICOM tag.

## 2. Review of the constraints against the measurements

**2.1 Wavelet domain with a stock encoder.** OpenJPEG and Kakadu accept pixels, not coefficients. So the wavelet-domain route is: convert the frame's coefficients (LL and detail bands) to the target colour components, IDWT them to the frame's pixels, and let the encoder's own forward 5/3 reproduce those coefficients before EBCOT. OpenJPEG's forward transform matches our NumPy conventions, so the bounce is exact **[verified, section 8g of the Transcoding doc]**. The cost is a redundant IDWT and DWT per frame, small next to Tier-1. Exact alternatives, only if the overhead ever matters: a patched OpenJPEG that skips its DWT stage, or an own Tier-1 coder. **open** (decision 2 in section 9).
- The encoder input is 8-bit: pixels outside 0..255 must be clamped (up to 2951 of 1 M pixels in one 1024-px edge frame **[verified]**; Y reaches 257 in the pixel domain and 266 in the top LL). The Philips SDK itself does not clamp below 0 (Color_Conversion C1), so a few pixels can never match it.

**2.2 Colour.** The Philips SDK's inverse is truncating (`/2`) and matches libisyntax's scalar path in 100 % of 22 tiles **[verified, Color_Conversion C1]**; the scanner's forward is therefore taken to be truncation (inferred, not observable). Use `q = trunc` in `Cb = -q(Co) - Cg`, `Cr = Cb + Co`, `Y' = Y`. The wavelet-domain conversion is biased (Cb mean -4.9, Cr -3.0 at D = 5; near-lossless at best) **[verified, C2]**. The colour step is one function behind a stable interface so that C2/C3 alternatives (pixel-domain exact, hybrid, RGB, per-tile LL correction) can be swapped in.

**2.3 Strict-range omission, layers 0 and 1 only.** A tile is omitted when its strict RGB range (max - min over the frame, worst channel) is <= S, i.e. no speck, no tissue. Savings as a share of the whole six-layer pyramid, computed from `transcoding_review/tile_stats_F*_L*.csv` **[verified]**:

| frame size | S = 8 | S = 16 | S = 32 |
|---|---|---|---|
| 256 px | 9.4 % (19.2 MB, 2590 tiles) | 12.8 % (26.1 MB, 3331 tiles) | 19.3 % (39.4 MB, 4523 tiles) |
| 512 px | 1.0 % | 3.0 % | 7.2 % |
| 1024 px | 0.0 % | 0.1 % | 1.0 % |

(pyramid totals: 204.4 / 198.6 / 196.4 MB; exact-filler tiles on layers 0 and 1 are 506 / 49 / 0 and save 0.08 / 0.01 / 0 MB.) The omission pays only at 256 px; at 512 px it nearly vanishes. **Frame size and omission therefore have to be decided together.** Proposal: parameters `frame_size` and `S`, defaults 256 and 16.
- Omitted frames must be filled by the viewer: one Recommended Absent Pixel colour about RGB (228, 230, 233) fits scanned glass, not the pure-white filler bands (8h). **open**
- Specks under the robust rule are up to 4 um across (Transcoding doc 8h), none with the strict rule; the strict rule deletes no speck by construction.

**2.4 Encoder interface.** One small interface: takes an 8-bit, 3-component frame (plus options: reversible 5/3, D levels, code-block size, progression order, MCT on) and returns a raw codestream (no JP2 header, SGcod = 1) **[spec: PS3.5, section 8c]**. Implementations: OpenJPEG now; Kakadu, Grok (HTJ2K needs a different transfer syntax), OpenJPH, DCMJP2K behind the same interface. Run one thread per encode and parallelise across frames; do not use the encoder's internal threading.

## 3. Review of open-source DICOM tools

| Tool | Use in this project |
|---|---|
| **DCMTK** (C++) | Metadata dataset: data dictionary (VR from group/element), UID generation, file-meta handling, writing the non-pixel part of the file. |
| **GDCM** (C++) | Alternative toolkit, ships an OpenJPEG-based JPEG 2000 codec; no reason to prefer it for the metadata. |
| **highdicom** (Python) | Reference for the full required-tag list of the VL Whole Slide Microscopy Image IOD; quick prototyping of the tag set. |
| **wsidicomizer** (Python) | Existing converter, reference for tiling and the tag set. |
| **dicom3tools `dciodvfy`** | Validation of the output IOD. |
| **Orthanc, QuPath or similar** | Read test of the files. |

Caveats on DCMTK:
- **JPEG 2000 codec: [unverified].** As I understand it the open-source DCMTK has no JPEG 2000 codec; OFFIS sells the DCMJP2K extension (reported to be Jasper-based). The 3.7.0 documentation has a `dcmcjp2k` page, but I could not read it (site bot check). Check dcmtk.org before relying on either statement. Third-party OpenJPEG-based DCMTK codecs exist (fmjpeg2koj).
- **Whole-slide IOD: [unverified].** I found no DCMTK class for the VL Whole Slide Microscopy IOD; plan to build it from `DcmDataset`/`DcmItem`.
- **Memory:** DCMTK holds an encapsulated pixel sequence in memory. Layer 0 of PAS (1.6 GB iSyntax) would be a similar size in codestreams. Plan: DCMTK for everything except Pixel Data, own streaming writer for Pixel Data.
- **Extended Offset Table (7FE0,0001) and Lengths (7FE0,0002)** are needed when the encapsulated Pixel Data exceeds 4 GB, since the Basic Offset Table cannot hold it; they require the Basic Offset Table item to be empty and each frame in one fragment **[spec: CP 2388, PS3.5 A.4]**. Whether DCMTK writes them: **[unverified]**, check `dcmdata/libsrc/dcpixel.cc`, `dcpixseq.cc`.

## 4. Architecture of the C++ pipeline

**Governing constraint: decode every codeblock once.** The windowed decoder of the Python prototype reads about 1 MB per 256-px layer-0 frame (22 codeblocks per colour from 5 chunks, Transcoding doc section 5). Done per frame, PAS (400 x 304 layer-0 tiles, about 121 600 frames) would read about 120 GB for a 1.6 GB file (extrapolation, **[hypothesis]**). So the pipeline decodes each chunk once and keeps decoded coefficients in a sliding window of rows; frames read from the window.

Stages:

1. **Open and index.** Parse the FIC (XML), block header table and chunk layout; memory-map the file. Implemented as a **new coefficient-level module** on top of the libisyntax headers, so the existing decoder is untouched (CLAUDE.md: regression tests compare decoded pixels; new code goes in new files). Both data models (v1 and v2) must be handled; the Python tools only cover v1.
2. **Decode.** A reader hands chunks, in the order the current strip needs them, to a pool that runs the Hulsken decode per codeblock into per-scale row-window buffers (i16, three colours).
3. **Strip assembler.** For each frame row: copy interior detail coefficients; recompute the halo ring (depth 2, use 4 for margin) from a windowed top-down IDWT of rectangles around the frame edges (Transcoding doc sections 3 and 4 **[verified]**). Never reconstruct the whole image.
4. **Frame task** (independent per frame): gather coefficients; convert colour (swappable); IDWT to pixels; clamp; strict-range omission test (layers 0 and 1 only); encode through the encoder interface.
5. **Ordered writer.** A reorder buffer per layer feeds one writer per layer; the six layers (one SOP Instance each, all in one series) are written to six files concurrently.
6. **File layout.** DCMTK builds the metadata dataset with a fixed-size placeholder for the Extended Offset Table (8 bytes per frame; the frame count is known beforehand, including omitted frames). The writer then streams the encapsulated Pixel Data (empty Basic Offset Table, one fragment per frame, sequence delimiter) and patches the offset table at the end. Sparse tiling (TILED_SPARSE) lists the position of every present frame in the per-frame functional groups (Plane Position (Slide)); omitted frames are simply absent.

Layer geometry and alignment (frame origins, `(i0 + 3) % 2^D == 0`, sample positions `ox + 2^k (i + 3) - 2`, pixel spacing 0.25 um x 2^k) are as in Transcoding doc 8a **[verified]**; the iSyntax-to-slide axis mapping and Image Orientation (Slide) are still **open**.

## 5. Parallelism on a large multi-core CPU

- **Work units:** chunk decode (I/O plus Hulsken), strip assembly, frame task. Frames are independent codestreams, so the frame task is embarrassingly parallel once the strip is available; Tier-1 encoding dominates the CPU time.
- **Scheduler:** oneTBB or a custom work-stealing pool with a bounded number of in-flight strips (a semaphore on strip rows) to keep memory fixed. Do not enable the encoder's internal threads.
- **Memory:** per-thread buffers and arenas reused across frames (libisyntax already ships `ltalloc` and arena helpers; mimalloc is an option); avoid allocator contention.
- **NUMA:** on multi-socket machines, one decoder pool per socket owning a row range, so strips stay local.
- **I/O:** memory-map the input, advise sequential access; write each layer file from a single writer thread.
- **Expected scale (order of magnitude, extrapolation, not a measurement):** the Python prototype encoded every tile of every layer of `1.isyntax` in about 4 minutes on 16 processes (`omission_size_estimate.py`); PAS is about 18 times larger, so roughly an hour on 16 cores in C++, less with more cores. Measure early with a prototype before trusting this.

## 6. UIDs and the FIC

- Match FIC attributes by `(Group, Element)`, not by name: (0008,0018) SOP Instance UID, (0020,000E) Series Instance UID. The VR comes from the DICOM dictionary (DCMTK), not from the FIC's `PMSVR` field (IString, IDouble, IStringArray).
- **None of the six example FICs in `example_FIC_files/` nor the new `PAS.fic` carries either UID as an attribute** (the only matches of "UID" are inside base64 payloads) **[verified]**. The owner reports that some files have them: the lookup must work when present and fall back when absent.
- Fallback: `2.25.<128-bit decimal>` derived from the SHA-256 of the input file plus a purpose string and the layer number, so that re-runs give identical UIDs.
- One SOP Instance UID **per layer** (each layer is its own instance). A single UID from the FIC can belong to at most one layer; the others get derived UIDs. The Series Instance UID and the Frame of Reference UID are shared by all layers. Media Storage SOP Instance UID in the file meta must equal the dataset's SOP Instance UID.
- Other FIC attributes with a group/element (e.g. `DICOM_MANUFACTURER` 0008,0070; `DICOM_ACQUISITION_DATETIME` 0008,002a; `DICOM_ICCPROFILE` 0028,2000; `DICOM_LOSSY_IMAGE_COMPRESSION*` 0028,2110/2112/2114) use the same mechanism and can be carried over later. In v1 of this plan they are not converted.

## 7. Scope for the first milestone

In: FIC and codeblock reading (v1 and v2), coefficient window, halo strips, wavelet-domain colour (swappable), IDWT bounce, clamp, strict-range omission on layers 0 and 1, OpenJPEG encoder behind the interface, streamed DICOM writer with Extended Offset Table, six layer files with placeholders for required tags, UID logic of section 6.
Out (deferred by the owner): tag conversion, LABEL and MACRO JPEG files. Out (not yet decided): exact colour routes, slide axis mapping, display of omitted frames, Kakadu/DCMJP2K implementations.

## 8. Test strategy

- **Bit-exactness of the geometry:** decode each written frame with OpenJPEG and compare with the Python windowed decoder (`dicom_tile_tests.py`, `windowed_idwt_halo.py`) for Y, Co and Cg: must be identical for the interior and the halo ring.
- **Colour:** compare against the libisyntax pixel output of the truncating build (`c_decoder_output`). Not equality: the wavelet-domain route is biased, so define metrics (share of pixels within 0, 1, 2 levels; max error; per-channel mean error) and compare them with the Part A numbers of C2/C3.
- **Structure:** `dciodvfy` on every file; read the files with an independent reader (highdicom/pydicom or a viewer); check that the Extended Offset Table points to each frame's first Item Tag.
- **Regression:** the same input must produce byte-identical files (deterministic UIDs, ordered writer).

## 9. Open decisions and risks

1. **Frame size and `S`** (256 px and S = 16 proposed); depends on which viewers must read the files.
2. **Pixel bounce versus exact route** (patched OpenJPEG or own Tier-1).
3. **Lossy Image Compression (0028,2110)** should be set truthfully even while other required tags are placeholders: transfer syntax `.90` is lossless-only, but the iSyntax source is lossy (about 15:1) and the wavelet-domain colour route and omission add further loss.
4. **Fill colour for omitted frames** (Recommended Absent Pixel CIELab value) and how viewers show gaps.
5. **Colour fidelity target** (C2/C3 decision) and the exhaustive review planned by the owner.
6. **Unverified facts to check first:** DCMTK JPEG 2000 support, DCMTK Extended Offset Table support, DCMTK WSI support.
7. **Risks:** read amplification if the sliding window is not implemented (about 75 times the file size per layer 0); v2 decoding is documented only by this repo's C code; one real v1 file and one v2 file (PAS) are the only test inputs; level 0 and 1 SDK exports for validation do not exist (the level-1 TIFFs are truncated).
