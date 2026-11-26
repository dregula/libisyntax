# iSyntax to DICOM WSI pipeline: plan v2

**Version:** v2 (2026-11-26). Supersedes [Pipeline_Plan_v1_Claude.md](Pipeline_Plan_v1_Claude.md), which stays unchanged. Sections marked "unchanged" are carried over from v1 by reference; everything else is restated here.
**Status:** plan plus one built and verified component (the OpenJPEG coefficient-input patch, section 2.1). No pipeline code yet.
Companion documents: [iSyntax_to_DICOM_WSI_Transcoding_Claude.md](iSyntax_to_DICOM_WSI_Transcoding_Claude.md), [Color_Conversion_Claude.md](Color_Conversion_Claude.md).
Tags: **[verified]**, **[spec]**, **[hypothesis]**, **[unverified]**, **open** as in v1.

## Changelog against v1

- **Decided by the owner:** frame size **256 x 256** and strict-range threshold **S = 16**; build the **exact route** that skips OpenJPEG's internal wavelet transform; **always set Lossy Image Compression (0028,2110)**; omitted frames are **displayed as solid white**; a **pure wavelet-domain colour conversion is adequate** for now (no further colour review blocks the pipeline).
- **Built and verified:** OpenJPEG 2.5.4 coefficient-input patch (section 2.1). Codestreams produced from coefficients are byte-identical to stock OpenJPEG's for the same frames.
- **New design item found while planning:** the strict-range omission test needs pixel values, but the pipeline is wavelet-domain only (section 2.3).
- **New consequence flagged:** "solid white" for omitted glass frames brightens the glass by about 22 to 27 levels (section 2.4).
- **New risk:** wavelet-domain colour coefficients must fit the code-block bit-depth the encoder allocates (section 2.1).

## 1. Constraints (updated)

1. Computation restricted to the wavelet domain (no full-image pixel decode).
2. Colour: YCoCg to YBR_RCT **in the wavelet domain, decided adequate**; the colour function stays behind an interface so a later review can replace it.
3. Strict-range omission of glass frames on **layers 0 and 1 only**, `S = 16`, frame size 256.
4. OpenJPEG (patched) now; other encoders possible later behind an interface.
5. C++ pipeline, parallel on a large multi-core CPU; tags are placeholders except those listed in section 5; LABEL and MACRO JPEGs deferred.
6. SOP Instance UID and Series Instance UID from the FIC when present, else placeholders (section 5, unchanged from v1 section 6).

## 2. Changed design points

### 2.1 Exact route: encoder fed with coefficients **[verified]**

OpenJPEG 2.5.4 is patched (new file `patches/openjpeg-2.5.4-coefficient-input.patch`, 4 files, 43 added lines, none removed) with one new call, `opj_set_coefficient_input(codec, TRUE)`, for a raw J2K compressor (`OPJ_CODEC_J2K`), called after `opj_setup_encoder()` and before `opj_start_compress()`. In that mode the data given to `opj_write_tile()` are **int32 coefficients, all components, planar**, and the encoder skips DC level shift, MCT and DWT and goes straight to Tier-1. Everything else (codestream header, SGcod = 1, Tier-1, rate allocation, Tier-2) is stock code.

Input contract for a 256 x 256 frame (one J2K tile, D = 5, 6 resolutions, reversible 5/3):
- Layout per component: OpenJPEG's Mallat layout, row stride 256: LL5 at top-left (8 x 8); for level j (1 = finest) the sub-bands of size `256 >> j` at HL = top-right, LH = bottom-left, HH = bottom-right of that level's `2 * (256 >> j)` square. HL is the band that is high-pass horizontally and low-pass vertically (the `full_pyramid` output order HL, LH, HH).
- Values: already DC-level-shifted and colour-transformed. For 8-bit unsigned samples the shift is **-128 on the Y' LL5 band only** (a constant added to every sample of Y' changes only LL by exactly that constant in the reversible 5/3, and leaves Cb and Cr untouched, as the RCT differences cancel the shift). The detail bands need no shift.
- The codestream header is unchanged: the decoder applies the inverse RCT and adds the shift back.

**Test** (`tools/transcode_tests/coeff_input_test.c`, data from `tools/v1_analysis/make_coeff_test_frames.py`): for three real 256-px layer-0 frames of `1.isyntax` (tissue, scan edge, glass) the coefficients computed with our NumPy transforms equal OpenJPEG's own post-DWT tile data bit for bit (all 3 components), and the coefficient-input codestream is **byte-identical** to the stock codestream of the same pixels (49 978 / 47 357 / 28 846 bytes), and decodes with the stock decoder back to the exact RGB. Reproduce with `tools/transcode_tests/build_patched_openjpeg.sh` (clones v2.5.4 into `ext/`, applies the patch, builds, compiles the test; verified from a clean checkout). One bug found and fixed on the way: the 4-byte sample size must be forced in `opj_tcd_copy_tile_data` as well as in the buffer-size function (a first attempt patched the wrong function and produced garbage that the size check did not catch).

What this buys: no IDWT back to pixels, no forward DWT, no forward MCT per frame. It also makes the wavelet-domain colour route independent of any pixel bounce. What it costs: a patched library to maintain (small), and **it ties the pipeline to OpenJPEG** unless another encoder offers the same entry point. Kakadu and others would need the pixel-bounce fallback of v1 section 2.1, kept behind the encoder interface (the interface has two capabilities: "accepts pixels", "accepts coefficients").

**Risk (new, open):** Tier-1 allocates the number of bit-planes of each band from the component precision (8 bits plus the RCT's extra bit for Cb and Cr, plus guard bits). The test coefficients are those of a true RCT of 8-bit pixels. Wavelet-domain converted coefficients (biased, section 2.2) can in principle exceed the range a true RCT of 8-bit data can produce, and Tier-1 would then mis-code them silently. The pipeline must check `max |coefficient|` per band against the allocated bit-planes and fall back (clamp or pixel bounce) with a logged warning. To be tested on real frames before relying on it.

### 2.2 Colour **[decided]**

Pure wavelet-domain conversion with the truncating halving found in Color_Conversion C1: `Cb = -trunc(Co/2) - Cg`, `Cr = Cb + Co`, `Y' = Y`, applied to LL and every detail band, then the -128 shift on LL5 of Y'. Known behaviour (C2, measured): the decoded Cb is biased by about -0.6 per level (mean -4.9 at D = 5), Cr by about -0.4 per level (-3.0), Y' exact in the interior; RGB after the decoder's inverse RCT differs from the Philips SDK/libisyntax pixels by up to 6-8 levels. The owner accepted this for now. The bias is a property of the output and should be stated in the file's (future) derivation description. Alternatives (per-tile LL correction, pixel-domain conversion, RGB) remain in the colour document and plug in at the colour function.

### 2.3 Omission test needs pixels **[new design item]**

Strict range (max - min of RGB over the frame, worst channel) is a pixel statistic; with a pure wavelet-domain pipeline there are no pixels. Proposal:
1. **Coefficient-domain prefilter (cheap, conservative):** reject a frame as "not glass" if any stored detail coefficient of Y at scale 0 or 1 (layer 0 and 1) exceeds a bound derived from S, or if the LL5 plane (8 x 8) of any colour varies by more than S/2 after conversion. Most tissue frames fail immediately.
2. **Exact check on the survivors only:** run the IDWT for the surviving candidates (they are nearly flat, so this is cheap) and compute the strict range in RGB with the truncating inverse and clamp; omit if `<= S`.
Because the exact check runs on the real pixels the result is identical to the Python analysis (`omission_size_estimate.py`, strict range) for every omitted frame; the prefilter must be proven never to reject a frame that the exact test would omit (property test on the tile CSVs: 3331 of 8217 layer-0/1 frames at S = 16, 12.8 % of the pyramid **[verified, v1 section 2.3]**). Bound derivation: **open**.

### 2.4 Omitted frames shown as solid white **[decided; consequence flagged]**

TILED_SPARSE: omitted frames are simply absent. The viewer's fill is the **Recommended Absent Pixel CIELab Value** (PS3.3 **[spec]**, Transcoding doc 8b) set to white: L* = 100, a* = b* = 0, stored as the standard's 16-bit unsigned encoding (L*: 0xFFFF; a*, b*: 0x8080). Consequences:
- The exact-filler frames (YCoCg = 255, 0, 0, the right and bottom margin bands) are exactly white, so omitting them changes no pixel.
- **Scanned glass is not white:** its median colour is about RGB (228, 230, 233) (Transcoding doc 8h). Omitting a glass frame at S = 16 and showing white makes that area about 22 to 27 levels brighter than the scan, in flat patches the size of whole 256-px frames (3331 frames on layers 0 and 1 of `1.isyntax`). The owner has decided this; viewers that show absent frames as white will therefore draw a slightly brighter background where tiles were dropped, with visible frame-aligned borders against retained glass frames. Retained glass frames next to omitted ones keep their scanned colour. If this proves objectionable, the options are: lower S, or encode a near-white fill into the file by keeping the frames but flattening them (loses the saving). **open as a visual-review item.**
- Layers 2 to 5 are not omitted (decided in v1).

### 2.5 Lossy Image Compression tag, always set **[decided]**

Always write (0028,2110) Lossy Image Compression = `01`, and (0028,2114) Lossy Image Compression Method and (0028,2112) Ratio. Values for v2: carry over the source's entries from the FIC when present (`DICOM_LOSSY_IMAGE_COMPRESSION_METHOD` `PHILIPS_DP_1_0`, ratio 15 in the sample FICs) and append `ISO_15444_1` with a computed ratio (uncompressed size over codestream size) for this file's own codestreams; the placeholders from the FIC are by group/element as in section 5. The transfer syntax stays `1.2.840.10008.1.2.4.90` (lossless-only syntax carrying data that are lossy relative to the scan: the tag is how that is declared). The exact method value list and the ratio rule are **open details**; the tag itself is always present.

## 3. Unchanged from v1

- Section 2.4 (encoder interface, extended by the coefficient-input capability above), section 3 (DICOM tool review; DCMTK for metadata, own streaming writer for Pixel Data, Extended Offset Table), section 4 (architecture: decode each codeblock once, sliding window, strip assembler, frame task, ordered writer, patched fixed-size offset table), section 5 (parallelism), section 6 (UID rules), section 8 (tests, plus the new ones below).
- The frame task in section 4 becomes: gather coefficients, convert colour (section 2.2), shift LL5 of Y', arrange in Mallat layout, omission test (section 2.3), `opj_write_tile` in coefficient mode. **No IDWT bounce** for frames that are kept.

## 4. Parameters fixed for the first implementation

| Parameter | Value |
|---|---|
| Frame size | 256 x 256 |
| Decomposition levels | D = 5 (6 resolutions) |
| Code-blocks | 64 x 64 (as measured in 8g) |
| Transfer syntax | 1.2.840.10008.1.2.4.90, photometric YBR_RCT, SGcod = 1 |
| Omission | strict range S = 16, layers 0 and 1 only |
| Absent-frame fill | white (CIELab L* 100, a* 0, b* 0) |
| Lossy Image Compression | always `01` |
| Colour | wavelet-domain, truncating halving |
| Layer files | six (one SOP Instance each), one series |

All are command-line parameters with these defaults except the transfer-syntax triple.

## 5. Placeholders and UIDs (unchanged in substance)

Required tags other than those above and below carry placeholder values. SOP Instance UID (0008,0018) and Series Instance UID (0020,000E) come from the FIC by group/element when present; otherwise `2.25.<128-bit decimal>` derived from the SHA-256 of the input plus a purpose string and the layer number; one SOP Instance UID per layer, series and frame-of-reference UIDs shared; Media Storage SOP Instance UID equals the dataset's.

## 6. Additional tests (v2)

- `coeff_input_test` must pass for every frame class (tissue, scan edge, glass) **and** on frames produced by the pipeline's own wavelet-domain colour route: there the reference is the pixel image obtained by the stock decoder from the produced codestream, compared with the Python wavelet-domain route (`convert_pyramid`, `q = trunc`), expecting identical pixels.
- Bit-depth guard test: `max |coefficient|` per band versus the allocated bit-planes on 1000 random frames of each file.
- Omission property test: prefilter never rejects a frame the exact strict-range test omits (both layers, S = 16, `tile_stats_F256_L0.csv` and `_L1.csv`).
- White-fill test: a reader that fills absent frames with the Recommended Absent Pixel value shows exactly (255, 255, 255) there.

## 7. Open items after v2

1. Prefilter bound for the omission test (section 2.3).
2. Bit-depth guard behaviour (clamp or fallback) after measuring real frames (section 2.1).
3. Visual review of white-filled glass areas (section 2.4); S may need lowering.
4. Lossy-compression method list and ratio rule (section 2.5).
5. Everything listed as unverified in v1 (DCMTK JPEG 2000 / WSI / Extended Offset Table support).
6. Kakadu or another encoder with a coefficient entry point; until then those encoders use the pixel bounce.
7. v2 data-model decoding in C++ (the Python tools read v1 only); PAS (data model 100.4) is the test file.
8. Axis mapping and Image Orientation (Slide).
