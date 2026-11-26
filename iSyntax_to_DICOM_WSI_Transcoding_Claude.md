# iSyntax v1 to DICOM WSI: wavelet-domain transcoding exploration

Date: 2026-11-25. Status: exploration; the halo algorithm is verified on `1.isyntax`, everything marked "open" is not.
Tags: **[verified]** = measured in this repo on the real file or by exhaustive/numerical test; **[derived]** = follows from the
geometry in `tools/v1_analysis/reconstruct_np.py` and was then checked; **[hypothesis]**; **[spec]** = from the Philips PDFs.
Companion notes: [Color_Conversion_Claude.md](Color_Conversion_Claude.md) (colour conversion, split out 2026-11-26), [Third_Party_Documentation.md](Third_Party_Documentation.md), [Example_iSyntax_File_Claude.md](Example_iSyntax_File_Claude.md).

## 1. Target (as given by the project owner)

- DICOM WSI, sparse tiles (TILED_SPARSE) in both source and target; photometric interpretation YBR_RCT; transfer syntax JPEG 2000
  lossless (reversible 5/3, 1.2.840.10008.1.2.4.90). Frames are independent codestreams. (Checked against PS3.3 and PS3.5 in sections 8b and 8c.)
- **One resolution layer = one SOP Instance**, all layers under one Series Instance UID. A layer must be self-contained: every frame carries
  its own LL band, with no reference to any lower-resolution layer.

## 2. Core finding: the two "overlaps" are decoder conventions, and only a thin ring differs

- iSyntax stores one seamless 5/3 transform of the whole padded image. Its "neighbour overlap" (4 coefficients, `ISYNTAX_IDWT_PAD_L/R`) exists only in the decoder.
- A JPEG 2000 frame transforms only its own pixels and mirrors them at its edges. Neither format stores redundant data.
- Interior coefficients are identical in both. Only about 2 coefficients per edge, per level, per band differ. **[verified]**
  - 1-D synthetic test (`halo_1d.py`): at most 2 coefficients per edge and level, at every level, for 256/512-px frames with 3 or 5 levels.
  - Real file, Y channel (`frame_interior_test.py`), 256-px frame, 5 levels: 0.7 / 3.2 / 9.9 / 22 / 44 % of detail coefficients differ from a plain copy at levels 1..5,
    all within 2 coefficients of an edge. Overall 2.0 % (256 px) and 0.9 % (512 px) of all detail coefficients.

## 3. Geometry **[derived, then verified]**

- iSyntax IDWT output `out_s` (scale s, 0 = full resolution) has sample i at padded level-0 position 2^s (i+3) - 2. Padded = valid + 189 for this file.
- That is the standard JPEG 2000 grid with the canvas shifted by +2: high-pass at 2^s * odd, low-pass at multiples of 2^(s+1).
- **Frame alignment:** for D decomposition levels, the frame's first sample `out_k[i0]` must satisfy `(i0 + 3) % 2^D == 0`. At layer 0 in valid-image
  coordinates that is x, y = 1 (mod 2^D) for D up to 6. Shifting the frame by one pixel breaks it: 38 to 99 % of coefficients differ. This means the DICOM
  pixel matrix must be offset by one pixel relative to the valid image, or each frame codestream must carry an image offset (the latter may hurt interoperability: **hypothesis**).
- **Index maps** (per axis, local level j = 1..D uses iSyntax scale s = k + j - 1, `q = i0 + 1`): high-type index `n_h = m - 1 + ((q+2) >> j)`,
  low-type `n_l = n_h - 1`. HL uses (rows low-type, cols high-type), LH the opposite, HH both high-type. Getting this wrong looks like a 13 % mismatch (my first test did).
- **Layer k of the DICOM pyramid = the pixel plane `out_k`.** Its J2K levels 1..D use iSyntax scales k..k+D-1. Scale 5 is the top stored scale here (LL stored
  only there), so for k + D - 1 > 5 the extra levels are an ordinary local forward transform of the frame's own LL. Nothing is read from another layer in the output. **[verified]** for layers 0 to 5.
- Per-layer pixel spacing and origin: samples sit at 2^k (i+3) - 2 in level-0 padded units, so layers have sub-pixel offsets relative to each other. Not worked out; needed for DICOM geometry metadata. **open**

## 4. Algorithm (strip halo recomputation) **[verified]**

1. Copy the interior detail coefficients straight from the stored H blocks (index maps above).
2. For each level, recompute only a ring of depth h along each frame edge from strips of width W = 2h + 4 with mirror extension (vertical then horizontal, T.800 order). h = 2 is sufficient, h = 4 is a margin; h = 1 fails.
3. Input strips: the pixel ring of `out_k`, rings of the global LL arrays at each level, and the small full LL_D. They come from a **windowed top-down IDWT**
   (`windowed_idwt_halo.py`, class `Windowed`): `out_s[i]` needs L[n], H[n] for n in [(i0>>1)-1, ((i1-1)>>1)+3), and L = `cur_s[n] = out_{s+1}[n-1]` is requested recursively up to the stored top LL.
4. Everything the builder does not read is filled with random garbage, so exactness proves only rings were read.
5. If the ring covers a whole level (n/2 <= 2h) a plain local forward transform is used there.

Results (`windowed_idwt_halo.py`, Y, Co and Cg, one tissue frame of 256 px per layer 0 to 5, plus a 512-px tissue frame and tissue-edge, blank-glass and scan-edge frames at layer 0, D = 5):
windowed IDWT equals the full decode on random rectangles at every scale (0 mismatches); built coefficients equal a full per-frame forward transform (0 mismatches); a standard J2K inverse returns the iSyntax-decoded pixels bit for bit.
Strip work is 14 to 39 % of a full forward transform (9 % for 512 px at h = 2). The saving is modest for small frames.

## 5. Read footprint of the windowed approach (h = 4, one frame) **[verified]**, versus the earlier conservative closure

| frame | codeblocks per colour | chunks | compressed, 3 colours |
|---|---|---|---|
| layer 0 tissue 256 | 22 | 5 | 1.02 MB |
| layer 0 tissue 512 | 34 | 5 | 1.20 MB |
| layer 0 scanned glass | 18 | 2 | 0.23 MB |
| layer 3 tissue 256 | 12 | 2 | 1.17 MB |
| layer 5 tissue 256 | 8 | 4 | 0.49 MB |

`count_chunks.py` (all neighbours whole, full LL chain) gave 8 to 37 chunks and a median of about 300 codeblocks (3.2 MB) for one scale-0 tile. That was an upper bound; for a 256-px frame at layer 0 the windowed decoder reads 66 codeblocks (22 per colour, about a fifth) and about a third of the bytes. One frame per layer only; a systematic count is open.

## 6. Colour space (YCoCg-R to YBR_RCT) **[moved]**

Moved to [Color_Conversion_Claude.md](Color_Conversion_Claude.md), section C1. Summary: the YCoCg-R to YBR_RCT algebra (Y' = Y, Cb = -Co/2 - Cg, Cr = Co/2 - Cg); the halving convention is truncation (inferred from the Philips SDK output, C1); a wavelet-domain conversion is biased, not "off by 1". RGB as target: section C3 there.

## 7. What the Philips documents say about chunks

- "Chunk" is libisyntax's word; Philips' v2 header says **cluster**, and the format PDF uses neither. Hulsken 2016 describes clusters as "regions of the wavelet transformed image data spanning both space as well as wavelet scale" for "fast streaming storage of images from the scanner", with tiles inside for "fast random access by the viewer". **[spec]**
- The format PDF (Figure 9 text) says codeblocks are written left to right, top to bottom, top level to base level, then colour. The real file groups 63-block chunks instead. Unresolved. **open**
- A chunk is a storage unit, not a decode unit: reconstruction crosses chunk boundaries (neighbour margins and the LL chain), so one frame touches several chunks. The AI note claiming "no cross-chunk dependencies" is true only for entropy decoding.

## 8. Not done

- A decision on colour fidelity (Color_Conversion_Claude.md section C2), now with real sizes in 8g.
- Sparse-frame rules beyond the blank-candidate test (section 8a), and what a frame whose neighbour is omitted should assume (the decoder uses black H blocks and a dummy LL for missing neighbours). Edge tiles and the matrix-corner tiles are covered in 8a.
- A real JPEG 2000 codestream (Tier-1/Tier-2 coding, guard bits, OpenJPEG decode; also T.800 itself, to confirm the grid, extension and order conventions) and the 8-bit range question (Y reaches 257 in the pixel domain, 266 in the top LL).
- Systematic footprint and timing over many frames; the minimal depth h = 2 is observed, not proven.
- DICOM metadata: the centre-versus-corner choice, the iSyntax-to-slide axis mapping and Image Orientation (Slide) (section 8b), UID structure, Dimension Index Sequence content.

## 8a. Tile grids, absolute slide coordinates and tile tests (256, 512, 1024) **[verified unless marked]**

`tools/v1_analysis/dicom_tile_tests.py` (decodes all three channels for ground truth; the tiles themselves are built by the windowed decoder).

- **Anchor to the file.** The header gives dimension unit MicroMeter, scale factor 0.25 um per unit, origin x = 77715, y = 20614 units (= 19428.75 um, 5153.5 um; the v2 FIC states 19428.8 / 5153.5 um), so block-header coordinates are already absolute level-0 positions.
  Sample i of layer k is at `ox + 2^k (i+3) - 2` units (centre of that sample); layer pixel spacing is exactly 0.25 um x 2^k. For every tested tile the position computed from the tile index equals the one implied by the file's own block header coordinates (`header-anchor ok`; this checks my index-to-position mapping against the file, it does not check physical placement).
- **Layer-to-layer offsets.** Co-located samples satisfy i = 2 i' + 3 (layer k to k+1). Numerically: the LL band of a layer-k tile's forward transform equals the layer k+1 pixels at the index predicted from absolute position (`cross-layer LL ok`, ring excluded). So the layer offsets are exact, not approximate.
- **Per-layer matrix and tile origin.** Tile origin `i0 = i_base + n F` with `(i_base + 3) % 2^D == 0`, `i_base` the largest such value not after the first sample overlapping the image. Resulting total pixel matrices (D = 5), layer 0: 25253 x 21157 (31 black padding pixels precede the image at top and left; the first valid pixel is column 32), layer 1: 12643 x 10595, layer 2: 6322 x 5298, layer 3: 3178 x 2666, layer 4: 1589 x 1333, layer 5: 795 x 667. For layers 3 to 5 `i_base = -3`, i.e. the first three samples do not exist in the decoded array.
- **Per-tile DICOM values computed:** 1-based column/row in the total pixel matrix (`1 + tx F`, `1 + ty F`) and the slide offset in mm of the top-left sample centre and of its corner (centre minus 2^k/2 level-0 px). Which of the two DICOM wants, the orientation of the slide axes relative to the iSyntax x/y, and whether the iSyntax origin is the DICOM slide origin remain open after reading PS3.3 (section 8b). **open**
- **Tile counts** (cols x rows tiles): F = 256: 99x83 (layer 0), 50x42, 25x21, 13x11, 7x6, 4x3; F = 512: 50x42, 25x21, 13x11, 7x6, 4x3, 2x2; F = 1024: 25x21, 13x11, 7x6, 4x3, 2x2, 1x1.
- **Tests per tile size and layer** (sampled tiles: origin, last, top edge, left edge, tissue, tissue edge, blank glass, scan edge): interior tiles built with the strip algorithm (h = 2 and 4) give 0 coefficient mismatches against a full per-tile forward transform in all three channels and decode standalone from their own coefficients to exactly the iSyntax pixels; a 2x2 block of independently decoded tiles reproduces the iSyntax decode exactly (0 mismatching samples). All of this passed for 256, 512 and 1024, layers 0 to 5, with the origin tile of layers 0 to 2 at the matrix corner included.
- **Edge tiles** (partly outside the decoded array: the last row and column of tiles, and every tile of layers 3 to 5 touching the array start) have their outside pixels filled with 0 (black padding), then are transformed locally. Their exactness is by construction: what the test shows is that the in-array pixels equal the iSyntax decode. Counts: F = 256 layer 0: 0 edge of 8217; layer 1: 91 of 2100; layer 2: 45 of 525; layer 3: 44 of 143; layer 4: 22 of 42; layer 5: 10 of 12. With F = 1024 every tile of layers 3 to 5 is an edge tile.
- **Sparse tiles.** This file has no empty codeblocks, so sparseness is a transcoder policy. A tile is a blank candidate if all stored detail coefficients are zero in all channels and the last LL band is constant. F = 256, layer 0: 457 of 8217 (5.6 %), F = 512: 49 of 2100 (2.3 %), F = 1024 and every layer 1 to 5: none. All candidates have LL = (Y 255, Co 0, Cg 0), i.e. white, and **all 457 lie in two margin bands: the last two tile columns (about 421 px at the right edge) and the last three tile rows (about 677 px, 3.2 % of the height); non-blank tiles span rows 0-79 and columns 0-96 of the 83 x 99 grid.** Real scanned glass is Y about 230 (mean 230.5, std 0.6 in the 'glass' test region), so exact (255, 0, 0) is synthetic filler, not scanned background. The 'blank glass' regions used in the size and footprint tests are scanned glass inside the valid-data envelope (x 189..25411, one rectangle; y range not checked) and contain no filler tile. 12 of 12 sampled candidates per case were pixel-constant (the rest are unverified). **Omitting a blank frame is only lossless if the viewer fills with exactly that white** (YBR_RCT Y 255, Cb 0, Cr 0); PS3.3 does not define the display of gaps in TILED_SPARSE but offers Recommended Absent Pixel CIELab Value (section 8b). **open**
- Only sampled tiles were tested, not every tile. Running every tile at F = 1024 (525 tiles at layer 0) would be feasible.

## 8b. What PS3.3 says (read 2026-11-25) **[spec]**

Sources: the original Supplement 145 text (read raw from its PDF, 2010 draft, so it may differ from the current edition) and current-edition PS3.3 pages A.32.8, A.1.2.30, C.7.6.3.1.2, C.7.6.17, C.7.6.17.3, C.8.12.2, C.8.12.4, C.8.12.6, C.8.12.14, **read through a summarising fetch tool, not verbatim**. Not read: PS3.5 (JPEG 2000 transfer syntax constraints), PS3.4, the exact figures C.8-16/17.

Settled:
- **Regular tiling, frames may overhang.** Column Position In Total Pixel Matrix (0048,021E) "shall be an integer multiple of Columns plus a constant" (same for rows); frames may extend beyond the matrix, and pixels beyond the edge may be a fixed padding value. So edge tiles with filler are allowed, and a per-layer constant offset (our `i_base`) is allowed.
- **Sparse:** some tiles may be absent; frame order is unconstrained; each frame carries its position in Plane Position (Slide) (required when not TILED_FULL; X, Y in mm, Z in um, Z required even for a single focal plane); TILED_SPARSE also requires the Dimension Index Sequence. Recipients may not assume anything about positions or overlap. For TILED_SPARSE the standard does not say how a gap is displayed, but the optional **Recommended Absent Pixel CIELab Value (0048,0015)** carries the colour to show for unencoded pixels (a recommendation only). That is the hook for white blank tiles.
- **Pyramid:** all instances of one Multi-Resolution Pyramid share the Frame of Reference and the Series, one pyramid per Series; each layer has uniform Pixel Spacing (a layer may even be split over several instances); VOLUME image type may be used for each layer.
- **Pixel Spacing is centre to centre** (mm), which matches our sample spacing 0.25 um x 2^k = 0.00025 x 2^k mm. Imaged Volume Width/Height (mm) are required for VOLUME images (for this file 25222 x 0.25 um = 6.3055 mm by 5.2815 mm).
- **Total Pixel Matrix Origin Sequence** gives the location of the top-left pixel ("pixel 1\1") of the matrix; **Image Orientation (Slide)** (0048,0102, direction cosines of the first row and column, required with Plane Position (Slide)) relates the matrix to the slide axes. Total Pixel Matrix Columns/Rows describe the extent as if the volume were encoded non-sparse.
- **Slide Coordinate System:** origin nominally at the slide's lower-left corner; X increases toward the right edge, Y toward the label edge (or top edge), Z in um from the glass surface; X, Y offsets in mm. Supplement 145 notes that this slide (X, Y) system is **rotated 180 degrees from the conventional image (row, column) orientation with the label on the left**, so Image Orientation (Slide) is generally not the identity.
- **YBR_RCT (C.7.6.3.1.2):** Y = floor((R + 2G + B)/4), Cb = B - G, Cr = R - G; inverse G = Y - floor((Cr + Cb)/4), R = Cr + G, B = Cb + G. Encapsulated data only, Samples per Pixel 3, Planar Configuration 0, valid at any Bits Allocated; for JPEG 2000 the DC level shift is applied **before** the forward colour transform and the transformed components may be signed. This matches section C1 of Color_Conversion_Claude.md (the level shift only changes Y).

**Not settled by the text I read (still open):**
1. Whether the X/Y offsets (per frame and for the matrix origin) refer to the **centre or the corner** of the top-left pixel. Slide Coordinates Note 2 says there is no a priori correspondence between pixels and slide coordinates, and the Image Center Point text speaks of the centre of the image pixel plane. The harness reports both; choosing is a project decision.
2. **Whether the iSyntax x/y axes and origin are the DICOM slide axes and origin.** The iSyntax origin (19.43 mm, 5.15 mm) is plausible as a slide position, but nothing in the file says which slide corner it refers to or how its axes map; this is needed to fill Image Orientation (Slide). Unknown without scanner documentation.
3. The statement that a JPEG 2000 frame must be self-contained (PS3.5) was not read.

## 8c. What PS3.5 says about JPEG 2000 (read 2026-11-25) **[spec, via a summarising fetch tool, not verbatim]**

Pages: 8.2.4 (JPEG 2000 transfer syntaxes), A.4 (encapsulation), A.4.4 (JPEG 2000 encapsulation). Not read: ITU-T T.800 itself, the Extended Offset Table definition, PS3.5 8.2 note details beyond these pages.

- **Transfer syntax:** 1.2.840.10008.1.2.4.90 is "JPEG 2000 (Lossless Only)": reversible transform, no quantization. 1.2.840.10008.1.2.4.91 permits lossy. (This confirms the UID I had from memory.)
- **Raw codestream only:** "The optional JP2 file format header shall NOT be included"; the DICOM attributes play the header's role. Where the stream's own characteristics disagree with the DICOM attributes, the stream controls decoding, and the attributes must be consistent with it.
- **Colour:** if the reversible multi-component transform is used (signalled by SGcod = 1 in the COD marker), Photometric Interpretation **shall be YBR_RCT**; with no component transform the components must match the photometric interpretation, and RGB is allowed in that case (the standard notes it as unusual and poorly compressing, citing whole slide microscopy). Colour images: Samples per Pixel 3, Planar Configuration 0, Pixel Representation 0, Bits Allocated 8 (or larger multiples), Bits Stored up to 38; signedness comes from the SIZ marker.
- **One frame, one codestream:** for multi-frame objects "each Frame shall be encoded separately", the codec runs per frame, and each fragment holds data from a single frame. **The standard never states "independently decodable" in so many words**; independence follows from separate encoding. Whether a frame may span several fragments is not clear: the A.4.4 text I got allows it for Part 1 syntaxes, while a search snippet says the 2025c edition requires exactly one fragment per frame. **Check the current A.4.4 text before writing the encapsulation.** **open**
- **Offset tables:** the first item of the encapsulated Pixel Data is a Basic Offset Table (may be empty); a Basic Offset Table value is not permitted if an Extended Offset Table (7FE0,0001) is present (64-bit offsets; needed for files beyond 4 GB: general knowledge, not checked).
- **Not specified at all:** number of decomposition levels, wavelet-filter parameters beyond "lossless only", tiling and tile-parts inside a frame, progression order. So D = 5, one JPEG 2000 tile per frame, and the choice of progression order are free.

Consequences for this plan: every frame is a raw reversible-5/3 codestream with SGcod = 1 and photometric interpretation YBR_RCT, bits allocated and stored 8, and Pixel Representation 0. Our per-frame independence is exactly what PS3.5 requires. **What this does not verify** is the JPEG 2000 Part 1 convention I assumed (low-pass at even canvas coordinates, whole-sample symmetric extension at the frame edge, vertical-then-horizontal forward order): that lives in ITU-T T.800, which I have not read, and my tests use my own NumPy transforms. OpenJPEG's forward transform was later checked against these conventions (section 8g); supplying our own coefficients to a codec is not possible with standard encoders, so the sizes in 8g encode pixel images that carry the intended coefficients. **partly closed**

## 8d. Colour conversion strategies **[moved]**

Moved to [Color_Conversion_Claude.md](Color_Conversion_Claude.md), section C2. Summary: wavelet-domain conversion is sparse but biased (Cb -4.9 at D = 5); pixel-domain conversion is exact but dense (Cb level-1 non-zero 0.12 % to 63 %); fidelity target is an open decision.

## 8e. Quantizer and deadzone versus EBCOT **[verified in the stored coefficients; real EBCOT sizes not measured]**

Measured per scale over the whole file (`colour_ebcot_probe.py`; detail bands together; "empty LSB" = number of always-zero low planes; zero-order entropy is a crude proxy and counts the padded margin):
- Y: scale 0 (header Quantizer 1, Deadzone 2): 1.55 % non-zero, 1 empty LSB plane, 22.9 % of 64 x 64 blocks all zero; scale 1 to 5: no quantization, 51-64 % non-zero, 2.4-4.0 b/coefficient, almost no empty blocks.
- Co, Cg: scale 0 (4, 16): **0.01 % and 0.00 % non-zero**, 4 empty LSB planes, 96-99 % of blocks all zero; scale 1 (2, 8): 1.6 % and 0.6 % non-zero, 2 empty planes, 47-61 % of blocks all zero; scale 2 to 5: no quantization, 58-69 % non-zero, 3.3-3.6 b/coefficient. The top LL band is never quantized (Y 4.4 b, Co 3.8 b, Cg 3.9 b).
- So the quantizer and deadzone empty out the finest scale (75 % of all coefficients) and part of the next, while the remaining bits sit in Y scales 1 to 5 and chroma scales 2 to 5. Zero-order totals: Y 51.7 MB, Co 16.7, Cg 14.3, i.e. about 83 MB against the 88.9 MB file (Hulsken is a lightly modelled bit-plane coder; EBCOT would be expected to do somewhat better, the white paper says about 10 %: not measured here).
- **Colour transform cost in bits:** Cb + Cr (wavelet domain) costs +10-11 % over Co + Cg (34.1-34.5 MB against 31.1 MB zero-order), about +4 % of the whole file: Cr - Cb = Co, so the pair is correlated. The transform also removes one empty LSB plane from the chroma (scale 0: 4 to 3, scale 1: 2 to 1) because Cb = -Co/2 - Cg halves the quantizer step.
- **What EBCOT gets for free:** all-zero codeblocks (empty packets), and missing most-significant planes (signalled per codeblock). That covers the sparse fine scales of all channels (before any pixel-domain re-transform) completely.
- **Deleting the empty least-significant planes:** see section 8f (T.800 Annex E read). Short version: legal and exact only with a decoder that reconstructs with r = 0, which OpenJPEG does not; the saving is at most a few percent anyway.
- **Deadzone effect:** the smallest non-zero coefficient equals the deadzone value (2, 16, 8), so non-zero chroma coefficients at scales 0 and 1 are isolated multiples of a large step: cheap to locate, but each costs a sign bit.
- Real encoder measurements: section 8g.

## 8f. T.800 Annex E, the midpoint, and dropping empty least-significant planes **[spec read from the ITU PDF; OpenJPEG source read through a summarising tool; measurements ours]**

**What T.800 (08/2002) Annex E says.** Each coefficient is decoded to `N_b(u,v)` of `M_b = G + eps_b - 1` bit-planes. For the reversible transform the step is 1 (E.1.2.1). If all planes are decoded (`N_b = M_b`) the coefficient is recovered exactly (E-7). If `N_b < M_b`, E-8 gives `R = (q +- r * 2^(M_b - N_b)) * 1` for non-zero `q` and 0 for `q = 0`, where **r is "a reconstruction parameter, which can be arbitrarily chosen by the decoder"**, typically 0 <= r < 1 with r = 1/2 common (note under E.1.1.2; it also applies to the reversible case). T.803 adds that r = 1/2 is what conformance testing normally assumes.
So: (1) r = 0 is legal; the standard does not fix the value. (2) With r = 0 and a truncation at a bit-plane boundary, dropped planes reconstruct as zero, i.e. exactly, **if** the dropped bits really are zero.

**What OpenJPEG does** (`t1.c`, as reported by the fetch tool, not verbatim): the decoder works with one extra fractional bit; every newly significant or refined coefficient receives the `half` term of the last decoded plane (`oneplushalf`, `poshalf`), and the reversible output is divided by 2 at the end (truncation toward zero). Decoding to plane 0 leaves no offset (exact); stopping above plane 0 leaves a midpoint offset of +-2^(Q-1) for Q dropped planes. I found no setting for r in the code I read (the public API was not checked), so with OpenJPEG r = 1/2 is effectively fixed.

**Re-evaluation if the midpoint were not required (r = 0, e.g. a patched or different decoder):**
- *Exactness:* dropping the Q always-zero low planes of a code-block becomes lossless **only for code-blocks whose coefficients are all untouched stored values.** Every coefficient recomputed in the halo ring (section 4) is an arbitrary integer, so every code-block that touches a frame border must be coded to plane 0. Ring-free share of code-blocks (64 x 64): F = 256: 0 %; F = 512: 25 % at level 1, 0 % at level 2; F = 1024: 56 % at level 1, 25 % at level 2 (32 x 32 blocks: 25 / 56 / 77 % at level 1 for F = 256 / 512 / 1024). Only J2K levels 1 and 2 carry a quantizer structure at all.
- *Saving:* the planes to skip are Y level 1 (1 plane), wavelet-domain Cb/Cr level 1 (3) and level 2 (1). Counting only ring-free, non-empty code-blocks (all-zero blocks already cost nothing): about 0.47 N skippable coefficient-plane cells for F = 1024 with 64 x 64 blocks (0.19 N for F = 512, 0.69 N for F = 1024 with 32 x 32 blocks), N = 550 M padded pixels. At an assumed cost of 0.01 to 0.05 bit per all-zero cell under adaptive arithmetic coding that was roughly 0.3-1.6 MB for F = 1024 and 0.13-0.64 MB for F = 512, against about 80 MB. **Measured with a real encoder in 8g: the empty planes cost at most 0.6-1.2 % of the stream (upper bound).**
- *Risk:* any decoder that uses r = 1/2 reads such a codestream wrongly. Simulated on a 1024-px tissue frame (`midpoint_probe.py`): 7.6 % of pixels change, Y by up to 3, Cb and Cr by up to 3, RGB by up to 4 (512 px: 2.2 % of pixels, up to 3; 256 px: nothing to skip). The transfer syntax 1.2.840.10008.1.2.4.90 is "Lossless Only", so a stream whose exactness depends on the decoder's choice of r is, in my reading, outside its intent (not a statement from PS3.5).
- *Conclusion:* with stock decoders every plane stays coded, and the cost of that is small (8g: at most about 1 %). Dropping planes would need a non-standard decoder for under about 1 %. **Not recommended.**
- A side finding: the ring-free analysis also says frame size matters for efficiency, since border code-blocks (with the halo) are the ones that cannot exploit the quantizer structure.

## 8g. Real EBCOT sizes (OpenJPEG 2.5.4 through Pillow 13.0) **[verified; scripts `openjpeg_convention_check.py`, `openjpeg_ebcot_sizes.py`, `openjpeg_empty_plane_cost.py`]**

**Conventions confirmed against a real codec.** An image built as the inverse of a sparse coefficient set, using the NumPy transform of every other test here, codes to almost nothing in OpenJPEG (F = 256: 107 non-zero coefficients, 372 bytes; 512: 751 bytes; 1024: 2308 bytes), while the same image shifted by one pixel codes to 3.7 to 83 KB. So OpenJPEG's reversible 5/3 forward transform uses the conventions assumed throughout (low-pass at even local coordinates, mirror at the frame edge, vertical then horizontal). All round trips were lossless, and `mct=1` really applies the reversible colour transform. (This confirms OpenJPEG, not the T.800 text.)

**Sizes** (2048 x 2048 level-0 regions; each frame an independent raw codestream, reversible 5/3, RCT, D = 5, 64 x 64 code-blocks; routes feed OpenJPEG the pixels that carry exactly the intended YBR_RCT coefficients; iSyntax size = area-weighted compressed codeblocks of the region):

| region | iSyntax | route | 256 px frames | 512 | 1024 | 2048 |
|---|---|---|---|---|---|---|
| tissue | 1.193 MB | wavelet-domain (W-floor) | 0.99x | 0.95x | 0.93x | 0.92x |
| | | exact pixel-domain, floor (E-floor) | 2.14x | 2.10x | 2.09x | 2.09x |
| | | exact pixel-domain, fix (E-fix) | 2.66x | 2.63x | 2.62x | 2.61x |
| tissue edge | 1.109 MB | W-floor / E-floor / E-fix at 1024 | 0.91x / 2.10x / 2.61x | | | |
| scanned glass (Y about 230, inside the envelope, not filler) | 0.453 MB | W-floor / E-floor / E-fix at 1024 | 0.88x / 3.12x / 4.02x | | | |

- The wavelet-domain route costs about the same as iSyntax (0.87-0.99x); larger frames help (tissue: -4.5 % at 512, -6.3 % at 1024, -7.2 % at 2048 relative to 256). This matches the white paper's "about 10 % better than Hulsken" only roughly (0.88-0.93x at 1024 and up).
- The exact pixel-domain route is **2.1 to 3.2 times larger** (floor), **2.6 to 4.0 times** (fix): the chroma sparsity loss of C2 (Color_Conversion_Claude.md), now measured with a real encoder.
- 32 x 32 code-blocks add 3 to 4 % over 64 x 64.
- Clipping: an 8-bit codestream cannot hold pixels outside 0..255, so those were clipped (tissue: 1.4-2.5 k of 4.2 M pixels; blank glass: none). This is inherent to 8-bit DICOM and changes the image there.
- The whole-slide file would be roughly 0.9-1.0x of 88.9 MB for the base layer in the wavelet-domain route (extrapolated from three regions). **A complete DICOM pyramid is far larger than the base layer alone:** each layer is its own image, and layer k carries the bits of iSyntax scales k and coarser, which hold most of the data (zero-order estimate from 8e: layers 0..5 hold about 1.21, 1.11, 0.62, 0.20, 0.05, 0.01 bits per padded pixel, about 2.6x the iSyntax size for all six layers; a coarser set of layers costs less). Rough, not encoded. **open**

**Cost of the empty least-significant planes, measured** (planes removed by dividing the stored multiples by 2^Q, which leaves the earlier planes' coding unchanged; applied also to border code-blocks, so an upper bound for what plane deletion could save): tissue 0.79 % (1024 px) to 1.18 % (256 px) of the stream; blank glass 0.34-0.56 %; Y planes account for about two thirds. With r = 1/2 decoders (OpenJPEG) the planes must be coded anyway; even with an r = 0 decoder the gain would be under about 1 %. Not worth a non-standard decoder (section 8f).

Net effect on the open decision (Color_Conversion_Claude.md C2): the wavelet-domain conversion keeps the file at iSyntax size but biases chroma (near-lossless at best with a per-tile correction); exact pixel-domain conversion is exact against a chosen convention at 2 to 4 times the size.

## 8h. Omitting tiles: size savings and review material **[verified: real OpenJPEG sizes of every tile; scripts `omission_size_estimate.py`, `omission_summary.py`, `omission_review_images.py`; files in `$SOURCEROOT/transcoding_review/`]**

All six layers, all tiles, wavelet-domain colour route (W-floor), reversible 5/3, RCT, D = 5, 64 x 64 code-blocks, frames as independent codestreams (bytes of the codestreams only). Whole pyramid encoded: **204.4 MB at 256 px frames, 198.6 MB at 512, 196.4 MB at 1024** (layer 0 alone: 84.0 / 79.4 / 77.7 MB, against 88.9 MB for the iSyntax file; the 6-layer pyramid is 2.3x, below the 2.6x zero-order estimate of 8g).

Two strategies kept in play (definitions in `omission_summary.py`):
- **A: omit only exact filler tiles** (YCoCg = 255, 0, 0 everywhere; the right and bottom margin bands). A filler frame encodes to about 170 bytes, so omission saves **0.08 MB (0.04 %) at 256 px, 0.01 MB at 512, nothing at 1024** (506 / 49 / 0 frames over all layers; the per-frame DICOM metadata saved is perhaps 0.1 MB at an assumed 250 B/frame). A is about tidiness and frame count, not size, and changes no pixel.
- **B: also omit scanned-glass tiles that are nearly constant.** "Nearly constant" has two readings, which behave very differently:
  - *Robust range* (99.5 % minus 0.5 % percentile of RGB, worst channel) <= T: T = 3: 3.0 %; **T = 4: 23.1 % (47.1 MB, 5627 frames)**; T = 5: 32.3 %; T = 8: 33.7 %; T = 16: 34.1 % of the 6-layer pyramid at 256 px frames (512 px: 1.5 / 16.1 / 24.2 / 26.3 / 26.8 %; 1024 px: 0.5 / 10.8 / 15.3 / 16.9 / 17.1 %). The step between T = 3 and 5 is the sensor-noise range of glass (about 4 levels). At layer 0 and 256 px, 4844 tiles (35 % of the layer's bytes, about 6 KB = 0.75 bit/px each) are omitted at T = 4: scanned glass is a large share of the file because its noise is stored losslessly.
  - *Strict range* (max - min) <= S, i.e. no speck: S = 8: 9.5 %; S = 16: 13.1 %; S = 32: 21.6 %; S = 64: 34.6 % (512 px: 1.0 / 3.0 / 7.2 / 23.3 %; 1024 px: 0 / 0.2 / 1.0 / 8.7 %).
  - The difference is dust: tiles omitted by the robust rule have strict ranges up to 167 levels (median 8, 95th percentile 59) in layer 0, so the robust rule deletes isolated specks (visible as dots in `diff_R4_x16.png` and the crops).
  - Savings fall with frame size because large frames contain tissue: at T = 4, 23 % / 16 % / 11 % for 256 / 512 / 1024 px. Layer 0 and layer 1 carry most of the saving (T = 4, 256 px: 35 % and 25 % of those layers; layer 2: 0.1 %, T = 8: 30 %; layers 3 to 5: small).
  - **Speck size (measured 2026-11-26, `speck_size_analysis.py`, `transcoding_review/speck_components_T4_F256.csv`).** Layer 0, 256-px frames, the 4387 non-filler tiles omitted at T = 4; pixel = 0.25 um. Specks = 8-connected pixels deviating from the tile's median (= fill) colour by more than 8 (16) levels in any channel (glass noise is about 4). 1563 of the 4387 tiles (35.6 %) contain at least one speck (2036 specks; 29 % of tiles and 1596 specks at > 16 levels). Area: median 12 px, p99 98 px, max 199 px (about 12 um2). Equivalent diameter: median 1.0 um, p99 2.8 um, **max 4.0 um**; longest bounding-box side: median 1.0 um, p99 3.7 um, max 6.5 um (the largest are elongated, 26 x 10 px, fibre-like). 97 % are under 2 um, all under 5 um. This is bounded by the rule itself: at most 0.5 % of a 256-px tile (328 px, about 20 um2, 5 um diameter) per tail can lie outside the robust range. Largest single-pixel deviation per tile: median 5, p95 55, max 163 levels. Not measured: layers 1-5 (the same pixel bound is 4x larger in area per layer step, about 10 um at layer 1), other frame sizes, other slides, and what the specks are (dust, debris or stain; the data cannot say).
- **Fill colour:** the median colour of the omitted glass tiles is RGB (228, 230, 233); their tile medians lie within 1 level (95th percentile) and 2 levels (max) of it at T = 4 (up to 3 to 6 levels at T = 8 to 16). A single Recommended Absent Pixel colour can match either the glass or the pure-white filler bands, not both (the bands lie outside the scan, so showing them as glass colour is arguably natural, but that is a visual-review item).
- B is lossy: omitted tiles lose the sensor noise (about 4 levels) and any specks. Visual review material is in `transcoding_review/` (README.txt there): masks for A, R4 and S16, single-fill previews, x16 difference images and full-resolution crops of borderline and speck cases.
- Not included: DICOM per-frame functional-group metadata (assumed 250 B/frame, not measured), the exact-pixel colour route (about 2 to 4 times larger, 8g), and any decision on which tiles are acceptable to omit.

## 8i. Target RGB instead of YBR_RCT **[moved]**

Moved to [Color_Conversion_Claude.md](Color_Conversion_Claude.md), section C3. Summary: every RGB channel carries Y; exact route 3.0x (YBR_RCT 2.1x) the iSyntax size, wavelet-domain route 1.6-1.8x and biased in all channels.

## 9. Tools (`tools/v1_analysis/`, see its README.txt)

`halo_1d.py`, `frame_interior_test.py`, `strip_halo.py`, `windowed_idwt_halo.py`, `dicom_tile_tests.py`, `colour_ebcot_probe.py`, `colour_bias_by_levels.py`, `colour_sparsity_cost.py`, `colour_dc_fix.py`, `rgb_route_probe.py`, `midpoint_probe.py`, `openjpeg_convention_check.py`, `openjpeg_ebcot_sizes.py`, `openjpeg_empty_plane_cost.py`, `omission_size_estimate.py`, `omission_summary.py`, `omission_review_images.py`, `count_chunks.py`, `rct_ycocg_vs_rct.py`, `ll_probe.py`. All read-only; the decoder is untouched.
