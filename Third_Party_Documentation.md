# Third-party documentation: review

Companion to [CLAUDE.md](CLAUDE.md) and [Example_FIC_Files_Claude.md](Example_FIC_Files_Claude.md). It reviews everything in `../third_party_documentation/` and says what each document gives you, how far to trust it, and where it disagrees with the code or the sample headers.

**Evidence labels:** **[doc]** = stated in the document; **[code]** = checked in this repo (`file:line`); **[FIC]** = checked in the sample headers; **[inferred]** = my deduction.

## 0. How the review was done, and its limits

- **No PDF tools are installed**, and installing one was denied. I extracted text from the four text-based PDFs with a stdlib-only Python script (zlib plus text operators). Prose, tables and MATLAB listings came through. **Figures, and the equations in the format spec, did not**, or came through scrambled. Page counts below are approximate.
- **`HulskenPatentUS20180255321A1.pdf` is unreadable by me.** It is a 12-page scan stored as CCITT fax images (made by `tiff2pdf`, 2018-08-21) with no text layer, and I have no OCR. Nothing in this file is reported from the patent itself. If you want it analysed, supply the text (for example from Google Patents for US 2018/0255321 A1).
- **I did not test the AI-assisted notes' mathematics.** I checked only the claims that touch this repo's code (section 5).
- **Licence note.** The Philips documents grant copying and redistribution **of the whole document only**, not parts. Quote them sparingly, and don't copy excerpts into the repo.

## 1. Inventory

| File | Type | Size | What it is | Trust |
|---|---|---|---|---|
| `2020_04_24_Pathology_iSyntax_image-format.pdf` | Philips spec (doc 4555 207 43941) | ~33 pp | **File format spec**: data model, XML header, codeblocks, seektable, padding, reading label/macro/WSI | Authoritative, but written for the v1 header (section 3.1) |
| `2020_04_24_30pp_..._free-online-conversion.md` | Markdown conversion of the above | 786 lines | Same text. Equations and figures are lost. | Same, plus conversion loss |
| `4522 207 43951_2020_04_24 Pathology iSyntax image viewing.pdf` | Philips spec | 22 pp | **Display post-processing**: CLAHE contrast enhancement and sharpening, with MATLAB code. Also gives the YCoCg-R formulas. | Authoritative |
| `Hulsken2016_Philips-iSyntax-for-Digital-Pathology.pdf` | Philips white paper | 4 pp | Why iSyntax exists; the compression idea; the cluster/tile storage model | Authoritative, but marketing level |
| `DP-200666_iSynta_Hulsken_Webinar.pdf` | Philips slides, 2020-07-01 | 13 slides | The compression algorithm in six bullet points; benchmarks | Authoritative on the algorithm, vendor claims on speed |
| `HulskenPatentUS20180255321A1.pdf` | US patent application | 12 pp | Presumably the Hulsken entropy-coding method. **Not read.** | n/a |
| `YCoCg-R_to_YBR_RCT_Transform_Claude.md` | AI chat export | 387 lines | Colour matrix YCoCg-R → DICOM YBR_RCT in the wavelet domain | Unverified (section 5) |
| `iSyntax_codeblock_distribution_reconstruction_order_Claude.md` | AI chat export | 2,855 lines | Largest note. Code review, GPU design, DC offset, review of the Philips Python and MATLAB reference code, 5/3 filter match, chunk layout. | Unverified; several self-corrections |
| `iSyntax_to_DICOM_Transcoding_Claude.md` | AI chat export | 459 lines | Feasibility of a wavelet-domain iSyntax → DICOM JPEG 2000 path | Unverified; first answers made without the code |
| `iSyntax_to_JP2000_wavelet_transcoding_Claude.md` | AI chat export | 1,285 lines | Edge-coefficient correction, lossy effects, AWS cost, file size, DICOM .91 | Unverified; contains models and guesses |
| `iSyntax_to_JP2000_wavelet_transcoding_Claude (1).md` | AI chat export | 635 lines | An **earlier, shorter** run of the same topic (halo convolution, "boundary-only" correction) | Unverified; superseded by the file above |
| `iSyntax_vs_OME_XML_Comparison_Claude.md` | AI chat export | 555 lines | iSyntax XML vs OME-XML, `tifftools`, `tifffile` | General knowledge only; no sample files seen |

The `*_Claude.md` files are exports of conversations, so every answer appears 2–3 times in them. Read the first copy only.

## 2. The Philips documents

### 2.1 Format specification (the core reference)

Everything below is **[doc]** unless marked.

- **Data model.** Three images per file: label, macro, WSI. Plus acquisition metadata (DICOM attributes).
- **File layout.** `XML header | EOT 0D 0A 04 | [pixel-data tag 0x7FE0,0x0010 + length, usually 0xFFFFFFFF] | [seektable] | codeblocks`.
  - The seektable starts with tag `0x301D,0x2015` and holds, per codeblock, offset `0x2010` (u64) and size `0x2011` (u64).
  - Each codeblock is a DICOM item `0xFFFE,0xE000` + u32 size + payload.
  - All values are little-endian.
- **XML grammar.** Three node types: leaf `Attribute` (Name, Group, Element, PMSVR), branch `DataObject` (ObjectType), array `Array`. Attribute names starting `DICOM_` are standard tags. The rest are Philips private tags in groups `0x101D` and `0x301D`.
- **Codeblocks.** After the recursive DWT, HL, LH and HH are "spatially aligned and merged into one block [HL, LH, HH]". LL and [HL, LH, HH] are cut into codeblocks per colour channel. Each is Hulsken-compressed. Typical size is 128×128, and the size is read from the header (it "may vary from scanner to scanner").
- **Quantizer and deadzone** (Tables 11 and 12):
  - "Typical 2 objects: first for luminance, second for the other color components."
  - Levels listed from level 0.
  - Quantizer: coefficients are "**rounded to 2^quantizer**".
  - Deadzone: coefficients are "pulled to zero for absolute values **greater or equal** to dead zone". This wording looks like a typo for "less than".
  - Interpretation, supported by the sample files, is in section 3.4 and in [Example_FIC_Files_Claude.md](Example_FIC_Files_Claude.md) section 6.
- **Compression method** (Table 8): `16` = Hulsken compression, `19` = Hulsken2. Matches `isyntax.c:727-732` **[code]**.
- **Padding.** Uniform on all sides, `PerLevelPadding × (2^levels − 1)`, with PerLevelPadding = 3 for LeGall 5/3. Then more on the right and bottom, so the base image is a multiple of the codeblock size at every level.
  - Equations 1 and 2 are scrambled in both my extraction and the markdown. The formula `3·(2^L − 1)` **[inferred]** is confirmed by the code and by every sample's valid-data envelope **[FIC]**.
  - The grid-size formula appears to be `ceil(width / (BlockWidth·2^L)) · 2^(L−1)`. The AI note quotes the same from the Philips Python reference.
- **Block ID** is counted left to right, top to bottom, "from the base level to the top level and finally from the first to the last color channel". Section 5 of the spec says the codeblocks are *written* "from top level to base level". Both statements are in the spec, and they describe ID order and file order separately.
- **Missing blocks.** Tissue gaps are not written. Their block-header coordinates are zero. The number of codeblocks addressed can exceed the number stored.
- **Reconstruction recipe.** Hulsken-decompress; stitch spatially related blocks; inverse DWT per level (the lifting scheme is "`rbio2.2` (legall5/3)"); convert YCoCg to RGB. The Matlab reference calls are `hulskendecompress(…, 16, [128 128], 1)`, `wavelet2dilift`, `ycocg2rgb`, `reconstruct_image`.
- **Reference code.** The spec refers to Python samples (`extract_macro_label_image.py`, `extract_codeblocks.py`) and Matlab/Octave files from the Open Pathology Portal. That code is not in `third_party_documentation/`, but it is now available and reviewed: the Python and the `reconstruct_image` MATLAB in [Image_Format_Code_Samples_Claude.md](Image_Format_Code_Samples_Claude.md), and the codec in [Philips_Compression_Code_Samples_Claude.md](Philips_Compression_Code_Samples_Claude.md). The spec's recipe matches that code, but **that code only reads the v1 layout**, not the v2 headers in the sample FICs.

### 2.2 Image-viewing specification

It covers step 3 of the pipeline (display post-processing), and it is not about decoding. Still, it matters for pixel equivalence with Philips' viewer.

- **Pipeline:** `enhanced = sharpening(applyCLAHE(image), level)`. Philips states that following these steps gives a display "pixel equivalent" to the IMS viewer for UFS-scanned slides.
- **CLAHE:**
  - Luma only, in YCoCg-R space.
  - 256 bins, contextual regions of **40×40** px, clip limit **1.2**, bilinear histogram interpolation.
  - Region size is independent of level. The requested patch is padded with a 40 px frame.
- **Sharpening:**
  - 5×5 integer kernel `[0 0 1 0 0; 0 5 15 5 0; 1 15 44 15 1; 0 5 15 5 0; 0 0 1 0 0] / 128`.
  - Gain 2.0, **halved per level above the base**. Two border pixels are discarded.
- **YCoCg-R formulas** (Figures 8 and 14), the vendor's own definition **[doc]**:
  ```
  forward: Co = R − B;  tmp = B + (Co >> 1);  Cg = G − tmp;  Y = tmp + (Cg >> 1)
  inverse: tmp = Y − (Cg >> 1);  G = Cg + tmp;  B = tmp − (Co >> 1);  R = B + Co
  ```
  The listing uses `bitshift(·, -1)`, i.e. a right shift (floor behaviour on negatives, *if* it is arithmetic). Section 5.1 explains why this matters.
- **Header attributes that belong here** are `DP_IMAGE_POST_PROCESSING` (0x1014), `DP_SHARPNESS_GAIN_RGB24`, `DP_CLAHE_CLIP_LIMIT_Y16`, `DP_CLAHE_NR_BINS_Y16`, `DP_CLAHE_CONTEXT_DIMENSION_Y16` (0x1015–0x1018). They appear as stubs in the code (`isyntax.c:624-628`). **Neither the format spec nor any sample FIC contains them**, so the values for these parameters are not available to us. The Philips reference code hard-codes sharpening gain 2.0, clip factor 1.2, 256 bins and region size 40 (see [Viewing_Code_Samples_Claude.md](Viewing_Code_Samples_Claude.md) §2). The listings in the PDF also differ slightly from the `.m` files (see §7 there).
- libisyntax does **no** CLAHE or sharpening **[code]**, so its output will not match the Philips viewer pixel for pixel.

### 2.3 Hulsken 2016 white paper (4 pages, text extracted)

- iSyntax uses "the same recursive wavelet transformations as JPEG 2000" but replaces EBCOT with "a far simpler entropy coder, based on exploiting correlation of data in a local neighborhood combined with a fast modified Huffman coder".
- Cost: about **10% larger** than JPEG 2000. Gain: "an order of magnitude faster" compression and decompression.
- Figure 2 (not extracted) is described: DWT, then per-coefficient-type bit-plane serialisation, then entropy coding "in groups of **4x2 neighboring pixels**", stored "in order of significance".
- **Storage model:** a three-level hierarchy of wavelet-transformed image, **clusters** (regions spanning space and scale), and **tiles** (single scale, small). Clusters are for scanner streaming, and tiles for viewer random access.
  - This is exactly the cluster/codeblock structure in the v2 headers (Example_FIC_Files_Claude.md section 5.2).
- Other claims: whole slides can be stored across tiers by resolution level (the paper's example: 40× scan with levels up to 5× on fast storage gives a 64× saving). Data signing uses MD5, and encryption is at rest and in transit. These are not visible in the header or in libisyntax.

### 2.4 Webinar (DP-200666, 2020-07-01)

The most useful short statement of the algorithm, slides 6–8:

1. Wavelet coefficients are cut into blocks for random access and converted to bit-planar form.
2. Each bit plane of each block is split into **4×2 pixel areas**.
3. Areas are serialized in a **snake order**.
4. Each 4×2 area of 1-bit pixels becomes one symbol in **[0, 255]**.
5. Symbols are Huffman-coded, and **zero runs get a zero-run symbol plus a length count**.
6. Within a block, coefficient bits are stored **most significant to least significant**, which gives progressive resolution and quality.

The "only novelty" is generalising constant-area coding from three states (all black, all white, other) to all 256 patterns of a 4×2 area.

Benchmarks (vendor claims): within 3.4% of JPEG 2000 size, 17.7% smaller than PNG; compression 250% faster than Kakadu; decoding 79% faster than OpenJPEG. Reference code is described as publicly available and free to use.

**How this clarifies the code.** One 4×2 area is 8 bits, which is **one byte** of a bitplane in `isyntax_hulsken_decompress`. The byte's 8 bits map to 8 consecutive coefficients, and the reshuffle that follows (`isyntax.c:3014-3036`) places 16 consecutive values (two bytes, i.e. two 4×2 areas) into a 4×4 area **[code]**. So the earlier note in `CLAUDE.md` ("comment says 4x2 snake but the code copies a plain 4×4 tile") was misleading. The code's 4×4 reshuffle is two stacked 4×2 areas, which is consistent with the webinar. Whether the area order really follows a snake path remains **[inferred]**: the code does not reverse any rows.

**Order of planes.** The webinar says MSB first. Compressor **v2** stores sign then MSB→LSB, matching it. Compressor **v1** stores sign then LSB→MSB (`isyntax.c:2975-2977`). All six sample files use compressor v1, so the webinar's order describes a format the samples do not use **[FIC]**.

### 2.5 Patent

Not read (section 0). It likely contains the formal claims on the constant-area bit-plane Huffman scheme and the hierarchical data model ("patent pending" in the 2016 paper). Useful when you need the *intended* definition of the Huffman tree serialization, the zero-run convention, and any quantization or deadzone language.

## 3. What the Philips documents say about the four goals

### 3.1 File format: spec versus the six samples

The spec describes the **v1 style header**. All six sample FICs are data model 100.4 (v2). They differ from the spec as follows, all **[doc]** versus **[FIC]**:

| Item | Spec says | Samples show |
|---|---|---|
| Block count and layout | Block header table 0x2014 plus seektable | Absent. A **cluster header table** 0x201f (binary, intact) and **cluster header templates** 0x2016. Neither is in the spec. |
| Valid-data envelopes (0x2023–0x2026), dimension origin 0x2027, cluster dims 0x2017/0x2021 | Not mentioned | Present in every file |
| `DICOM_BITS_STORED` / `HIGH_BIT` | "iSyntax: fixed value 16" / "15" | **8** / **7** |
| Wavelet coefficient names | `["LL" "LH" "HL" "HH"]` | `"HL" "LH" "HH" "LL"` |
| Block template waveletcoef range | LL `[0 0 0]`, otherwise `[1 1 3]` | LL `3 1 3`, H group `0 1 2` |
| Block template scale range | `[dwt_level 0 dwt_level]` | `s 1 s` (step 1) |
| `DIMENSIONS_IN_BLOCK` | `[1 0 4]` (y x waveletcoef) | `0 1 4` (x y waveletcoef) |
| Colour ranges | `[1 0 1]`, `[2 0 2]` | `1 1 1`, `2 1 2` |
| Block coordinate `coefficient` | 0 for LL, 1 for HL/LH/HH | Cluster templates use waveletcoef 3 for LL; the code maps it to `coefficient 0` (`isyntax.c:970`) |
| Lossy ratio and method | Not specified | 15 / `PHILIPS_DP_1_0` (WSI); 26 / `ISO_10918_1` (label, macro) |

Treat the spec as the v1 reference and the code as the only description of the v2 header. The notes in `Example_FIC_Files_Claude.md` document v2 from the samples.

One more point **[doc]**: the spec says the iSyntax lossy colour values are `["Y" "Co" "Cg"]` and lossless are `["R" "G" "B"]`. Only lossy files exist in the samples.

### 3.2 5/3 wavelet

- **Spec:** the lifting scheme is `rbio2.2` (LeGall 5/3), and the vendor MATLAB uses `liftwave("rbio2.2","integer2integer")` (per the AI note about `idwt_at_level.m`). The base level is padded 3 px per level.
- **Code:** `isyntax_dwt.c` is OpenJPEG's reversible 5/3 inverse, adapted **[code]**. The AI note's finding that it is "bit for bit" the JPEG 2000 reversible 5/3 inverse is consistent with the code, but the note did not run it.
- **The code always calls it with `cas = 1`** (`isyntax.c:1841`, `:1861`). The first half of each row or column (LL/LH or LL/HL) is therefore at **odd** positions, which is not the usual JPEG 2000 tile convention. The edge-correction formulas in the AI notes assume an **even-origin (`cas = 0`) tile**. That mismatch is not addressed in any of the notes. It matters if you transcode to JPEG 2000 (section 5.3).
- 16-bit coefficients: the SIMD paths wrap at 16 bits while the scalar vertical paths widen to 32 bits. This is raised in the AI notes as a caveat for inputs near ±32767. On the real file the largest intermediate magnitude over all levels is 282 ([Example_iSyntax_File_Claude.md](Example_iSyntax_File_Claude.md) §9), so the caveat does not arise for ordinary 8-bit slides.

### 3.3 Hulsken compression

The authoritative statements are the webinar (4×2 areas, snake order, Huffman plus zero-run, MSB first) and the white paper (10% larger than JPEG 2000, much faster). **Nothing in the readable documents gives the bit-level format** (tree serialization, header, bitmask layout, v1 vs v2 differences). That definition comes only from the code (`isyntax.c:2551`) and from the Philips MATLAB reference (`hulskendecompress`), which isn't in the folder. The patent may contain it.

### 3.4 Quantizer and deadzone

- Defined in spec Tables 11 and 12 (section 2.1).
- The viewing document says nothing about them.
- The webinar's "bitplane compression supports any bit depth" and "MSB first" are compatible with the idea that quantization drops low planes. They do not state it.
- **Update:** Philips' `quantanddeadzone.m` now defines both parameters exactly (deadzone `abs(c) <= dz → 0`, then `fix(c/2^Q)*2^Q`; see Philips_Compression_Code_Samples_Claude.md section 3), so the spec's "greater or equal" wording is confirmed wrong.
- **Do not rely on the AI note "Effect of iSyntax Q2-lossy"** (in `iSyntax_to_JP2000_wavelet_transcoding_Claude.md`, lines 461-546). It models the decoder as `ĉ = sign(q)·⌊(|q|+r)·Δ_b⌉` with a dequantization step and a reconstruction offset `r`, and states that JPEG 2000's deadzone quantizer "is the same as iSyntax's". Neither claim comes from a Philips document. The spec's "rounded to 2^quantizer" and libisyntax's lack of any dequantization **[code]** point to a different model: low bits zeroed in place, no scale factor. See [Example_FIC_Files_Claude.md](Example_FIC_Files_Claude.md) section 6.
- The DICOM lossy-flag advice in that note (0028,2110 = 01, ratio 15) is sound, and consistent with what the header already carries.

## 4. AI-assisted notes

These record earlier work and indicate your goal: **transcode iSyntax to DICOM WSI with JPEG 2000 lossless (5/3), staying in the wavelet domain**, possibly on GPU. They are chat exports, not references, and they include several retractions. Use them as a to-do list and a source of hypotheses.

### 4.1 Reliability problems visible in the notes themselves

- Several first answers were written **without the attached files** and later retracted. Examples: the "seven bugs" in `iSyntax_to_DICOM_Transcoding_Claude.md` (line ~235: "Please disregard them"), and a metadata summary in the codeblock note (line ~1271: "I made that summary up").
- The DICOM note and the OME note say they work "from memory of the library". The OME note never saw the sample headers.
- Some claims are later contradicted by the same file. Two examples:
  - The codeblock note says a wavelet-domain transcode "will not work" because DICOM frames need their own symmetric extension (§3 of the answer near line 2290).
  - The JPEG 2000 notes then describe how to do it with boundary-only corrections.

### 4.2 Content summary

| Note | Main content | Verified against this repo |
|---|---|---|
| `YCoCg-R_to_YBR_RCT` | 3×3 matrix `Y'=Y, Cb=−Co/2−Cg, Cr=Co/2−Cg`; Y differs from RCT luma by δ ∈ {−1,0}; transform not exactly commuting with lifting floors; lists code issues | Partly: see 5.1. Its claim that the scalar path "clamps negatives" is **wrong**. |
| `codeblock_distribution…` | Review of `isyntax_reader.c`/`isyntax_dwt.c`; GPU tile-parallel design; **no DC level shift is applied**; chunk = 1+4+16 H blocks per colour for 3 scales, index `0 / 1+(y%2)*2+(x%2) / 5+(y%4)*4+(x%4)`; Philips Python index formula; encoder-side MATLAB findings; "5/3 filter matches J2K" | Partly: no DC shift **[code]**; the 21-block layout matches the sample clusters **[FIC]** (16+4+1); the chunk index formula not checked |
| `iSyntax_to_DICOM` | Feasibility verdict: libisyntax is a good front end (parse, schedule, Hulsken decode), but its public API returns pixels only. A fork must expose coefficients. Edge-correction cascades across levels. | Consistent with `libisyntax.h` (pixels only) |
| `iSyntax_to_JP2000…` | 1-D/2-D edge-correction formulas; transpose identity (H-then-V vs V-then-H); lossy effect; **AWS cost and file-size models**; DICOM .91 with 5/3 | Not verified. The cost and size figures are models with stated 2-3× uncertainty. |
| `…Claude (1).md` | Earlier version of the edge-correction topic | Superseded |
| `iSyntax_vs_OME_XML` | Generic mapping of iSyntax to OME-XML; `tifftools` vs `tifffile` | General knowledge |

### 4.3 Facts from the notes that are worth keeping

All **[AI note]** unless a check is given.

- **Chunk layout.** One chunk holds three consecutive scales per colour, as 1 + 4 + 16 H blocks (plus the top LL only at the highest scale). **[FIC]**: of the 17 cluster templates across the six files, 11 have 63 headers (21 per colour), 3 have 66 (three scales plus LL), and 3 have 18 (a top cluster of two scales plus LL).
- **Row-major, H then V.** The IDWT buffer is row-major and the horizontal pass runs first **[code]** (`isyntax_idwt`, `isyntax.c:1843-1870`). Confirmed independently.
- **Subband order in a block:** HL, LH, HH **[code]** (`get_offsetted_coeff_blocks`).
- **3×3 neighbour requirement:** a tile's IDWT needs LL and H of all eight neighbours (or dummy coefficients if absent) **[code]**.
- **Philips Python indexing** (now read directly; see [Image_Format_Code_Samples_Claude.md](Image_Format_Code_Samples_Claude.md) §3.2): `offset = (3<<s) − 2` for H and `+ (1<<s)` for LL give the per-scale origin. The C block-ID code is a port of it (`isyntax.c:3300-3330`). This matches `get_first_valid_coef_pixel` / `get_first_valid_ll_pixel` **[code]** (`isyntax.c:3047-3055`), and the cluster range starts in the sample FICs **[FIC]**.
- **Hulsken v2 bitplane seektable.** The encoder-side MATLAB reportedly does not split zero runs at plane boundaries, so the recorded plane offsets may point past a plane's start. libisyntax reads these offsets into `bitplane_offsets[]` but **never uses them** **[code]** (`isyntax.c:2628-2641`), which is consistent with that advice not to depend on them.
- **DC offset:** none. Y white fallback is 255 and chroma is 0 **[code]**.
- **Colour:** the note states that, with a reversible YCoCg-R, Y is bit-identical to the RCT luma (`⌊(R+2G+B)/4⌋`), and that only chroma needs the parity term. The vendor viewing document supplies the YCoCg-R definition the note assumes **[doc]**.

## 5. Findings from checking the notes against this repo

### 5.1 YCoCg-R: rounding and clamping in `isyntax.c`

- **Which path runs depends on the build:** CMake uses `-march=native` (SSSE3, floor path); the Meson build on x86 has no `-march` flag (only SSE2), so it runs the scalar truncating path for every pixel (see CLAUDE.md, "Build / run").
- **Measured effect (real file, [Example_iSyntax_File_Claude.md](Example_iSyntax_File_Claude.md) §9.1):** floor vs truncation changes **64.5% of pixels** (R/B by 0-2, G by 0-1 levels), not just tail pixels; Y never negative and no RGB below 0 on that slide.
- **Resolved: Philips has two versions.** The viewing PDF and `ViewingCodeSamples/*.m` use `bitshift(x,-1)`; `CompressionCodeSamples/rgb2ycocg.m` and `ycocg2rgb.m` use `fix(x/2)` (truncation toward zero). The scalar C path matches the compression samples, and the SIMD path matches the viewing samples. See Philips_Compression_Code_Samples_Claude.md section 6. How MATLAB/Octave `bitshift` treats negative `int16` is unverified.
- **The two code paths differ [code].** The SIMD paths use arithmetic shifts (`_mm_srai_epi16`, `vshrq_n_s16`), which floor (`isyntax.c:1693-1695`, `:1720`). The scalar fallbacks use `Cg/2` and `Co/2`, which truncate toward zero (`isyntax.c:1664-1667`, `:1672-1675`). They differ for negative odd Co or Cg.
- **The scalar path does not clamp negatives [code].** It applies `ATMOST(255, x)` (a `MIN`) and stores into `u8`. A negative result wraps modulo 256 instead of becoming 0, whereas the SIMD paths saturate with `packus`/`vqmovun`. The YCoCg note's statement that "scalar clamps negatives then caps at 255" is incorrect. The scalar path is used for the tail pixels and for any build without SSSE3 (`isyntax.c:1679`). Tile width is always a multiple of 8, so the tail is rarely used on x86.
- **Y is folded with an absolute value** after the IDWT (`isyntax.c:2399`), only for the output tile and not for the LL passed to child tiles. Confirmed **[code]**. It is nonlinear, so it has no wavelet-domain equivalent. Whether the vendor viewer does the same is not stated in the documents.

### 5.2 Possible bugs raised in the notes

- **`isyntax_is_parent_tile_missing` has swapped arguments.** The signature is `(wsi, scale, tile_x, tile_y)` (`isyntax.c:1959`), but all eight calls pass `(…, tile_y±k, tile_x±k)` (`isyntax.c:2038-2221`). It is symmetric only for tiles on the diagonal. It affects whether a neighbour's missing LL is treated as invalid, not the colour or wavelet maths. **Verified, upstream issue, still present in this fork.** Not fixed, per the "do not change decoder behavior" rule in `CLAUDE.md`.
- **The `isyntax_copy_cstring` off-by-one** reported in a note was already fixed by the last upstream merge (commit `3146a9a`, `92471db`) **[git]**. The notes were written against an older copy.

### 5.3 The cas = 1 parity issue (for the JPEG 2000 plan)

The JPEG 2000 edge-correction notes derive their formulas for a tile whose origin is an even coordinate, with the low band at even positions (`cas = 0`). libisyntax inverts with `cas = 1` (section 3.2). One of the notes raised "whether cas is ever 1" as a caveat and then did not close it. Before using those formulas, work out the mapping between iSyntax's parity (`first valid pixel` = `(3<<s) − 2` for H, `+ (1<<s)` for LL) and JPEG 2000's subband origins, and test the result against the bit-exact oracle in `isyntax_dwt.c`.

### 5.4 Samples versus the notes

- The codeblock note describes the Philips reference file `1.txt` as "the smallest reference, **version-1**, iSyntax WSI file". The sample `1.fic` has `PIM_DP_UFS_INTERFACE_VERSION = 100.4`, which libisyntax classes as **v2** (data model ≥ 100) **[FIC]**. I ran the Philips Python metadata reader on copies of `1.fic` and `PAS.fic`: it reads the dimensions correctly but fails at `block_header_table()` because there is no `UFS_IMAGE_BLOCK_HEADER_TABLE` (see [Image_Format_Code_Samples_Claude.md](Image_Format_Code_Samples_Claude.md) §3.1). So either that SDK file is not the same as `1.fic`, or "version-1" meant compressor v1 (true of all six samples). Resolved: `example_isyntax_files/1.isyntax` is the v1 (5.0) file and `1.fic` is the v2 export of the same scan; see [Example_iSyntax_File_Claude.md](Example_iSyntax_File_Claude.md) §1.
- The "Q2-lossy" note says the *v1* schema achieves 15:1. The samples show ratio 15 as the declared value, and observed ratios of 11.6–64 depending on tissue content (Example_FIC_Files_Claude.md section 1).

## 6. Gaps: what is referenced but not in the folder

| Missing | Why it matters | Where it is mentioned |
|---|---|---|
| ~~Philips **MATLAB reference code**~~ **Now available and reviewed** at `/mnt/d/Whole_Slide/iSyntax_Reference_Implementation/CompressionCodeSamples` (35 files). See [Philips_Compression_Code_Samples_Claude.md](Philips_Compression_Code_Samples_Claude.md). `ImageFormatCodeSamples` is reviewed in [Image_Format_Code_Samples_Claude.md](Image_Format_Code_Samples_Claude.md); `ViewingCodeSamples` (CLAHE, sharpening) is reviewed in [Viewing_Code_Samples_Claude.md](Viewing_Code_Samples_Claude.md). | Bit-level definition of Hulsken, quantization and deadzone | Spec section 8; codeblock note |
| ~~Philips **Python reference**~~ **Now available and reviewed** (`ImageFormatCodeSamples/Python`, 4 files): see [Image_Format_Code_Samples_Claude.md](Image_Format_Code_Samples_Claude.md). It handles the v1 layout only. | Header parsing and codeblock indexing | Spec section 8; codeblock note |
| Hulsken "**Pathology iSyntax Compression**" document | The spec says to use it for reconstruction | Spec section 8 |
| **Patent text** (readable) | See 2.5 | n/a |
| A real **`.isyntax` file** (or just some codeblocks) | Needed to test the bitplane and mask predictions for quantizer and deadzone | Example_FIC_Files_Claude.md section 8 |
| Values for CLAHE/sharpening header attributes (0x1014–0x1018) | Pixel equivalence with the Philips viewer | Section 2.2 |

## 7. Suggested reading order

1. Format spec (md is enough for text; use the PDF for figures and equations).
2. Webinar, then Hulsken 2016, for the algorithm intent.
3. `quantanddeadzone.m` and `hulskendecompress`, once you add them. They settle goals 3 and 4.
4. Viewing document, if pixel equivalence with the Philips viewer matters.
5. The AI notes last, as a source of hypotheses: start with the "codeblock distribution" note's chunk and index sections, and the two YCoCg and 5/3 answers. Skip the AWS cost and file-size sections unless you need them.
