# Philips viewing code samples (CLAHE and sharpening): review

Source folder: `/mnt/d/Whole_Slide/iSyntax_Reference_Implementation/ViewingCodeSamples/` (8 `.m` files, 539 lines, Copyright 2020 Koninklijke Philips N.V.). This is the MATLAB/Octave reference for step 3 of the iSyntax pipeline, display post-processing. It is the code behind the viewing PDF (`4522 207 43951 … iSyntax image viewing.pdf`, summarised in [Third_Party_Documentation.md](Third_Party_Documentation.md) §2.2).

Companions: [Philips_Compression_Code_Samples_Claude.md](Philips_Compression_Code_Samples_Claude.md) (the codec), [Image_Format_Code_Samples_Claude.md](Image_Format_Code_Samples_Claude.md) (reading the file), [CLAUDE.md](CLAUDE.md).

**Evidence labels:** **[M]** = read in these `.m` files; **[doc]** = the viewing PDF; **[code]** = this repo's C; **[derived]** = computed by hand or by a short Python check; **[inferred]** = my deduction, untested.

## 0. Scope and limits

- **Nothing here was executed as MATLAB or Octave**, because neither is installed. Constants and formulas were checked with a short Python script (section 8). The behavioural claims come from reading the code.
- **This stage is not decoding.** It is image enhancement applied to the already-decoded RGB pixels. It has no connection to the Hulsken codec, the wavelet transform, or the quantizer and deadzone. It matters only if you want output that matches the Philips viewer.
- **libisyntax does none of this [code]**: it returns the decoded "source image". Its pixels will therefore differ from what the Philips viewer shows.

## 1. Inventory and call graph

| File | Lines | Role |
|---|---|---|
| `imageProcessingPipeline.m` | 98 | Entry point: read a PNG/RGB patch, run CLAHE then sharpening, write `<name>_IMAGE_PROCESSING_APPLIED.png` |
| `applyCLAHE.m` | 30 | RGB → YCoCg → CLAHE on Y → RGB |
| `calculateHistogram.m` | 87 | Per-region luma histograms, clipped and turned into equalisation maps |
| `clipHistogram.m` | 76 | Clip a histogram and redistribute the excess |
| `applyHistogramMap.m` | 155 | Bilinear interpolation of the four neighbouring region maps |
| `sharpening.m` | 41 | 5×5 integer unsharp mask on RGB |
| `rgb2ycocg.m`, `ycocg2rgb.m` | 26, 26 | YCoCg-R, using `bitshift(·,-1)` |

```
imageProcessingPipeline
  -> applyCLAHE:  rgb2ycocg -> calculateHistogram (-> clipHistogram) -> applyHistogramMap -> ycocg2rgb -> uint8
  -> [optional crop of the CLAHE/sharpening borders]
  -> sharpening(level)
  -> PNG
```

## 2. Pipeline and parameters **[M]**

- **Order:** CLAHE first, then crop, then sharpening. The header comment names it the "image processing pipeline of Pixel Engine", and it is meant to reproduce the output of Philips' `ISyntaxPatchExtrator.exe` (its `--ranges` option corresponds to the optional `xstart, xend, ystart, yend` arguments).
- **`level`** (default 0) is the pyramid level of the input: 0 is the highest resolution. It only affects the sharpening gain (section 5). In this repo's terms it is the scale index.
- **All the numbers are hard-coded in the code**:

  | Parameter | Value | Where |
  |---|---|---|
  | Contextual region size | **40** px (square) | `applyCLAHE.m`, `calculateHistogram.m` default |
  | Histogram bins | **256** | `calculateHistogram.m` |
  | Pixel range | 0–255 | `calculateHistogram.m`, `applyHistogramMap.m` |
  | Clip limit factor | **1.2** | `clipHistogram.m` default |
  | Sharpening gain | **2.0** at level 0, halved per level | `sharpening.m` |
  | Sharpening kernel | 5×5, sums to 128 | `sharpening.m` |
  | Border ignored by CLAHE | 40 px (`invalidTop`, `invalidLeft`) | `applyCLAHE.m` |
  | Border dropped by sharpening | 2 px | `sharpening.m` |

- **Header attributes.** The data model has attributes for exactly these settings: `DP_SHARPNESS_GAIN_RGB24` (0x1015), `DP_CLAHE_CLIP_LIMIT_Y16` (0x1016), `DP_CLAHE_NR_BINS_Y16` (0x1017), `DP_CLAHE_CONTEXT_DIMENSION_Y16` (0x1018). They are empty stubs in `isyntax.c:625-628` **[code]**, and **they do not appear in the spec or in any sample FIC**. The values above are what the reference code uses. That the header attributes would carry the same values (gain 2.0, clip 1.2, 256 bins, context 40) is **[inferred]**, from the names and the numbers.

## 3. Colour handling

- **`rgb2ycocg.m` and `ycocg2rgb.m`** use `int16` and `bitshift(x, -1)`, with the comment "as seen in H.Malvar and G. Sullivan 2003". `ycocg2rgb` asserts that its input class is `int16` and that it has three channels. These are the **viewing** versions; the compression folder has different files that use `fix()` (Philips_Compression_Code_Samples_Claude.md §6). How `bitshift` treats negative `int16` in Octave or MATLAB is unverified.
- **CLAHE changes only the luma Y**, and it leaves Co and Cg as they were. Because the YCoCg-R inverse is `tmp = Y − (Cg>>1); G = Cg + tmp; B = tmp − (Co>>1); R = B + Co`, a new luma `Y'` adds **the same amount `Y' − Y`** to R, G and B **[derived]**, before the final saturation to `uint8`. So the enhancement is a luminance shift with chroma preserved, as the PDF says ("hue and saturation unchanged").
- **The cast back to `uint8` saturates** (MATLAB and Octave integer casts clamp to 0..255). That is where out-of-range results become clamped.

## 4. CLAHE in detail

### 4.1 Histograms (`calculateHistogram.m`)

1. Convert to YCoCg and take the **Y channel only**.
2. **LUT:** `LUT(i+1) = floor((i + 1) · 256 / 256) = i + 1`, so the "bin" of a luma value `v` is simply `v + 1` (1-based). With 256 bins and range 0–255 the LUT is the identity.
3. **Grid:** the image is cut into 40×40 contextual regions: `rows = floor((H + 39)/40)`, likewise columns. Last-row and last-column regions may be smaller.
4. Each region's histogram counts its luma values.
5. **Full regions** (40×40): the histogram is **clipped** (4.2) and then turned into an **equalisation map**: the running sum of the clipped counts is scaled by `255 / (rows·cols)` and floored, capped at 255:
   `map[b] = min(255, floor(cumulative_count(b) · 255 / 1600))`. Since clipping keeps the total at 1600, the last entry is 255.
6. **Partial regions** (image edges) are **copied from a full neighbour** rather than computed: a partial region in the last column takes the histogram map of the region to its left, and one in the last row takes the one above. (A third branch, which would only clip a partial region without equalising it, cannot be reached, since a partial region is always in the last row or column **[inferred]**.)

### 4.2 Clipping (`clipHistogram.m`)

- `clipLimit = floor(1.2 · floor(numEntries / numBins))`. For a 40×40 region, `numEntries = 1600` and `numBins = 256`, so `floor(1600/256) = 6` and **`clipLimit = floor(7.2) = 7` counts per bin** **[derived]**. For partial regions the code still passes the **full-region** count (1600), "otherwise there could be an infinite loop".
- Counts above the limit are cut to it. The total excess is spread evenly (`binIncr = floor(excess / bins)` to each bin), bins near the limit are topped up to it, and the remaining excess is handed out one count at a time with a changing stride until none is left, or until no bin can take more (the `excessPixelsStored` check).
- A comment notes that rounding in `binIncr` "is particularly visible for cliplimit == 1".

### 4.3 Applying the maps (`applyHistogramMap.m`)

- **Region maps are interpolated bilinearly** between the four nearest region centres. The grid of interpolation cells is shifted by half a region (20 px), so a pixel is blended between the maps of the four regions whose centres surround it.
- **The first 40 rows and columns are not modified** (`invalidTop`, `invalidLeft` = 40). The first interpolation cell starts at row 20, but processing begins at row 40, which is exactly halfway between the centres of region 0 (row 20) and region 1 (row 60). The code computes this as `histRegRowOffset = 20` **[derived]**.
- **Edge cells** use fewer maps: the top and left half-cells use a single map in that direction, and the bottom and right edges clamp to the last map. The corners therefore use one map.
- **Weights** are `(h − ridx)·((w − cidx)·LU + cidx·RU) + ridx·((w − cidx)·LL + cidx·RL)`, divided by `h·w`, rounded with `floor(x + 0.5)`. The pixel's own luma selects the entry in each of the four maps.
- A Philips modification clamps the luma used as a map index to the valid range (`pixelMin`, `pixelMax`).
- **Guard:** if the image is no larger than a region (`≤ 40`) in either dimension, the function asserts and stops. The processing therefore needs a patch of more than 40 px on each side.

## 5. Sharpening (`sharpening.m`)

The code is a fixed-point unsharp mask applied to **each of R, G and B separately** (not to luma):

```
sum  = imfilter(wsi, kernel)            % kernel: [0 0 1 0 0; 0 5 15 5 0; 1 15 44 15 1; 0 5 15 5 0; 0 0 1 0 0]
gain = 2.0 * 2^5 / 2^level
temp = floor( (gain*(wsi*128 - sum) + wsi*4096 + 2048) / 4096 )
ret  = uint8(temp), with the outer 2 pixels cropped on every side
```

- **The kernel sums to 128** **[derived]**, so `128·x − sum` is 128 times the high-pass signal `x − LP(x)`, where `LP` is the 5×5 low-pass filter. The expression reduces to
  `out = round( x + (gain/32) · (x − LP(x)) )`, an **unsharp mask with gain 2.0 · 2^(−level)**:

  | Level | 0 | 1 | 2 | 3 |
  |---|---|---|---|---|
  | Effective gain | 2.0 | 1.0 | 0.5 | 0.25 |

  **[derived]**, for the code's constants. The gain is halved per level "to prevent over-enhancement", as the PDF says, because the wavelet downscale already sharpens.
- `imfilter` pads with zeros by default, so the outer two pixels are wrong, and the code discards them (`b = 2`). The saturating `uint8` cast clamps the result.
- **It runs after CLAHE**, on the contrast-enhanced RGB.

## 6. The patch borders (`imageProcessingPipeline.m`)

To reproduce the Philips viewer exactly you must feed the pipeline a padded patch and crop afterwards. The crop function (`cropImage`) shows how the pipeline pads:

- The original requested range `[xstart, xend] × [ystart, yend]` is first widened by `claheBorder + sharpeningBorder = 40 + 2 = 42` pixels on every side, and then **rounded outward to a multiple of 40** (start rounded down, end rounded up), because CLAHE regions sit on a 40 px grid anchored to the patch origin.
- After CLAHE and sharpening, the result is cropped back to the requested range.

Consequence: the region grid is anchored at the **top-left of the image that is passed in**. Because the padded patch always starts at a multiple of 40 in the slide's level coordinates, the regions then coincide with a fixed 40 px grid over the whole level, and the result for a given pixel does not depend on which patch it was extracted in **[inferred from `cropImage`]**. A patch that starts anywhere else would use different region boundaries and give slightly different pixels. This is why a pixel-exact imitation of the Philips viewer needs the same padded, 40-aligned extraction, not just the same filter.

## 7. Differences between the `.m` files and the viewing PDF

I compared the PDF's listings (text extracted earlier) with the `.m` files. They are not identical:

| Item | PDF listing | `.m` file | Note |
|---|---|---|---|
| `endRowRegion`, `endColRegion` in `applyHistogramMap` | `idivide(invalidTop + imHeight - 1 + halfCRegionSize, contextualRegionSize, floor)` | the same expression **minus 1** | Could change which cells are processed at the far bottom and right edge. Those pixels are in the zone that the 42 px border and the crop remove, so I expect no visible effect in a correctly padded patch **[inferred]**, but I did not run either version. |
| Partial-region branch in `calculateHistogram` | refers to `colRegion`, `numColRegions_`, `rowRegion`, `numRowRegions_` (not defined anywhere) | uses `colContext`, `nrContextCols`, `rowContext`, `nrContextRows` | The PDF listing has stale variable names. The `.m` is the working version. |
| Excess redistribution in `clipHistogram` | `histogram(idx) += binIncr` | `histogram(idx) += floor(binIncr)` | `binIncr` is already an integer (`idivide`), so the same result. |
| Small-image guard in `applyHistogramMap` | absent | assert if image ≤ 40 px | Added in the code. |
| Everything else I compared (kernel, gain, clip factor, region size, bins, YCoCg-R formulas, interpolation weights) | | | Same |

Philips states that following the PDF gives output "pixel equivalent" to the IMS viewer. Given the table above, **the `.m` files, not the PDF listing, should be treated as the reference**, and the two should be reconciled before anything is claimed to be pixel-exact.

## 8. Checks I ran

A short Python script **[derived]**:

- `floor(1.2 · floor(1600/256)) = 7` (the clip limit per bin).
- The kernel sums to 128.
- Effective sharpening gain `gain · 128 / 4096` is 2.0, 1.0, 0.5, 0.25 for levels 0 to 3.
- The first interpolation cell: `startRowRegion = 1`, `regionStartRow = 20`, `histRegRowOffset = 20`.

## 9. Code quality and portability notes

- **The code is written in Octave syntax**: `+=`, `size(x)(1)`, `printf`, `endfunction`. MATLAB rejects `size(yCoCgImage)(1)` (chained indexing of a function result), so `calculateHistogram.m` and `applyHistogramMap.m` would need edits for MATLAB **[inferred; not run]**, even though the PDF says the code was verified with Matlab 9.8. `rgb2ycocg.m` and `ycocg2rgb.m` carry a "revised to satisfy MatLab R2020a" comment, so the files were partly adapted.
- **Needs an image toolbox:** `imfilter` (Image Processing Toolbox in MATLAB, `image` package in Octave) and `imread`/`imwrite`.
- **Leftover path handling:** `addpath('./CLAHE')` and `addpath('..')` refer to directories that do not exist in this flat folder, and the output name is built with a Windows `\` separator.
- **Speed:** the histogram and interpolation loops are pixel-by-pixel scalar code. This is a reference for correctness, not a fast implementation (the PDF says the same).
- **Single-patch only:** it processes one image at a time and holds it in memory.

## 10. What this means for your goals and for libisyntax

| Goal | Relevance |
|---|---|
| 1. File format and FIC header | None directly. It explains the three unused header attributes 0x1015–0x1018 (section 2). |
| 2. 5/3 wavelet | None |
| 3. Hulsken compression | None |
| 4. Quantizer and deadzone | None |
| Pixel match with the Philips viewer | **High.** A match needs CLAHE (40×40, 256 bins, clip 7 per bin, bilinear interpolation) on luma, then level-dependent sharpening, on patches padded to a multiple of 40 plus a 42 px frame. libisyntax provides none of it. |
| Wavelet-domain transcoding (your earlier notes) | The post-processing is applied at display time and is not part of the stored data, so a faithful DICOM transcode should **not** bake it in. The sharpening gain is tied to the pyramid level, which is worth remembering if DICOM pyramid levels are viewed through a Philips-style pipeline. |

## 11. Open questions

1. Do the header attributes 0x1015–0x1018 carry the values hard-coded here (2.0, 1.2, 256, 40) in real files? None of the sample headers has them. A newer or different scanner file might.
2. Which `applyHistogramMap` is the canonical one: the PDF (no `−1`) or the `.m` file (`−1`)? Running both on a test patch would settle it, and would show whether the interior and cropped result differs at all. This needs Octave.
3. Exact semantics of `bitshift(·, −1)` on negative `int16` in Octave and MATLAB, which decides how the viewing colour transform rounds (Philips_Compression_Code_Samples_Claude.md §6).
4. Whether any production viewer applies CLAHE and sharpening **per tile** (each tile padded and extracted on its own 40 px grid) or per larger patch. The reference only shows the patch function. This affects visible seams between tiles.
