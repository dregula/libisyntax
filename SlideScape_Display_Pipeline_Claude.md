# Do SlideScape and QuPath apply CLAHE or sharpening? (source reviews, 2026-11-27)

Written for the project owner. Source: `https://github.com/LUMCPathAI/slidescape`, commit `f6e2461` (2026-11-26), cloned to a temporary directory and searched; owner's visual observation: SlideScape looks sharper than OmniHisto showing the transcoded DICOM. Tags: **[verified in code]**, **[measured]**, **[hypothesis]**.

## Findings

1. **No CLAHE, no sharpening, no contrast equalisation anywhere in the viewer [verified in code].** The only occurrences of "CLAHE" are three empty `{}` cases in `src/isyntax/isyntax.c:627-629` (`DP_CLAHE_CLIP_LIMIT_Y16` 0x1016, `DP_CLAHE_NR_BINS_Y16` 0x1017, `DP_CLAHE_CONTEXT_DIMENSION_Y16` 0x1018): the header values are parsed and discarded, as in libisyntax. A search for sharpen, unsharp, equaliz, convolution, kernel and gamma in `src/core`, `src/isyntax` and `src/platform` finds nothing relevant.
2. **The shaders do no filtering [verified in code].** `shaders/basic.frag` samples the tile texture and applies an optional black/white level stretch `(c - black) / (white - black)`; this is used only when `use_image_adjustments` is on, otherwise the uniforms are 0 and 1 (`src/core/renderer_opengl.c:347-353`). `finalblit.frag` only cross-fades two textures.
3. **Texture filtering [verified in code].** `src/core/renderer_opengl.h:84-85`: magnification `GL_NEAREST`, minification `GL_LINEAR_MIPMAP_LINEAR`. When zoomed in beyond 1:1, SlideScape therefore shows crisp, blocky pixels without interpolation; most viewers interpolate (bilinear or bicubic) when magnifying, which looks smoother. A software renderer (`use_fast_rendering`) forces nearest for both.
4. **Colour conversion in this fork [verified in code].** `ycocg_to_rgb_8bit` (`src/isyntax/isyntax.c:1672`) uses truncating `/2` (same as libisyntax's scalar path) and **adds a lower clamp at 0** that libisyntax's scalar path lacks (relevant for the few negative pixels found on files 12 and 4).
5. **The header does carry CLAHE parameters [measured].** `1.isyntax` contains `DP_CLAHE_CLIP_LIMIT_Y16 = 1.2`, `DP_CLAHE_NR_BINS_Y16 = 256`, `DP_CLAHE_CONTEXT_DIMENSION_Y16 = 40` (group 0x301D, elements 0x1016 to 0x1018); the six sample `.fic` files do not. The clip limit 1.2 differs from the constant 7 hard-coded in Philips' viewing reference (`Viewing_Code_Samples_Claude.md`); **[open]** whether the header value is a different normalisation or the reference code uses a different default. Neither SlideScape nor libisyntax applies it, so any CLAHE in Philips' own viewer is not in the images being compared.

## Is the transcoded DICOM less sharp? [measured]

No. For three tissue frames of layer 0 (`1.isyntax`, DICOM frames decoded with OpenJPEG) the luma detail equals the iSyntax Y plane: Laplacian variance of the DICOM luma 22.5 / 29.4 / 15.0, of the stored Y plane 22.5 / 29.4 / 15.0 (identical), and the RCT luma of the DICOM RGB differs from the SlideScape-style RGB by at most 1 level. The SlideScape-style decode (truncating YCoCg to RGB, then luma) measures 15 to 30 % higher (26.7 / 34.2 / 19.5), which is the +-1 rounding noise of the RGB conversion, not extra detail. The wavelet-domain colour conversion changes only chroma (mean absolute RGB difference 0.7 / 2.6 / 2.3 levels in R / G / B).

## Likely causes of the perceived difference [hypothesis]

- The viewers' **magnification filters** (SlideScape nearest, OmniHisto probably smoothing), and which layer each shows at a given zoom: an off-by-one layer choice looks twice as blurry.
- OmniHisto may apply colour management with the embedded sRGB profile or its own display processing.
- A fair test: compare both viewers at exactly 100 % on layer 0 (one image pixel per screen pixel), where no interpolation is involved. If SlideScape is still sharper there, look at colour management and display sharpening in OmniHisto.

## QuPath v0.7.0 (tag `v0.7.0`, commit `04ccfa4`, 2026-02-25) **[verified in code]**

- A search of the whole source tree finds **no CLAHE, sharpen, unsharp or equalisation code** (no hits for `clahe`, `sharpen`, `unsharp`, `equaliz`). QuPath's processing code has Gaussian blur and similar filters (`ImageOps.Filters`) for analysis, not for display; the viewer shows the pixels with an optional per-channel brightness/contrast (and gamma) setting.
- **Magnification interpolation is nearest-neighbour by default**: `PathPrefs.viewerInterpolateBilinear` defaults to `false` (`PathPrefs.java:1411`), and `QuPathViewer.java:2108-2110` selects bilinear only when that preference is on. QuPath should therefore look crisp like SlideScape, which makes it a good second viewer for the sharpness question.
- DICOM whole slide images are read through **Bio-Formats 8.4.0** (`gradle/libs.versions.toml`; the changelog lists DICOM WSI support since the Bio-Formats 6.11 era). Whether Bio-Formats assembles our six layer files into one pyramid when `1_L0.dcm` is opened is **[open]** (the owner will report).
- Launching (owner-granted, outside the sandbox): `QuPath-0.7.0.exe C:\Users\Public\iSyntax_review\dicom_1\1_L0.dcm`; the program started (JVM process present), whether the image opened was not checked by me.
