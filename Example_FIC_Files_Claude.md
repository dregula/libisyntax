# Example FIC headers: analysis

Companion to [CLAUDE.md](CLAUDE.md). It covers the six sample iSyntax XML headers ("FIC") in `../example_FIC_files/`: `1.fic 2.fic 4.fic 9.fic 10.fic PAS.fic`.

**Evidence labels** used throughout:

- **[FIC]**: read directly from the sample headers.
- **[code]**: read from this repo's source (`file:line`).
- **[spec]**: stated in `../third_party_documentation/2020_04_24_30pp_Pathology_iSyntax_image_format_free-online-conversion.md`.
- **[inferred]**: my deduction. Treat it as a hypothesis until tested against codeblock data.

## 0. What these files are, and what they are not

- Each file is the **XML header only**: one `DPUfsImport` root, no `\r\n\x04` EOT marker, no seektable, no codeblocks. They cannot be decoded to pixels. Nothing here tests the Hulsken decoder or the IDWT; only the header and the cluster table can be analysed.
- Base64 payloads were deleted by hand **[per user]**: the LABELIMAGE and MACROIMAGE JPEGs (`PIM_DP_IMAGE_DATA` 0x1005 is empty) and the ICC colour tables (`DICOM_ICCPROFILE` 0x0028,0x2000 is empty). The barcode (0x1002) is intact.
- **The cluster header table (0x201f) is intact**, and it holds the compressed size of every codeblock. This is the only quantitative evidence about the bitplane payload in these files (section 7).
- Absolute file offsets in the table refer to the **original** file. For example, `1.fic` has its first cluster at byte 4,390,072 while the stripped file is 529,582 bytes. Offsets cannot be used to seek in these files.
- All six files parse as well-formed XML (checked with Python's ElementTree).

## 1. Inventory

All values **[FIC]** except the derived columns, which I computed from the ranges (rule in section 4).

| | 1 | 2 | 4 | 9 | 10 | PAS |
|---|---|---|---|---|---|---|
| Acquired | 2019-06-12 23:16:59 | 23:12:09 | 23:15:40 | 23:28:14 | 23:26:24 | **2020-02-03 18:42:34** |
| Scanner rack/slot | 1/9 | 1/6 | 1/8 | 1/16 | 1/15 | **2/4** |
| Serial | FMT0132 | FMT0132 | FMT0132 | FMT0132 | FMT0132 | **FMT0425** |
| Last calibration | 20190612 230908 | same | same | same | same | **20200203 175508** |
| Barcode (base64 → text) | 7 spaces | 7 spaces | 7 spaces | 7 spaces | 7 spaces | `39V79D` |
| Padded size (px, level 0) | 25600×21504 | 206848×76800 | 165888×81920 | 89088×73728 | 99328×74752 | 103424×78848 |
| Levels (scales 0..L−1) | 6 | 9 | 9 | 8 | 8 | 8 |
| Level-0 padding / side | 189 | 1533 | 1533 | 765 | 765 | 765 |
| Valid size (px) | 25222×21126 | 203782×73734 | 162822×78854 | 87558×72198 | 97798×73222 | 101894×77318 |
| Valid size (mm at 0.25 µm/px) | 6.3×5.3 | 50.9×18.4 | 40.7×19.7 | 21.9×18.0 | 24.4×18.3 | 25.5×19.3 |
| `UFS_IMAGE_NUMBER_OF_BLOCKS` | 33,867 | 304,251 | 228,762 | 391,077 | 442,422 | 498,357 |
| Clusters in 0x201f table | 537 | 4,829 | 3,631 | 6,214 | 7,029 | 7,919 |
| Cluster header templates | 2 | 3 | 3 | 3 | 3 | 3 |
| Block header templates | 21 | 30 | 30 | 27 | 27 | 27 |
| Valid-data envelopes | 1 | 5 | 2 | 1 | 1 | 1 |
| Sum of codeblock bytes | 84.3 MB | 801 MB | 600 MB | 1,641 MB | 1,370 MB | 1,528 MB |
| Valid RGB bytes ÷ codeblock bytes | 19.0 | 56.3 | 64.2 | 11.6 | 15.7 | 15.5 |

Notes:

- **Constant across all six** [FIC]:
  - `PIM_DP_UFS_INTERFACE_VERSION = 100.4`, so these are **iSyntax v2** (data model ≥ 100).
  - Software versions `"1.8.6614" "20180906_R51"`.
  - `DICOM_DERIVATION_DESCRIPTION = PHILIPS UFS V1.8.6614 | Quality=2 | DWT=1 | Compressor=16`.
  - Lossy compression ratio 15 and method `PHILIPS_DP_1_0`.
  - Block compression method 16, in every one of the block templates.
  - Quantizer and deadzone tables (section 6).
  - 0.25 µm per pixel, with `Y Co Cg` colour components.
  - Block templates declare 3 samples per pixel, 16 bits allocated, 8 stored, high bit 7.
- **No v1 sample is present.** Every file is data model 100.4. The v1 path in the code (`UFS_IMAGE_BLOCK_HEADER_TABLE` 0x2014) is never exercised here. The attribute is absent from all six files. Philips' own Python reference reads only that v1 table, so it cannot process these headers past the metadata step; I ran its metadata reader on `1.fic` and `PAS.fic` and it reproduced the valid sizes in the table above exactly (see [Image_Format_Code_Samples_Claude.md](Image_Format_Code_Samples_Claude.md) §3.1).
- **Compressor v1 inside a v2 data model.** Method 16 means Hulsken compressor v1 (`isyntax.c:729`) even though the container is v2. No sample uses method 19 (compressor v2).
- **The declared ratio of 15 is nominal.** Tissue-dense slides (9, 10, PAS) land at 11.6–15.7 on a valid-area RGB basis. Slides with large empty areas (2, 4) reach 56–64, because blank regions have no stored data.
- **Padded dimensions are multiples of 256** in every file. A 256 px tile is the output of one IDWT block, i.e. 2 × 128.

## 2. Document structure [FIC]

```
DPUfsImport
├─ DICOM_* scalars (0x0008/0x0018 groups): datetime, manufacturer PHILIPS, model "UFS Scanner",
│     serial, software versions, calibration date/time
├─ PIIM_DP_SCANNER_RACK_NUMBER 0x1007, _SLOT_NUMBER 0x1008, _CALIBRATION_STATUS 0x100a ("OK")
├─ PIM_DP_UFS_INTERFACE_VERSION 0x1001 = 100.4,  PIM_DP_UFS_BARCODE 0x1002 (base64)
├─ PIM_DP_SCANNED_IMAGES 0x1003  →  3 × DPScannedImage, always in this order:
│    1. WSI
│    2. LABELIMAGE
│    3. MACROIMAGE
└─ PIM_DP_SCANNER_RACK_PRIORITY 0x1010 = 0
```

**WSI `DPScannedImage`**, in document order:

| Element | Content |
|---|---|
| 0x0008,0x2111 / 0x0028,0x2110 / 0x2112 / 0x2114 | Derivation description, lossy flag `01`, ratio 15, method `PHILIPS_DP_1_0` |
| 0x1004 | `WSI` |
| 0x1013 `DP_COLOR_MANAGEMENT` | One `DPColorManagement` with empty ICC (stripped) |
| 0x1019 `DP_WAVELET_QUANTIZER_SETTINGS_PER_COLOR` | 2 × PerColor → PerLevel[] → {0x101b, 0x101c}. See section 6. |
| 0x2000 `UFS_IMAGE_GENERAL_HEADERS` | One `UFSImageGeneralHeader` (see below) |
| 0x2009 `UFS_IMAGE_BLOCK_HEADER_TEMPLATES` | 21 / 27 / 30 `UFSImageBlockHeaderTemplate` |
| 0x2016 `UFS_IMAGE_CLUSTER_HEADER_TEMPLATES` | 2 or 3 `UFSImageClusterHeaderTemplate`, each with a list of `UFSImageBlockHeader` (63 or 66 or 18) |
| 0x201f `UFS_IMAGE_CLUSTER_HEADER_TABLE` | base64 binary table, 0.3–4.7 MB decoded |
| 0x2023 `UFS_IMAGE_VALID_DATA_ENVELOPES` | `UFSImageValidDataEnvelope` → 4 × `UFSImageOppExtremeVertex` (0x2025), plus 0x2026 dims `0 1` |

`UFSImageGeneralHeader` contains:

- 0x2001 number of blocks.
- 0x2002 `DIMENSIONS_OVER_BLOCK = 0 1 2 3 4`.
- 0x2017 `DIMENSIONS_OVER_CLUSTER = 0 1 3`.
- 0x2003 five `UFSImageDimension` objects: `x` and `y` (spatial, MicroMeter, scale factor 0.25, origin), `component` (`"Y" "Co" "Cg"`), `scale`, and `waveletcoef` (`"HL" "LH" "HH" "LL"`).
- 0x200a five dimension ranges.

**LABELIMAGE and MACROIMAGE** carry only: lossy ratio 26, method `ISO_10918_1` (JPEG), image type, empty image data, and a two-dimension general header at **32 µm/px**. The macro image also has an (empty) colour management object.

| | Label | Macro |
|---|---|---|
| Size in file 1 | 783×733 (range 0 1 782 / 732) | 1814×819 |
| Origin (µm) | x 58016, y 0 | 0, 0 |

Label and macro sizes and origins vary by file. I only compared them for file 1.

## 3. Attribute facts the code does not use

These are present in the files but `isyntax.c` ignores them or handles them as `{} break;`.

- **0x1019 / 0x101a / 0x101b / 0x101c**: quantizer settings, empty stubs at `isyntax.c:629-632` **[code]**. This is the core of section 6.
- **0x2027 `UFS_IMAGE_DIMENSION_ORIGIN`**: x/y origin in µm, for example 19428.8 / 5153.5 in file 1. It is an empty stub (`isyntax.c:1030`), so the slide-stage position is not exposed by the library.
- **0x2002 / 0x2017 / 0x200c**: dimension-ordering lists. The cluster one (0x2021 `IN_CLUSTER`) **is** used to order raw coordinates (`isyntax.c:1008`). The others are not needed because the order is fixed.
- **0x1007, 0x1008, 0x100a, 0x1010**: scanner rack, slot, status, priority.
- **Label/macro compression ratio (26) and method (`ISO_10918_1`)** are read, but the library decodes those images through the JPEG payload anyway.

## 4. Coordinates, ranges, and padding

Dimension ranges are `start step end` (`UFS_IMAGE_DIMENSION_RANGE`, 0x200b). The code computes `numsteps = ((end + step) − start) / step` (`isyntax.c:665`). For the WSI general header, the five ranges are:

| Dim | File 1 | Meaning |
|---|---|---|
| x | `0 1 25599` | 25600 padded px |
| y | `0 1 21503` | 21504 padded px |
| component | `0 1 2` | Y, Co, Cg |
| scale | `0 1 5` | **6 levels** (files 2/4: `0 1 8` = 9; files 9/10/PAS: `0 1 7` = 8) |
| waveletcoef | `0 1 3` | HL, LH, HH, LL |

**Padding rule, verified on all six files.** The code derives `level0_padding = (3 << level_count) − 3` (`isyntax.c:692`). That gives 189 for 6 levels, 765 for 8 and 1533 for 9. These are exactly the minimum coordinates of the valid-data envelopes in every file:

- File 1 envelope: x 189..25411. 25600 − 189 = 25411, so the high edge is symmetric.
- File 4 first envelope starts at `1533 1533`.
- Files 9, 10 and PAS start at `765 765`.

So the code's padding formula **[code]** agrees with the file **[FIC]**. The valid size in section 1 is padded size minus 2 × padding.

**Valid-data envelopes** are rectangles given as 4 vertices (all corners, called "opposite extreme vertices"). Files 1, 9, 10 and PAS have one envelope covering essentially the whole padded area. File 2 has five disjoint tissue regions and file 4 has two. Envelope bounds are in padded level-0 pixel coordinates. The code parses and stores them (`valid_data_envelopes[16]`, `isyntax.c:1019`, `:1411`), but a grep finds no other use, so they are not used to skip tiles.

**Block coordinates in templates use level-0 pixel units with a step that doubles per scale.**

- Block template at scale `s` has x/y range `0 2^(s+1) (128·2^(s+1)−2)`, i.e. 128 steps. Scale 0: `0 2 254`. Scale 5: `0 64 8128`.
- One 128×128 coefficient block at scale `s` therefore covers 256·2^s level-0 pixels per side. At scale 0 that is the 256 px tile.

## 5. Templates

### 5.1 Block header templates (0x2009)

One template per (colour, scale, band group). All of them have `BITS_ALLOCATED 16`, `BITS_STORED 8`, `HIGH_BIT 7`, `SAMPLES_PER_PIXEL 3` and compression method 16.

- Per colour `c ∈ {0,1,2}` and scale `s ∈ [0, L−1]`, there is one template with waveletcoef range `0 1 2`. That is the 3-band group (HL, LH, HH) in a single codeblock.
- At the top scale there is one more template per colour with waveletcoef range `3 1 3`, i.e. **LL**. Its `DIMENSIONS_IN_BLOCK` lists only x and y, not all five.
- Count: 3 × (L + 1). That is 21 for L = 6, 27 for L = 8 and 30 for L = 9 **[FIC]**, matching the table above.
- Template ids run colour-major. In file 1, ids 0–6 are Y (scales 0–5 and LL), 7–13 are Co, 14–20 are Cg.

**A discrepancy in the code**, worth fixing in comments when convenient:

- `isyntax.h:178` comments `waveletcoeff` as "either 1 for LL, or 3 for LH+HL+HH".
- The code at `isyntax.c:674` sets `waveletcoeff = (range.start == 0) ? 1 : 3`. In the FIC, start 0 is the HL/LH/HH group and start 3 is LL.
- The cluster path agrees with the FIC, not the comment: `codeblock->coefficient = (waveletcoeff == 3) ? 0 : 1` (`isyntax.c:970`), and `coefficient == 0` marks LL (`isyntax.c:971`, `isyntax_reader.c:110`).
- So the real convention is **`waveletcoeff` 3 = LL, 1 = the 3-band H group**, and the comment is inverted.

### 5.2 Cluster header templates (0x2016)

A **cluster** is a group of codeblocks stored contiguously: 3 consecutive scales × 3 colours. Cluster templates list each block's relative coordinate and its block-header-template id, in the order the blocks appear in the data.

| Cluster template | Scales | Blocks per colour | Headers | Contains LL | Present in |
|---|---|---|---|---|---|
| 0 | 0, 1, 2 | 16 + 4 + 1 | 63 | no | all files |
| 1 | 3, 4, 5 | 16 + 4 + 1 | 63 (66 in file 1) | file 1 only | all files |
| 2 | 6, 7, 8 | 16 + 4 + 1 (9 levels: 66, with LL) or 4 + 1 (8 levels: 18, with LL) | 66 or 18 | yes | files 2, 4, 9, 10, PAS |

The 16/4/1 pattern follows from each scale halving the tile grid in x and y. The last cluster holds whatever scales remain, plus the LL at the top. In file 1, cluster template 1 is the top cluster, scales 3–5 plus LL, so 66 headers.

Cluster range `ranges` examples **[FIC]**:

- Template 0: `(1,2,1023)` in x and y. This spans 1024 level-0 px, which is 4×4 tiles at scale 0.
- Template 1: `(22,16,8198)`.
- Template 2 (files with ≥ 8 levels): `(190,128,65598)` for 9 levels, `(190,128,32830)` for 8.
- **The range starts are the first valid coefficient pixel of the cluster's base scale**: 1, 22 and 190 equal `(3 << s) − 2` for s = 0, 3 and 6, which is the code's `get_first_valid_coef_pixel` (`isyntax.c:3047`). Verified in all six files. The step is `2^(base_scale+1)`, so 2, 16 and 128.

`UFS_IMAGE_DIMENSIONS_IN_CLUSTER` (0x2021) is `0 1 3 2` for clusters that have no LL, and `0 1 3 4 2` for clusters that do. In the order used by the code's `dimension_order`, this is x, y, scale, [waveletcoef], colour. That matches the order in which `UFS_IMAGE_BLOCK_COORDINATE` values are listed **[FIC]**: for example `515 3 1 0` for a block with colour 1 at scale 1 (x 515, y 3).

Relative coordinates within a cluster: in cluster template 0, scale 0 blocks have offsets 0, 256, 512, 768. Scale 1 blocks sit at 3 and 515. The scale 2 block sits at 9. In template 1 the same pattern appears multiplied by 8 (24, 72). A pattern of `3·(2^r − 1)` for scale offset `r` from the base scale, times the base-scale unit, fits these offsets **[inferred]**; I only checked it on file 1's first blocks.

## 6. Quantizer and deadzone

### 6.1 What the header says [FIC]

Identical in all six files, with 2 `DPWaveletQuantizerSettingsPerColor` objects:

| PerColor entry | Level 0 (Q, DZ) | Level 1 (Q, DZ) | Level 2 (Q, DZ) |
|---|---|---|---|
| 1st | (1, 2) | (0, 0) | none |
| 2nd | (4, 16) | (2, 8) | (0, 0) |

### 6.2 What the spec says [spec]

- **Table 11, PerColor**: "Per color component a list of quantizer settings. First entry belongs to first color component etc… The typical order for colors is either Y-Co-Cg or R-G-B". Range note: "**Typical 2 objects: first for luminance, second for the other color components.** Only applicable to WSI."
- **Table 11, PerLevel**: "Per level a list of quantizer settings. First settings belong to level 0, second to level 1 etc."
- **Table 12, `DP_WAVELET_QUANTIZER`** (0x101b): "Wavelet coefficients are **rounded to 2^quantizer**", range 0..n.
- **Table 12, `DP_WAVELET_DEADZONE`** (0x101c): "Wavelet coefficients are **pulled to zero for absolute values greater or equal to dead zone**", range 0..n.

Reading these together:

- The first PerColor entry is **Y** (levels 0 and 1).
- The second PerColor entry applies to **Co and Cg** (levels 0, 1 and 2).
- This resolves the "two groups, three colours" question left open in `CLAUDE.md`.
- Levels beyond the listed ones have no setting. Presumably they are not quantized, and an explicit (0, 0) is the same as "none".
- **The deadzone wording is wrong, confirmed by Philips' own encoder.** `quantanddeadzone.m` zeroes coefficients with `abs(c) <= deadzone`, and truncates the rest toward zero to a multiple of `2^Q` (`fix(c/2^Q)*2^Q`). See [Philips_Compression_Code_Samples_Claude.md](Philips_Compression_Code_Samples_Claude.md) section 3. That file also confirms column 1 of the parameter matrices is Y and column 2 is Co and Cg, that `Quality=2` is the "Q2" set whose values appear in these headers, and that level 0 is the first (finest) DWT level.

### 6.3 How it connects to the code [code] + [M] (the Q/DZ semantics are now confirmed by the Philips encoder)

- `isyntax.c` never reads these values (stubs at lines 629–632), and the decode path has no dequantization.
- That fits "rounded to 2^quantizer" if the **rounded values keep their original scale**: the encoder zeroes the low `Q` magnitude bits rather than dividing. Then there is nothing for the decoder to multiply back. **[inferred]**
- In the Hulsken bitplane layout (`isyntax.c:2975`, v1): mask bit 0 is the sign plane and mask bits 1…15 are magnitude bits from LSB to MSB. Under the rounding reading, a level with quantizer `Q` should show:
  - **mask bits 1…Q cleared** (the `Q` lowest magnitude planes absent), for every H-band codeblock of that colour and level;
  - **deadzone zeros as zero runs** in the planes that remain.

**Both of these are predictions.** The stripped files hold no codeblocks, so they cannot be checked here. A decoded `.isyntax` file is needed. The check is cheap: read the 6 mask bytes of a Y/Co/Cg codeblock (v1: they are at the end of the serialized data, `isyntax.c:2889`).

### 6.4 Supporting evidence from codeblock sizes [FIC]

The cluster header table gives the stored size of every codeblock. The code treats a block of ≤ 8 bytes as empty (`isyntax.c:2577`), and the median size of the sparse groups below is exactly 8. Percentage of H-band codeblocks that are empty, by scale and colour:

| Scale | Setting (Y / chroma) | File | Y | Co | Cg |
|---|---|---|---|---|---|
| 0 | (1,2) / (4,16) | 1 | 5.5% | **81.7%** | **92.6%** |
| 0 | | 2 | 4.9% | 77.3% | 92.9% |
| 0 | | 4 | 4.2% | 83.1% | 91.9% |
| 0 | | 10 | 2.1% | 87.3% | 96.2% |
| 0 | | PAS | 2.8% | 72.0% | 90.9% |
| 1 | (0,0) / (2,8) | 1 | 4.3% | 8.6% | 33.1% |
| 1 | | 2 | 3.4% | 6.0% | 26.9% |
| 1 | | 4 | 2.6% | 7.2% | 34.6% |
| 1 | | 10 | 1.7% | 3.2% | 15.0% |
| 1 | | PAS | 2.3% | 8.0% | 39.7% |
| 2 | (none) / (0,0) | 1 | 0.0% | 0.0% | 0.0% |
| 2 | | 4, 10, PAS | 0.0% | 0.0% | 0.0% |

Mean H-block size at scale 0, file 1: **Y 1001 B, Co 16 B, Cg 11 B**. At scale 1: Y 14,786 B, Co 1,356 B, Cg 558 B. At scale 2: 17.6–20.4 KB for all three colours.

What this shows:

- **Chroma at scale 0 is almost entirely zeroed** (72–96% empty blocks), consistent with the strongest setting (Q = 4, DZ = 16).
- **Chroma at scale 1** is much less empty, consistent with the milder (2, 8).
- **At scale 2 nothing is empty**, consistent with (0, 0).
- **Y is barely affected by its (1, 2) setting** at scale 0: its empty rate (2–5%) is close to Y at scale 1 (2–4%, setting (0, 0)). Blank glass produces genuinely zero coefficients, so a small baseline is expected anyway.
- **Mean size at scale 0 drops by about two orders of magnitude from Y to Co/Cg** (about 1 KB vs about 10–30 B).

Caveats:

- Image content confounds this. Co and Cg typically carry less high-frequency energy than Y, and Cg less than Co here. That alone could explain part of the difference. Only a decode that shows *which bitplanes* are present can separate quantization from content.
- **Structural empties** (tiles with no data at image edges or tissue gaps) are identical across the three colours. At scale 3 in file 1, all three colours show 25.5%. These are not quantizer effects.
- The 8-byte block threshold counts empty blocks, but quantization also shrinks the non-empty ones. Only the empty count is tabulated here.

## 7. Cluster header table (0x201f), as decoded [FIC] + [code]

The base64 text decodes to a binary table. Layout, verified on all six files, which agrees with the parser at `isyntax.c:821-990`:

```
u32  total length of what follows  (= decoded size − 4)
repeat per cluster:
  u16 0xFFFE, u16 0xE000, u32 item length        // DICOM sequence-item header
    tag 0x301D,0x2019  len 12   u32 x, u32 y, u32 base_scale     // cluster origin
    tag 0x301D,0x201A  len 4    u32 cluster header template id
    tag 0x301D,0x201B  len 8    u64 cluster data offset in file
    tag 0x301D,0x201C  len 8    u64 cluster size in bytes
    tag 0x301D,0x201D  len 4·n  u32[n] codeblock offsets (relative to cluster offset)
    tag 0x301D,0x201E  len 4·n  u32[n] codeblock sizes
```

Each tag header is `u16 group, u16 element, u32 size`, the same 8-byte header as `isyntax_dicom_tag_header_t` **[code]**.

Facts:

- `n` is 63, 66 or 18, matching the cluster template for that cluster. The j-th block of a cluster is described by the j-th `UFSImageBlockHeader` in the template.
- **Cluster offsets are monotonically increasing** in all six files. Clusters are stored in file order.
- The first cluster's coordinates are `(1, 1, 0)` in files 1, 9, 10, PAS. The coordinate has a +1 bias (cluster base range starts at 1). File 2 starts at `(6145, 36865, 0)` and file 4 at `(151553, 8193, 0)`.
- **The sum of block counts equals `UFS_IMAGE_NUMBER_OF_BLOCKS`** in file 1 (33,867), which I checked per scale. The code's comment that the stored count may be lower than the header count (`isyntax.c:636`) is about other files.
- Cluster data offset range: file 1 from 4,390,072 to 88,898,174; PAS from 64,210,980 to 1,596,088,923. The gap between the header end and the first cluster is the stripped portion, which is mostly JPEG and ICC data in the original.

**Per-scale layout in file 1** (hence the 16:4:1 per cluster):

| Scale | H blocks per colour | LL |
|---|---|---|
| 0 | 8,400 (525 clusters × 16) | no |
| 1 | 2,100 | no |
| 2 | 525 | no |
| 3 | 192 | no |
| 4 | 48 | no |
| 5 | 12 | yes, 12 per colour |

Highest-scale LL blocks are 12 per colour in file 1 and PAS, 9 per colour in 9, 8 per colour in file 2.

## 8. Open questions

0. **Philips' reference code (Python and MATLAB) does not describe the v2 layout** (cluster templates, the 0x201f table, valid-data envelopes). The C code and these headers are the only sources, so section 5 and section 7 above are reverse-engineered, not vendor-documented.

1. ~~Do the predicted bitplane masks appear?~~ **Yes, verified on all codeblocks of the v1 export of this same scan** (see [Example_iSyntax_File_Claude.md](Example_iSyntax_File_Claude.md) §3). Original question: For Y at level 0 (Q = 1) mask bit 1 should be clear. For chroma at level 0 (Q = 4) mask bits 1–4 should be clear, and at level 1 (Q = 2) mask bits 1–2. This needs a full `.isyntax` file and a small tool that prints the bitmasks per codeblock. `isyntax_hulsken_decompress` is the place to hook in.
2. ~~What is the real deadzone rule?~~ **Resolved:** `abs(c) <= deadzone → 0` (Philips `quantanddeadzone.m`).
3. ~~How do Q and DZ relate?~~ **Resolved:** they are tuned table values (the Q1 and Q2 parameter sets in `quantanddeadzone.m`), not a formula.
4. ~~Does "level" in the PerLevel list mean scale?~~ **Resolved:** level 1 in the Philips code is `HL{1}`, the first DWT level of the full image, i.e. scale 0 (finest).
5. **Why does the header carry quantization for only 2–3 levels** when images have 6–9? The Philips code leaves levels beyond the table unquantized, so this is settled for the data. Why chroma lists a third (0, 0) entry when the Q2 matrices have two rows is still unexplained.
6. **Is the cluster `0x2019` third value always the base scale?** The code asserts it (`isyntax.c:917`). All sampled clusters showed 0 for the first cluster; I did not histogram it.
7. **Method 19 (Hulsken compressor v2)**: no sample uses it, so nothing in these headers shows what changes in the header between v1 and v2 files.

## 9. How this was produced

- Parsing was done with throwaway Python scripts (ElementTree for the XML, `struct` for the 0x201f table) in `/tmp/ficwork/`. They are **not kept in the repo**. Say so if you want them saved, for example as `tools/fic_summary.py`.
- Useful one-liners on the raw files (they are text, so use `grep -a`):

```bash
grep -a 'Element="0x101[bc]"' ../example_FIC_files/1.fic      # quantizer / deadzone
grep -a 'PIM_DP_IMAGE_TYPE\|NUMBER_OF_BLOCKS' ../example_FIC_files/*.fic
```
