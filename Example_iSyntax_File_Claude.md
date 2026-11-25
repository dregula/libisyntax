# Example iSyntax file (`example_isyntax_files/1.isyntax`): analysis

File: `/mnt/u/regula/source/dregula/example_isyntax_files/1.isyntax`, 88,928,078 bytes, dated 2020-05-11. This is the first real file with codeblocks, so it is the first time the Quantizer/Deadzone predictions in the earlier documents could be tested on data. Companions: [Philips_Compression_Code_Samples_Claude.md](Philips_Compression_Code_Samples_Claude.md) (the algorithm), [Example_FIC_Files_Claude.md](Example_FIC_Files_Claude.md) (the sample headers), [Image_Format_Code_Samples_Claude.md](Image_Format_Code_Samples_Claude.md) (Philips' file reader), [CLAUDE.md](CLAUDE.md).

**Evidence labels:** **[ran]** = I executed it; **[derived]** = computed from the results; **[inferred]** = deduction. Scripts are in [tools/v1_analysis/](tools/v1_analysis/). Nothing was written next to your file, and the decoder in `src/` was not touched.

## 1. What the file is

| Item | Value | Source |
|---|---|---|
| Interface version (`PIM_DP_UFS_INTERFACE_VERSION`) | **5.0**, i.e. data model **v1** | header [ran] |
| Derivation | `PHILIPS UFS V1.8.6614 \| Quality=2 \| DWT=1 \| Compressor=16` | header |
| Acquired / scanner | 2019-06-12 23:16:59, serial FMT0132 | header |
| Header (XML) size | 2,422,997 bytes before the `\r\n\x04` marker | [ran] |
| Seektable | tag `301D,2015`, 1,967,040 bytes, **49,176 entries** (40 bytes each) | [ran] |
| Block header table | attribute 0x2014, base64, 1,625,616 bytes decoded, 33,867 items | [ran] |
| First codeblock offset | 4,390,072 | [ran] |
| WSI | 25,600 × 21,504 padded, 6 levels, origin (77715, 20614), 128×128 blocks | [ran] |
| Quantizer/deadzone table | `(1,2) (0,0)` and `(4,16) (2,8) (0,0)`, identical to the sample FICs | header |

**This is the same scan as `1.fic`.** Acquisition time, serial number, block count (33,867), derivation string and the quantizer table all match. `1.fic` is the same slide written with the **v2 data model (100.4)**. Beyond the header version, the two describe identical data:

- the v1 seektable gives **33,867 blocks with offsets 4,390,072 … 88,927,897**, and the last block ends exactly at end-of-file (88,928,078);
- the total of all block sizes is **84,267,078 bytes**, exactly the total I got from `1.fic`'s cluster header table (Example_FIC_Files_Claude.md §6.4);
- the first offset, 4,390,072, is also the first cluster offset in `1.fic`.

So the stripped `1.fic` headers are a faithful description of this file's data, and the v2 cluster table and the v1 block-header-plus-seektable tables agree block for block in count, total size, and first offset. (I did not compare every individual offset between the two tables.)

This also resolves the open "version-1" question in the earlier notes: **this file is the version-1 file**, and `1.fic` is its v2 counterpart. The two differ in header layout, not in the pixel data.

## 2. Philips' reference Python on a real file **[ran]**

The Python in `ImageFormatCodeSamples` was run read-only (functions called directly, no output file written):

- `isyntax_metadata`: levels 6, valid size 25,222 × 21,126, codeblock 128 × 128, 33,867 blocks. Same as the v2 header.
- `extract_block_header` (`get_block_header_table`, `read_seektable`, `fill_block_header`): **33,867 blocks, all with a nonzero offset; offsets and sizes identical to my independent parse for every block.**
- Delimiter found at byte 2,423,000 (= 2,422,997 + 3), in 0.2 s for this file.

So the Philips samples work as designed on a v1 file, and the limitation recorded earlier (they cannot read v2 headers) is the only one found.

## 3. Quantizer and deadzone on real codeblocks **[ran]**

I decoded **all** 33,867 codeblocks with an independent pure-Python port of the compressor-v1 decoder (`tools/v1_analysis/hulsken_v1.py`, written from `isyntax.c:2551` and the Philips MATLAB). Consistency checks on the port itself:

- **0 decode failures**, and for every non-empty block the decoded length equals `planes × 2048 + mask bytes` (the same invariant the C code asserts).
- Block counts and sizes agree with the headers.

### 3.1 The prediction held for every block

Prediction (Philips_Compression_Code_Samples_Claude.md §3.3), in compressor-v1 mask numbering where bit 0 is the sign plane and bit `k` is magnitude bit `k − 1`: mask bits `1 … Q` must be clear.

| Group | (Q, DZ) | Non-empty blocks | Blocks with a forced-clear bit set |
|---|---|---|---|
| Y, scale 0 | (1, 2) | 7,938 | **0** |
| Y, scale 1 (and coarser) | (0, 0) | 2,009 (scale 1) | 0 (nothing forced) |
| Co, scale 0 | (4, 16) | 1,535 | **0** |
| Cg, scale 0 | (4, 16) | 618 | **0** |
| Co, scale 1 | (2, 8) | 1,919 | **0** |
| Cg, scale 1 | (2, 8) | 1,405 | **0** |
| Scale 2 and above, all colours | none | all | 0 |

**Total violations over all colours, scales and bands: 0**, across all 17,626 non-empty blocks (each checked on all of its subbands; the other 16,241 blocks are the 8-byte empty ones). The union of mask bits per group shows the pattern directly (12 sampled blocks per group; `1` = plane present, bit 0 first):

```
Y  s0 HL  1.111...........      bit 1 clear          (Q=1)
Y  s1 HL  1111111.........      nothing forced clear (Q=0)
Co s0 HL  1....11.........      bits 1-4 clear       (Q=4); only magnitude bits 4 and 5 are used
Co s1 HL  1..1111.........      bits 1-2 clear       (Q=2)
Co s2 HL  11111111........      nothing forced clear
Y  s5 LL  .111111111......      sign plane absent (no negative luma in LL), LL never quantized
```

### 3.2 The deadzone and the "multiple of 2^Q" property **[ran]**

For 6 randomly chosen non-empty blocks per group I unpacked the actual coefficient values:

| Group | (Q, DZ) | Smallest nonzero \|c\| | Largest \|c\| seen | Values not a multiple of 2^Q |
|---|---|---|---|---|
| Y s0 | (1, 2) | **2** | 14 | 0 |
| Co/Cg s0 | (4, 16) | **16** | 32 | 0 |
| Co/Cg s1 | (2, 8) | **8** | 72 | 0 |
| everything else | (0, 0) | 1 | up to 262 | n/a |

The smallest surviving magnitude equals the deadzone in each case, exactly as `abs(c) <= dz → 0` then `fix(c/2^Q)·2^Q` predicts (3 → 2, 17 → 16, 9 → 8). This is a sample of coefficients, not all of them. The mask test above is exhaustive.

**Conclusion:** the Philips encoder algorithm in `quantanddeadzone.m` is what was used to write this file, and it is visible in the bitplanes exactly as predicted. The decoder needs no dequantization, so nothing in `isyntax_hulsken_decompress` has to change. This answers goal 4.

### 3.3 Empty blocks, exactly **[derived from the seektable]**

Blocks of 8 bytes or less (zero coefficients in every subband):

| | Y | Co | Cg |
|---|---|---|---|
| scale 0 (of 8,400) | 462 (5.5%) | 6,865 (81.7%) | 7,782 (92.6%) |
| scale 1 (of 2,100) | 91 (4.3%) | 181 (8.6%) | 695 (33.1%) |
| scale 2 (of 525) | 0 | 0 | 0 |
| scale 3 (of 192), scale 4 (of 48), scale 5 (of 12) | 49, 6, 0 | 49, 6, 0 | 49, 6, 0 |

These match the earlier cluster-table numbers, which now stand confirmed. The scale 3–4 empties are the same across the three colours, so they are structural (missing tiles). The scale 0–1 chroma empties are the deadzone. In tissue regions it is the deadzone that makes up most of the chroma empties, so the MATLAB comment that 8-byte blocks are "non-tissue background" understates them (Image_Format_Code_Samples_Claude.md §4.3).

In a non-empty chroma scale-0 block, often only some subbands carry data: for example Cg scale 0 HH had an empty mask in the blocks I sampled while HL and LH did not.

## 4. Facts this file settles

- The earlier "version-1" confusion: this file is v1 (5.0), `1.fic` is v2 (100.4), same scan.
- Philips' Python works on real v1 data, with offsets identical to an independent parse.
- The v1 block-ID formula is validated end to end: all 33,867 lookups land on distinct, in-bounds seektable entries whose sizes sum to the file's data size.
- The compressor-v1 mask semantics (`1…Q` clear), the `(Q, DZ)` table, and the Y / chroma column mapping are correct on real data.

## 5. Limits, and what would add more

- **I could not run a compiled decoder.** I compiled a small C tool that links `libisyntax.so` and calls `isyntax_hulsken_decompress` (it built), but running the binary was blocked by the permission system, so I did not retry. The results above come from the independent Python port, which also acts as a cross-check of the C logic. If you want the C decoder confirmed against the same data, I can write that tool as a repository file and you can run it, or you can allow the run.
- **Only compressor v1** (method 16) is exercised. The v2 compressor (method 19: bitmasks in the header, MSB-first planes, bitplane seektable) is untested, as no sample uses it.
- **No v2 codeblocks yet.** `1.fic` and the other FICs are headers only, so the v2 cluster layout is still verified only through `1.fic`'s table totals against this file.
- **Pixels:** section 9 reconstructs the whole image independently (numpy) and checks it visually and statistically, but I could not compare it pixel-for-pixel with libisyntax's own output, because running the compiled decoder was blocked.

## 6. About truncated files

You asked whether truncated files would be needed and how large a file I can handle. For this file, none of that was necessary:

- I read files straight from disk with scripts, so **size is not a limit for me**. The 89 MB file was fully decoded in about 40 seconds with a pure-Python decoder (roughly 2 MB/s). What matters is that I do not load whole files into the conversation, and I never need to.
- A multi-gigabyte file also works for header, seektable and **sampled** codeblock analysis, which takes seconds. Decoding every block of a 1.5 GB file would take on the order of 10–15 minutes with the Python decoder.
- If you do cut a file, **keep a prefix**: the complete header up to the end-of-table marker, the whole seektable (or, for v2, the whole header, which holds the cluster table), and then as many whole codeblocks as you like. All offsets are absolute, so any block that lies entirely inside the prefix remains decodable. For this file the header plus seektable is the first 4,390,072 bytes. A prefix of 4.4 MB plus 20–50 MB of codeblocks would already be plenty.
- The most useful additional file would be a **v2 (data model 100.x) file with codeblocks**, ideally a small one, and one using compressor method 19 if you have it.

## 7. Byte layout of the file **[ran]**

```
0              XML header, 2,422,997 bytes
2,422,997      EOT  0D 0A 04 (3 bytes)                      -> 2,423,000
2,423,000      seektable: tag 301D,2015 + u32 size 1,967,040  (49,176 entries x 40 bytes) -> ends 4,390,048
4,390,048      FFFE,E0DD 00000000      sequence delimiter (8 bytes)
4,390,056      FFFE,E000 00000000      an empty item (8 bytes)
4,390,064      FFFE,E000 <u32 size>    item header of codeblock 0, payload at 4,390,072
   ...         every block is [8-byte item header][payload], back to back
88,928,078     end of file = end of the last payload
```

- There is **no `7FE0,0010` pixel-data tag** in this file (the spec lists it as optional).
- **The codeblock section is one contiguous run**: all 33,866 gaps between consecutive payloads are exactly 8 bytes (the next item header). The seektable's `size` and `offset` refer to the payload only, not the 8-byte header.
- **Chunks (clusters):** the file holds 525 chunks of 63 blocks (scales 0, 1, 2; 16 + 4 + 1 blocks per colour) followed by 12 chunks of 66 blocks (scales 3, 4, 5, with the top-level LL). These counts equal the 537 clusters of `1.fic` (templates 0 and 1).
- **Order inside a chunk, from the real offsets:** colour-major (Y, then Co, then Cg). Within a colour: the coarsest block first, then its 2×2 children (x varying fastest), then the 4×4 grandchildren (x fastest). Example, Y of the first chunk: scale 2 block `(77725,20624)`, then scale 1 blocks at x = 77719, 78231 for y = 20618, then y = 21130, then 16 scale-0 blocks stepping by 256 px. This matches the example table in `isyntax.c:2435` and the AI notes' claim about the 1 + 4 + 16 layout.
- **Block coordinates:** the first scale-0 block is at `(77716, 20615)`, i.e. origin `(77715, 20614)` + 1 = `(3<<0) − 2`, as for every sample FIC.

### 7.1 Complete cross-check of the v1 and v2 tables **[ran]**

I compared **every block** of this file's v1 seektable with the v2 cluster header table of `1.fic`: for all 33,867 blocks, **(colour, scale, coefficient, offset, size) is identical**, as multisets, with 0 blocks only in one table. This closes the gap noted in section 1 (where only totals and the first offset were compared). So `1.fic`'s v2 header, including its cluster templates, describes this file's data exactly.

## 8. Statistics of the Hulsken streams **[ran]**

All 17,626 non-empty blocks, decoded with the Python port (`tools/v1_analysis/stream_stats.py`, numpy for counting):

- **Overall size:** the non-empty blocks hold 1,730 MB of 16-bit coefficients and compress to 84.1 MB (**20.6 : 1**). Dropping the empty bitplanes alone brings that to 339 MB, so **the entropy stage (Huffman + zero runs) contributes a further 4.0 : 1**.
- **Where the bytes go (compressed MB):**

  | | Y | Co | Cg |
  |---|---|---|---|
  | scale 0 | 8.4 | 0.08 | 0.03 |
  | scale 1 | **31.0** | 2.9 | 1.2 |
  | scale 2 | 9.2 | 10.7 | 9.9 |
  | scale 3 | 2.6 | 2.7 | 2.6 |
  | scale 4 | 0.7 | 0.7 | 0.7 |
  | scale 5 + LL | 0.3 | 0.2 | 0.2 |

  The largest single item is **Y at scale 1 (37% of the file)**. Y at scale 0, though almost four times as many non-empty blocks (7,938 against 2,009), costs a quarter as much: its coefficients are tiny (largest magnitude 14 in the sampled blocks) and the deadzone and Q = 1 remove most of them. At scale 0 the chroma is almost free: 0.11 MB for 2 channels.
- **The entropy stage works hardest where the data is sparsest:** the ratio serialized/compressed is 15 : 1 for Y scale 0, 112 : 1 and 105 : 1 for Co and Cg scale 0 (99.7% of the serialized bytes are zero), and only about 2 : 1 for scales 2–5 (46–68% zero bytes).
- **Huffman trees:** 2 to 256 leaves per block (median 19). The **longest code is 15 bits** (the C decoder assumes at most 16: "a Huffman code of at most 16 bits"). The code-length histogram peaks at 8–12.
- **Zero-run parameters:** the zero-run symbol is the least frequent non-zero byte, which in practice is a small value: 3 (29% of blocks), 1 (14%), 5 (7%), 2, 11, 7. The counter width is mostly 3 or 4 bits (75% of blocks), with a second cluster at 11–12 bits for very sparse blocks.
- **Bit density per plane:** in the unquantized bands the density of 1 bits falls with plane significance: for Y at scale 2: 0.40 (magnitude bit 0), 0.17, 0.11, 0.07, 0.03, 0.005, … and the sign plane 0.20. At Y scale 0 the planes in use are magnitude bits 1 and 2 (densities 0.0065 and 0.013); everything above is below 0.002, and the sign plane is 0.007. At chroma scale 0, planes 0–3 are never stored (Q = 4), and the surviving planes (magnitude bits 4–6) have density about 0.0005. This is the quantizer and deadzone seen directly in the planes.
- **Planes stored per subband-block** (share of blocks that store the plane): in the unquantized bands the three lowest planes are almost always stored; the share falls to 0.5–0.6 around magnitude bit 4 and to 0 beyond bit 8. The largest magnitude in any intermediate array of the reconstruction is 282 (9 bits), so the upper magnitude planes are never stored.

## 9. Independent reconstruction of the image **[ran]**

`tools/v1_analysis/reconstruct_np.py` (numpy; about 3 minutes and under 15 GB of memory): decode every block, assemble the subband arrays by block coordinates, run the integer 5/3 inverse level by level, convert to RGB. Its geometry is taken from the Philips MATLAB, **not** from libisyntax: H coefficient `n` of scale `s` sits at padded position `(3<<s) − 2 + 2^(s+1)·n`, LL at that `+ 2^s`; the low band is at odd positions; each level inverts horizontally then vertically; the finer LL is the previous output shifted by one sample. The inverse is the standard JPEG 2000 reversible 5/3 lifting (`s = L − ((H + H' + 2) >> 2)`, `H + ((s + s') >> 1)`).

**Result: a clean image.** The full-size reconstruction is 25,222 × 21,126 pixels (the valid region, padded size minus 2 × 189) and shows an H&E tissue section with correct colour and shape on a grey glass background. A full-resolution crop (1024 × 768 px at 0.25 µm/px) shows individual nuclei, with no blocking, ringing or seams. Because the geometry came from a different source than libisyntax, this is an independent confirmation of: the block placement formulas, the parity (H first, low band at odd positions), the H-then-V order, the 8-neighbour halo handling being unnecessary for a whole-level inverse, and the decoder port.

**Quantitative checks on the full reconstruction:**

| Check | Result |
|---|---|
| **Tile seams** (mean \|ΔY\| of neighbouring pixels by position modulo the block pitch 256, 512, …) | max / median = 1.02 at pitch 256 and 1.03 at pitch 512, with no peak at any block boundary. The one apparent outlier (ratio 1.3) is a single column and a single row, each a step of about 24.5 grey levels from glass (Y ≈ 230) to a flat fill (Y ≈ 255), at the **edge of the scanned area** (valid-region column 24,751, row 20,341). It is real image content, not a decode artefact. |
| Largest intermediate magnitude in any level | **282**, so the 16-bit coefficient type (`DWT_COEFF_BITS 16`) never overflows on this slide; the SIMD 16-bit wrap noted in earlier documents does not occur here |
| Negative Y after the IDWT | **0 pixels** (Y ranges 48 to 257), so the `abs(Y)` step (`isyntax.c:2399`) changes nothing on this slide |
| RGB below 0 before clamping | **0 pixels**, so the scalar YCoCg path's missing lower clamp never matters here |
| RGB above 255 before clamping | 1,890,104 pixels (0.35%), clamped by every path |
| Co range, Cg range | −111 … 134, −132 … 61 |

### 9.1 Rounding in the YCoCg-R inverse: a measured, large effect

The decoder has two colour paths (Philips_Compression_Code_Samples_Claude.md §6): the vectorised path uses `>> 1` (floor) and the scalar path uses `/ 2` (truncation toward zero), and Philips ships a reference of each kind. On the full reconstruction:

- **64.5% of all pixels (343,625,443 of 532,839,972) differ between the two** when the same Y, Co, Cg values are converted.
- The differences are small but not tiny: R and B differ by 0 to 2 levels (mean 0.80), G by 0 to 1 (mean 0.37).
- The cause is how common negative odd chroma values are in this H&E slide: 53% of pixels have Cg < 0 (36.5% negative and odd), 75% have Co < 0 (43% negative and odd).

So the earlier statement that the scalar path matters only for "tail pixels" understates it. **Which path runs decides the output of most pixels, by up to two levels.** **Which path runs depends on the build system (checked in `build/compile_commands.json`):** the CMake build uses `-march=native`, so `__SSSE3__` is defined and the vectorised floor path handles every 256-wide row; the Meson build adds no `-march` flag on x86 (the existing `build/` used plain `-O3`), so only `__SSE2__` is defined and the **truncating scalar path converts every pixel**. The same file therefore decodes to different pixels with `meson` and with `cmake` builds. `CLAUDE.md` describes how to run both and store the output. This file cannot say which rounding Philips' own viewer uses, since the transform is a bijection for either rounding; that needs the Philips viewer or its reference output.

### 9.2 What the reconstruction does not show

- It is **not** compared pixel for pixel with libisyntax's output (the compiled decoder could not be run), so exact agreement of the C pipeline remains unconfirmed, apart from the shared geometry, which is now supported.
- It does not use the 8-neighbour halo scheme of the C code. A whole-level inverse with mirrored array ends gives the same interior, because the padding keeps edge effects out of the valid region.
- Only compressor v1 data.

## 10. Housekeeping notes

- `.venv/` (Python 3.14 with numpy) sits inside the repository folder; `git status` does not list it, because it carries its own `.gitignore`. The scratch arrays I wrote to `/mnt/u/regula/source/dregula/_scratch_recon/` (3 GB) were deleted at the end of this session.
- Tools added in [tools/v1_analysis/](tools/v1_analysis/): `layout.py`, `compare_v1_v2.py`, `stream_stats.py`, `reconstruct_np.py` (plus the earlier `parse_v1.py`, `hulsken_v1.py`, `exhaustive.py`). All are read-only with respect to the input files.
