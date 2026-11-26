# Read-only analysis tools for iSyntax files with the v1 data model (UFS_IMAGE_BLOCK_HEADER_TABLE + seektable).
# They never modify the decoder or the input file. Pure Python 3, no third-party packages.
#   python3 -I parse_v1.py <file.isyntax> <blocks.csv>        # block table + seektable -> CSV (follows Philips extract_block_header.py)
#   python3 -I exhaustive.py <file.isyntax> <blocks.csv>      # decode every block, check Quantizer/Deadzone bit patterns
#   python3 -I hulsken_v1.py <file.isyntax> <blocks.csv> <per_group> <seed>   # per-group bitmask union on a sample
# hulsken_v1.py is an independent port of the compressor-v1 decoder (isyntax.c:2551, Philips hulskendecompress.m).
#   python3 -I layout.py <blocks.csv>                         # byte layout, chunk structure, in-chunk order
#   python3 -I compare_v1_v2.py <blocks.csv> <file.fic>        # compare v1 seektable with the v2 cluster table of the matching .fic
#   <venv>/python3 -I stream_stats.py <file.isyntax> <blocks.csv>   # Hulsken stream statistics (needs numpy)
#   <venv>/python3 -I reconstruct_np.py <file.isyntax> <blocks.csv> <outdir> [downsample] [cx cy cw ch]   # full independent reconstruction (needs numpy, ~15 GB RAM)
#   (set SAVE_PLANES=1 to also write the Y/Co/Cg planes as .npy files in <outdir>)
# --- DICOM-WSI transcoding exploration (see ../../iSyntax_to_DICOM_WSI_Transcoding_Claude.md); all read-only, numpy needed unless noted ---
#   python3 -I halo_1d.py                                      # pure Python: width of the frame-edge halo, seamless vs mirrored 5/3 (1-D, synthetic)
#   <venv>/python3 -I frame_interior_test.py <file> <blocks.csv>   # real file: stored coefficients vs J2K per-frame forward transform (phase + interior-copy test)
#   <venv>/python3 -I strip_halo.py <file> <blocks.csv>        # strip halo recomputation, 3 channels, 5 frames (inputs sliced from a FULL decode)
#   <venv>/python3 -I windowed_idwt_halo.py <file> <blocks.csv>    # same, but strips come from a WINDOWED top-down IDWT; per DICOM layer 0..5; block/chunk footprint
#   python3 -I count_chunks.py <blocks.csv>                    # pure Python: chunks touched per scale-0 tile (conservative dependency closure)
#   python3 -I rct_ycocg_vs_rct.py                             # pure Python: YCoCg-R (fix / floor) vs DICOM YBR_RCT, rounding differences
#   python3 -I ll_probe.py <file> <blocks.csv>                 # pure Python: stitch the 12 top-level LL blocks per colour, rounding probe (inconclusive)
# count_chunks.py and ll_probe.py hard-code the geometry of 1.isyntax (origin 77715/20614, 6 levels, 128 x 96 tiles).
# Running any script that imports hulsken_v1 regenerates the tracked __pycache__/*.pyc; add -B (python3 -B -I ...) to avoid that (-I ignores PYTHONDONTWRITEBYTECODE).
#   <venv>/python3 -B -I dicom_tile_tests.py <file> <blocks.csv> [D=5] [sizes=256,512,1024]   # DICOM tile grids, absolute coordinates, per-tile tests (needs ~14 GB RAM, ~2 min per size set)
#   <venv>/python3 -B -I colour_ebcot_probe.py <file> <blocks.csv>    # wavelet-domain YCoCg->YBR_RCT error on real frames + per-scale sparsity/entropy of Y, Co, Cg, Cb, Cr (EBCOT-relevant)
#   <venv>/python3 -B -I colour_bias_by_levels.py <file> <blocks.csv> # how the wavelet-domain chroma error grows with the number of decomposition levels
#   <venv>/python3 -B -I colour_sparsity_cost.py <file> <blocks.csv>  # sparsity of Cb/Cr: wavelet-domain conversion vs exact pixel-domain conversion
#   <venv>/python3 -B -I colour_dc_fix.py <file> <blocks.csv>         # residual error after a constant correction of the LL band
#   <venv>/python3 -B -I midpoint_probe.py <file> <blocks.csv>        # pixel error if dropped (structurally zero) LSB planes meet a decoder with reconstruction r = 1/2; ring-free code-block shares
#   <venv>/python3 -B -I openjpeg_convention_check.py                      # (needs Pillow+OpenJPEG) does OpenJPEG's forward 5/3 match our NumPy convention? sparse-coefficient size test; mct check
#   <venv>/python3 -B -I openjpeg_ebcot_sizes.py <file> <blocks.csv> [2048]    # real EBCOT sizes of 3 colour routes x frame sizes 256..2048 on real regions, versus iSyntax bytes
#   <venv>/python3 -B -I openjpeg_empty_plane_cost.py <file> <blocks.csv>      # real cost of coding the always-zero low planes (upper bound of the saving from dropping them)
#   <venv>/python3 -B -I omission_size_estimate.py <file> <blocks.csv> [outdir] [sizes]   # per-tile stats + real OpenJPEG size of every tile of every layer (16 processes, about 4 minutes), writes tile_stats_F*_L*.csv
#   python3 -B -I omission_summary.py [dir]                                           # savings of omission strategies A / B(robust T) / strict S from those CSVs
#   <venv>/python3 -B -I omission_review_images.py <file> <blocks.csv> [dir] [F]      # PNGs for the visual review: masks, single-fill previews, diffs x16, full-resolution crops
#   <venv>/python3 -B -I rgb_route_probe.py <file> <blocks.csv> [A] [B] [C]   # DICOM Photometric Interpretation RGB (mct=0): wavelet-domain error (A), sparsity of R/G/B (B), OpenJPEG sizes versus YBR_RCT (C); needs blocks.csv from parse_v1.py
#   <venv>/python3 -B -I compare_sdk_tiff_tiles.py <c_decoder run folder> <iSyntax_test_materials> [PAS]   # C-decoder tile PNGs vs the Philips SDK LZW TIFF levels 2/3/5: which YCoCg->RGB inverse the SDK uses
#   <venv>/python3 -B -I speck_size_analysis.py <file> <blocks.csv> [stats_dir] [T=4] [F=256] [out.csv]   # size (px, um) of the specks inside tiles that the robust-range omission rule would delete
#   <venv>/python3 -B -I make_coeff_test_frames.py <file> <blocks.csv> <outdir>   # RGB + YBR_RCT coefficient (Mallat layout) test frames for the patched-OpenJPEG test, see tools/transcode_tests/
