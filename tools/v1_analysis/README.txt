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
