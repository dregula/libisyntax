# DpsBulk production collection: census, throughput and output capacity (2026-11-27)

Tags: **[verified]** (measured here), **[open]**.

## 1. Where it is

`$SOURCEROOT/DpsBulk` is a symlink to `/mnt/Box/Donald Regula Files/DPS_decommission/Bulk` (Box Drive, 9p drvfs mount, read-only inside the sandbox). The output folder `$SOURCEROOT/DpsDICOM` points to `.../DPS_decommission/DICOM` on the same Box volume.

## 2. Size census **[verified]** (metadata only, no content read)

Listing: `c_decoder_output/dpsbulk/dpsbulk_listing.tsv` (80 s to produce).

- 12,673 `.isyntax` files in 32 numbered directories (plus one `.tis` file, 12.6 MB).
- 18.79 TB in total; mean 1.48 GB, median 1.37 GB, min 7.4 MB, max 6.77 GB, p5 233 MB, p95 3.27 GB.
- GUID file names; modification dates 2017-01 to 2019-09. Nearly all are expected to be 40x scans with 0.25 um pixels.

## 3. Read throughput **[verified]**

Cold reads, `dd bs=4M`, three 1.3 GB files: 85 to 92 MB/s per stream; four parallel streams 104 MB/s aggregate. The link ceiling is about 100 MB/s (mount msize 64 KiB). Reading 18.8 TB takes about 2.1 days of pure reading.

## 4. Consequences for the pipeline

- **CPU:** the measured transcoder costs 1.1 to 1.4 core-seconds per MB of input (the earlier model said 0.5). For 18.8 TB that is about 5,700 to 7,300 core-hours, so the read link (2.1 days) and the CPUs (16 cores: about 15 to 19 days) are not balanced; the link is not the bottleneck on a 16-core machine.
- **Output size:** the DICOM pyramid is about 2 to 2.3 times the input, i.e. about 38 to 43 TB. The output volume (B:) is 977 GB with 663 GB free, so the output **cannot fit**; it must be written in stages and moved off, or go to another target. **[open]** decision for the owner.
- **Write throughput is unmeasured [open].** The write test in the sandbox failed twice with "Read-only file system": the settings paths contain the literal text `${SOURCEROOT}/DpsDicom` (not expanded, and the case differs from `DpsDICOM`). Absolute paths under `sandbox.filesystem.allowWrite` are needed. Separately the owner reports a Box Drive problem (writing to an offline folder needs elevated privileges); the test waits for that fix.
- Do not read content of this tree casually: it is cloud-backed (downloads, network traffic). Stage files locally and ask before bulk reads.

## 5. Access control advice given

The sandbox can read 55 top-level entries under `/mnt/Box`. To block the rest: `sandbox.filesystem.denyRead` / `allowRead` / `allowWrite` with real absolute paths, `sandbox.allowUnsandboxedCommands: false`, and `permissions.blockReadsOutsideWorkingDirectories`.

## 6. Top-level batch (started by the owner, 2026-11-27 evening)

Rule from the owner: read only the **top level** of `DpsBulk`; never enter a subdirectory until told otherwise. Output goes to `/mnt/d/DpsDICOM` (local disk D:, 3.9 TB, 3.7 TB free; **[verified]** write 110 MB/s, one stream, 1 GiB with fsync). The Box write path is not used.

- Top level: 400 `.isyntax` files, 505.7 GB (mean 1,264 MB, max 4.56 GB), plus one `.tis` file (skipped) and 31 subdirectories (not entered). List: `c_decoder_output/dpsbulk/top_level.tsv`; order by size: `order_by_size.tsv`. Estimated output about 1.2 TB (2.3 x input); estimated time about 11 h on 16 cores.
- Script: `tools/transcode_tests/batch_dpsbulk.sh LIST.tsv`: stages each input from Box to `c_decoder_output/dpsbulk/work/` (the next file is prefetched during the current transcode), transcodes there, copies the layers to `/mnt/d/DpsDICOM/<guid>/`, checks sizes, writes `DONE`, deletes its local copies. Restartable (skips folders with `DONE`); failures go to `c_decoder_output/dpsbulk/logs/failed.tsv`, per-file logs in `logs/<guid>.log`, progress in `c_decoder_output/dpsbulk/batch.log`. Must run in a foreground terminal (a backgrounded child of a tool shell is killed with the shell).
- Trial (3 smallest files): 54 MB file 6 s and 57 MB file 9 s, both complete (10 layers, 88 MB output for the first). **The 25 MB file `3ef7d5ea-a6dd-4f93-a3ba-6a04f3e5bfff.isyntax` aborted: "FATAL: tile grid smaller than the coefficient array (0 0 0)"** [open]: probably a degenerate or non-WSI file; investigate after the batch (check the header with `isyntax_example FILE`).
- Not yet done: validation of the batch output (dciodvfy / wsidicom on a sample), the final tally of failures, and the capacity plan for the other 31 subdirectories (about 18.3 TB input, 38 to 43 TB output; D: cannot hold it).
