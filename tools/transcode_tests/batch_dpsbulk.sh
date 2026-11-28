#!/bin/bash
# Batch iSyntax -> DICOM-WSI for the TOP LEVEL of DpsBulk only (never enters subdirectories).
# Usage: batch_dpsbulk.sh LIST.tsv   (lines: name<TAB>size; processed in the order given)
# Per file: stage input from Box to local disk (next file is prefetched during the current transcode),
# transcode into a local work dir, copy to OUT/<stem>/, verify sizes, write a DONE marker, delete the local copies.
# Restartable: a file with OUT/<stem>/DONE is skipped.
set -u
ROOT=/home/regula/source/dregula
SRC=$ROOT/DpsBulk
OUT=${OUT:-/mnt/d/DpsDICOM}
W=$ROOT/c_decoder_output/dpsbulk/work
LOGS=$ROOT/c_decoder_output/dpsbulk/logs
BIN=$ROOT/libisyntax/build-ssse3/isyntax_strip_transcode
LIST=$1
mkdir -p "$W" "$LOGS"
mapfile -t NAMES < <(cut -f1 "$LIST")
prefetch() { local n=$1; [ -e "$OUT/${n%.isyntax}/DONE" ] && return; cp "$SRC/$n" "$W/$n.part" && mv "$W/$n.part" "$W/$n"; }
prefetch "${NAMES[0]}"
for i in "${!NAMES[@]}"; do
  n=${NAMES[$i]}; stem=${n%.isyntax}; t0=$(date +%s)
  if [ -e "$OUT/$stem/DONE" ]; then echo "$(date +%T) skip $stem"; continue; fi
  nxt=${NAMES[$((i+1))]:-}
  if [ -n "$nxt" ]; then prefetch "$nxt" & pf=$!; else pf=; fi
  if [ ! -e "$W/$n" ]; then prefetch "$n"; fi
  mkdir -p "$W/$stem"
  if "$BIN" "$W/$n" --threads 16 --dicom "$W/$stem" > "$LOGS/$stem.log" 2>&1; then
    mkdir -p "$OUT/$stem"
    ok=1
    for f in "$W/$stem"/*.dcm; do cp "$f" "$OUT/$stem/" && [ "$(stat -c %s "$f")" = "$(stat -c %s "$OUT/$stem/$(basename "$f")")" ] || ok=0; done
    if [ $ok = 1 ]; then date -u +%FT%TZ > "$OUT/$stem/DONE"; echo "$(date +%T) done $stem $(( $(date +%s)-t0 )) s"; else echo "$(date +%T) COPYFAIL $stem"; printf '%s\tcopy\n' "$n" >> "$LOGS/failed.tsv"; fi
  else
    echo "$(date +%T) FAIL $stem"; printf '%s\ttranscode\n' "$n" >> "$LOGS/failed.tsv"
  fi
  rm -f "$W/$stem"/*.dcm "$W/$stem"/.spool_* "$W/$n"; rmdir "$W/$stem" 2>/dev/null
  [ -n "$pf" ] && wait $pf
done
echo "$(date +%T) batch finished"
