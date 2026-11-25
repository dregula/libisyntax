#!/usr/bin/env bash
# Run the C decoder (libisyntax example programs) on an iSyntax file and store everything Claude needs for later
# analysis in a self-describing run folder. MEANT TO BE RUN BY THE USER in a normal terminal: Claude's own sessions
# cannot execute compiled binaries (the permission system blocks it), but Claude can read the output folder.
#
#   tools/run_c_decoder.sh <variant> [isyntax_file]
#
#   variant:
#     meson-default  existing ./builddir (Meson, plain -O3). x86 gets only __SSE2__  => SCALAR YCoCg path (truncating "/2")
#     meson-ssse3    ./build-ssse3 (Meson, -mssse3)                               => SIMD   YCoCg path (floor ">>1")
#     cmake-native   ./cmake-build (CMake adds -march=native)                     => SIMD   YCoCg path (floor ">>1")
#   isyntax_file: default ../example_isyntax_files/1.isyntax (relative to the repository root)
#
# Environment:  C_DECODER_OUTPUT  output root (default: ../c_decoder_output next to the repository)
#               FULL_TIFF=1       also run isyntax-to-tiff (LZW, lossless; large and slow; off by default)
#               EXTRA_TILES="L X Y;L X Y"  extra tiles to write, e.g. EXTRA_TILES="0 40 40;1 20 20"
# Never modifies the input file or the decoder sources. Safe to re-run: each run gets its own timestamped folder.
set -euo pipefail

VARIANT="${1:?usage: $0 <meson-default|meson-ssse3|cmake-native> [isyntax_file]}"
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
FILE="$(readlink -f "${2:-$REPO/../example_isyntax_files/1.isyntax}")"
OUTROOT="$(readlink -m "${C_DECODER_OUTPUT:-$REPO/../c_decoder_output}")"
[ -f "$FILE" ] || { echo "input file not found: $FILE" >&2; exit 1; }

case "$VARIANT" in
  meson-default)
    BUILD="$REPO/builddir"
    ninja -C "$BUILD"
    EXAMPLE="$BUILD/isyntax_example"; TOTIFF="$BUILD/isyntax-to-tiff"
    EXPECT="scalar YCoCg->RGB path (truncating /2); no -march flag, __SSE2__ only"
    FLAGS="meson buildtype=release (-O3), no -march/-mssse3"
    ;;
  meson-ssse3)
    BUILD="$REPO/build-ssse3"
    [ -d "$BUILD" ] || meson setup "$BUILD" "$REPO" -Dc_args=-mssse3
    ninja -C "$BUILD"
    EXAMPLE="$BUILD/isyntax_example"; TOTIFF="$BUILD/isyntax-to-tiff"
    EXPECT="SIMD YCoCg->RGB path (floor >>1); __SSSE3__ defined via -mssse3"
    FLAGS="meson buildtype=release (-O3) -mssse3"
    ;;
  cmake-native)
    BUILD="$REPO/cmake-build"
    cmake -S "$REPO" -B "$BUILD" -DCMAKE_BUILD_TYPE=Release
    cmake --build "$BUILD" -j
    # CMakeLists.txt sets CMAKE_RUNTIME_OUTPUT_DIRECTORY to the source directory
    EXAMPLE="$REPO/isyntax_example"; TOTIFF="$REPO/isyntax-to-tiff"
    EXPECT="SIMD YCoCg->RGB path (floor >>1); -march=native (this CPU: ssse3+avx2)"
    FLAGS="cmake Release with -march=native"
    ;;
  *) echo "unknown variant: $VARIANT" >&2; exit 1 ;;
esac
[ -x "$EXAMPLE" ] || { echo "missing executable: $EXAMPLE" >&2; exit 1; }

STAMP="$(date +%Y%m%d_%H%M%S)"
RUN="$OUTROOT/${STAMP}_${VARIANT}_$(basename "$FILE" .isyntax)"
mkdir -p "$RUN/tiles"
echo "run folder: $RUN"

{
  echo "date:            $(date -Is)"
  echo "variant:         $VARIANT"
  echo "expected colour: $EXPECT"
  echo "build flags:     $FLAGS"
  echo "compiler:        $(cc --version | head -1)"
  echo "git commit:      $(git -C "$REPO" rev-parse --short HEAD)   (dirty src files: $(git -C "$REPO" status --short src | wc -l))"
  echo "cpu isa flags:   $(grep -o -w -E 'sse2|ssse3|sse4_1|avx|avx2' /proc/cpuinfo | sort -u | tr '\n' ' ')"
  echo "input:           $FILE"
  echo "input size:      $(stat -c %s "$FILE") bytes"
  echo "input sha256:    $(sha256sum "$FILE" | cut -d' ' -f1)"
  echo "example binary:  $EXAMPLE"
  echo "---- commands ----"
} > "$RUN/RUN.txt"

run() {  # run <logname> <command...> : log the command, capture stdout+stderr, record the exit code
  local name="$1"; shift
  echo "$*" >> "$RUN/RUN.txt"
  ( time "$@" ) > "$RUN/$name.log" 2>&1 || echo "  (exit status $? for: $*)" | tee -a "$RUN/RUN.txt"
}

# 1. metadata and level/tile geometry
run levels        "$EXAMPLE" "$FILE"
# 2. associated JPEG images
run label         "$EXAMPLE" "$FILE" label "$RUN/label.jpg"
run macro         "$EXAMPLE" "$FILE" macro "$RUN/macro.jpg"
# 3. tiles: "level tile_x tile_y". Level 0 is full resolution; tile = 256 x 256 px, RGBA, written as PNG.
#    Chosen from the real file (1.isyntax, 6 levels): tissue edge, tissue interior, blank glass, edge of the scanned area,
#    and the whole top level (level 5 = 4 x 3 tiles = the entire slide at 1/32 resolution).
TILES="0 35 35;0 36 35;0 35 36;0 36 36;0 60 40;0 10 10;0 97 50;2 15 9;3 7 4"
for ty in 0 1 2; do for tx in 0 1 2 3; do TILES="$TILES;5 $tx $ty"; done; done
[ -n "${EXTRA_TILES:-}" ] && TILES="$TILES;$EXTRA_TILES"
IFS=';' read -ra TL <<< "$TILES"
for t in "${TL[@]}"; do
  set -- $t
  run "tile_L$1_${2}_$3" "$EXAMPLE" "$FILE" "$1" "$2" "$3" "$RUN/tiles/L$1_x$2_y$3.png"
done
# 4. optional full-resolution lossless conversion (large)
if [ "${FULL_TIFF:-0}" = "1" ]; then
  run to_tiff "$TOTIFF" "$FILE" "$RUN/full_lzw.tif" --compression LZW
fi

( cd "$RUN" && find . -type f ! -name SHA256SUMS.txt -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS.txt )
echo "---- files ----" >> "$RUN/RUN.txt"; ( cd "$RUN" && ls -l --time-style=long-iso . tiles >> RUN.txt )
echo "done. Tell Claude: \"C decoder output is in $RUN\""
