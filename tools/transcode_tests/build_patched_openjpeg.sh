#!/usr/bin/env bash
# Build OpenJPEG 2.5.4 with the coefficient-input patch and the byte-identity test (see Pipeline_Plan_v2_Claude.md, section 2.1).
#   tools/transcode_tests/build_patched_openjpeg.sh            # clone (if needed), patch, build, compile the test
#   ext/test-bin/coeff_input_test <dir with *_rgb.raw / *_coef.raw> tissue edge glass
# Test data: .venv/bin/python3 -B -I tools/v1_analysis/make_coeff_test_frames.py <file.isyntax> <blocks.csv> <outdir>
# Everything lands in ext/ (git-ignored); the patch itself is patches/openjpeg-2.5.4-coefficient-input.patch (tracked).
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
SRC="$REPO/ext/openjpeg-2.5.4"; BLD="$REPO/ext/openjpeg-build"; PATCH="$REPO/patches/openjpeg-2.5.4-coefficient-input.patch"
mkdir -p "$REPO/ext"
[ -d "$SRC" ] || git clone --depth 1 --branch v2.5.4 https://github.com/uclouvain/openjpeg.git "$SRC"
if git -C "$SRC" apply --check "$PATCH" 2>/dev/null; then git -C "$SRC" apply "$PATCH"; echo "patch applied"
else echo "patch already applied (or does not apply cleanly: check git -C $SRC status)"; fi
cmake -S "$SRC" -B "$BLD" -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=OFF -DBUILD_CODEC=OFF -DBUILD_THIRDPARTY=OFF -DBUILD_TESTING=OFF -DCMAKE_POSITION_INDEPENDENT_CODE=ON
cmake --build "$BLD" -j
mkdir -p "$REPO/ext/test-bin"
gcc -O2 -Wall -o "$REPO/ext/test-bin/coeff_input_test" "$REPO/tools/transcode_tests/coeff_input_test.c" \
    -I"$SRC/src/lib/openjp2" -I"$BLD/src/lib/openjp2" "$BLD/bin/libopenjp2.a" -lm -lpthread
echo "built $REPO/ext/test-bin/coeff_input_test"
