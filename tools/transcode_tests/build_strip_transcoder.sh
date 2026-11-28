#!/usr/bin/env bash
# Build src/transcode/strip_transcoder.cpp -> build-ssse3/isyntax_strip_transcode (baseline build tree, see CLAUDE.md) against the patched OpenJPEG of
# tools/transcode_tests/build_patched_openjpeg.sh (ext/, git-ignored). Needs `ninja -C build-ssse3` first (libisyntax.so).
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"; B="$REPO/build-ssse3"; SRC="$REPO/ext/openjpeg-2.5.4"; BLD="$REPO/ext/openjpeg-build"
g++ -std=c++17 -O3 -mssse3 -Wall -Wno-unused-variable -Wno-deprecated-declarations -pthread -DDEFAULT_ICC="\"$REPO/src/transcode/srgb_lcms.icc\"" -o "$B/${OUT:-isyntax_strip_transcode}" "$REPO/src/transcode/strip_transcoder.cpp" "$REPO/src/transcode/dicom_layer_writer.cpp" \
    -I"$REPO/src" -I"$REPO/src/transcode" -I/usr/include/dcmtk -I"$SRC/src/lib/openjp2" -I"$BLD/src/lib/openjp2" -L"$B" -lisyntax -Wl,-rpath,"$B" "$BLD/bin/libopenjp2.a" -ldcmdata -loflog -lofstd -lz -lm
echo "built $B/${OUT:-isyntax_strip_transcode}"
