#!/usr/bin/env bash
# Fetch the vendored edlib clone used by `tests/bench --method edlib`.
#
# Edlib (https://github.com/Martinsos/edlib) is a C/C++ library for exact
# pairwise edit-distance alignment (Myers' bit-vector algorithm). It builds
# with its own CMakeLists, which bench consumes via add_subdirectory (see
# tests/CMakeLists.txt); the examples/utilities/test tree are switched off
# there, so no patching is needed here.
#
# Usage: tests/scripts/setup_edlib.sh
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
DEST="$ROOT/third_party/edlib"

if [ ! -f "$DEST/CMakeLists.txt" ]; then
	mkdir -p "$ROOT/third_party"
	git clone --depth 1 https://github.com/Martinsos/edlib.git "$DEST"
else
	echo "third_party/edlib already present; leaving it alone"
fi

echo "edlib ready at $DEST"
