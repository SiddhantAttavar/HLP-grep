#!/usr/bin/env bash
# Fetch the upstream BED-tree reference behind `bench --method bed_tree`.
#
# BED-tree (https://github.com/ZhangZhenjie/bed-tree) is a disk-based B+
# tree executable without a reusable library API or CMake build, so it
# cannot plug into bench directly: bench's `bed_tree` method is a native
# in-memory reimplementation of its gram-count order range query (see
# BedTreeIndex in tests/bench.cpp). This script clones the upstream
# sources into third_party/bed-tree (see .gitignore) for provenance and
# builds the reference `bedtree` binary for manual cross-checks.
#
# Usage: tests/scripts/setup_bed_tree.sh
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
DEST="$ROOT/third_party/bed-tree"

if [ ! -f "$DEST/src/Makefile" ]; then
	mkdir -p "$ROOT/third_party"
	git clone --depth 1 https://github.com/ZhangZhenjie/bed-tree.git "$DEST"
else
	echo "third_party/bed-tree already present; leaving it alone"
fi

# Build the reference binary with upstream's own Makefile (needs g++).
make -C "$DEST/src"

echo "BED-tree ready at $DEST (reference binary: $DEST/bin/bedtree)"
