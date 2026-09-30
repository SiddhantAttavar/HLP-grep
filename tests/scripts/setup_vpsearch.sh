#!/usr/bin/env bash
# Fetch the upstream vpsearch reference behind `bench --method vpsearch`.
#
# vpsearch (https://github.com/enthought/vpsearch) ships as Python/Cython
# on top of parasail alignment scores, so it cannot plug into the C++
# bench executable directly: bench's `vpsearch` method is a native
# reimplementation of its vantage-point tree algorithm over the
# Levenshtein edit distance (see VpSearchTree in tests/bench.cpp). This
# script clones the upstream sources into third_party/vpsearch (see
# .gitignore) for provenance and installs the Python package with its
# `vpsearch` command-line utility for manual cross-checks.
#
# Usage: tests/scripts/setup_vpsearch.sh
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
DEST="$ROOT/third_party/vpsearch"

if [ ! -f "$DEST/pyproject.toml" ]; then
	mkdir -p "$ROOT/third_party"
	git clone --depth 1 https://github.com/enthought/vpsearch.git "$DEST"
else
	echo "third_party/vpsearch already present; leaving it alone"
fi

# Installs the package (Cython extension over parasail) and the `vpsearch`
# CLI into the current Python environment.
python3 -m pip install "$DEST"

echo "vpsearch ready at $DEST"
