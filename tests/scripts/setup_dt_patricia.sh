#!/usr/bin/env bash
# Fetch (or refresh) the vendored DT-Patricia clone used by
# `tests/bench --method dt_patricia`.
#
# DT-Patricia is not tracked here (see .gitignore); this script does a
# shallow clone into third_party/DT-Patricia. It is header-only and its
# examples/tests are off when used as a subdirectory, so no patching is
# needed.
#
# Usage: tests/scripts/setup_dt_patricia.sh
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
DEST="$ROOT/third_party/DT-Patricia"

if [ ! -f "$DEST/CMakeLists.txt" ]; then
	mkdir -p "$ROOT/third_party"
	git clone --depth 1 https://github.com/kurrrru/DT-Patricia.git "$DEST"
else
	echo "third_party/DT-Patricia already present; leaving it alone"
fi

echo "DT-Patricia ready at $DEST"
