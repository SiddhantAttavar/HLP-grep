#!/usr/bin/env bash
# Fetch (or update) the vendored WFA2-lib clone used by `tests/bench --method wfa`.
#
# WFA2-lib is not tracked here (see .gitignore); this script does a shallow
# clone into third_party/WFA2-lib and removes WFA2's own ctest registration
# from its CMakeLists (its unit tests need tool binaries this repo does not
# build, and the stray wfa2lib ctest entry would fail every `ctest` run).
#
# Usage: tests/scripts/setup_wfa2.sh
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
DEST="$ROOT/third_party/WFA2-lib"

if [ ! -f "$DEST/CMakeLists.txt" ]; then
	mkdir -p "$ROOT/third_party"
	git clone --depth 1 https://github.com/smarco/WFA2-lib.git "$DEST"
else
	echo "third_party/WFA2-lib already present; leaving it alone"
fi

# Drop the unconditional `enable_testing(); add_wfa_test()` block so the
# vendored build only contributes the wfa2_static library.
python3 - "$DEST/CMakeLists.txt" <<'EOF'
import re
import sys

path = "/".join(sys.argv[1:2]) if len(sys.argv) > 1 else None
path = sys.argv[1]
with open(path) as handle:
    text = handle.read()
patched = text.replace("""
function(add_wfa_test)
  add_test(
    NAME wfa2lib
    COMMAND ./tests/wfa.utest.sh ${CMAKE_CURRENT_BINARY_DIR}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )
endfunction()

add_wfa_test()""",
"""# Disabled by HLP-grep: unit-test binaries are not built here.""")
if patched != text:
    with open(path, "w") as out:
        out.write(patched)
    print("patched: removed add_wfa_test in CMakeLists.txt")
else:
    print("no patch needed (already applied?)")
EOF
echo "WFA2-lib ready at $DEST"
