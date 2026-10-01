#!/usr/bin/env bash
# Fetch (and patch) the vendored parasail clone used by
# `tests/bench --method parasail`.
#
# Parasail (https://github.com/jeffdaily/parasail) is a SIMD C library for
# pairwise sequence alignment. Bench consumes it via add_subdirectory (see
# tests/CMakeLists.txt); the embedded Python below removes parasail's own
# ctest registration (unbuilt test binaries would otherwise fail every
# `ctest` run), mirroring tests/scripts/setup_wfa2.sh.
#
# Usage: tests/scripts/setup_parasail.sh
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
DEST="$ROOT/third_party/parasail"

if [ ! -f "$DEST/CMakeLists.txt" ]; then
	mkdir -p "$ROOT/third_party"
	git clone --depth 1 https://github.com/jeffdaily/parasail.git "$DEST"
else
	echo "third_party/parasail already present; leaving it alone"
fi

python3 - "$DEST/CMakeLists.txt" <<'EOF'
import sys

path = sys.argv[1]
with open(path) as handle:
    text = handle.read()
pairs = [
    ("ADD_TEST( NAME test_isa COMMAND test_isa )",
     "# Disabled by HLP-grep: test binaries are not built here.\n"
     "#ADD_TEST( NAME test_isa COMMAND test_isa )"),
    ("ADD_TEST( NAME test_basic COMMAND test_basic )",
     "#ADD_TEST( NAME test_basic COMMAND test_basic )"),
    ("ADD_TEST( NAME test_verify COMMAND test_verify -f ${CMAKE_SOURCE_DIR}/data/test_small_2.fasta )",
     "#ADD_TEST( NAME test_verify COMMAND test_verify -f ${CMAKE_SOURCE_DIR}/data/test_small_2.fasta )"),
]
for old, new in pairs:
    if new in text:
        continue
    if text.count(old) != 1:
        sys.exit("patch failed: pattern found %d times (want 1): %r"
                 % (text.count(old), old[:60]))
    text = text.replace(old, new)
with open(path, "w") as out:
    out.write(text)
print("patched: removed ADD_TEST registrations in CMakeLists.txt")
EOF
echo "parasail ready at $DEST"
