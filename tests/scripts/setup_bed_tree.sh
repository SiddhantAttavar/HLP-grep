#!/usr/bin/env bash
# Fetch (and patch) the vendored BED-tree clone used by
# `bench --method bed_tree`.
#
# BED-tree (https://github.com/ZhangZhenjie/bed-tree) is a legacy
# disk-based B+ tree codebase without its own CMake build, so bench
# compiles its library sources (everything except main.cpp) directly
# (see tests/CMakeLists.txt). The embedded Python below applies the
# small bench patch: SBPTree/Node RangeQuery gain an optional per-match
# (id, distance) result collector (upstream returns only a hit count),
# BufferManager's hardcoded "disk.dat" backing file becomes redirectable
# per bench child, and two chatty stdout prints on the bench path are
# silenced. main.cpp still compiles, so the reference `bedtree` binary
# keeps building for manual cross-checks.
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

python3 - "$DEST" <<'EOF'
import sys

dest = sys.argv[1]
MARK = "HLP-grep bench hook"


def apply(rel, pairs):
    path = dest + "/" + rel
    with open(path, "rb") as handle:
        text = handle.read().decode("utf-8")
    text = text.replace("\r\n", "\n")  # upstream uses CRLF line endings
    if MARK in text:
        print("already patched: " + rel)
        return
    for old, new in pairs:
        if text.count(old) != 1:
            sys.exit("patch failed in %s: pattern found %d times (want 1): %r"
                     % (rel, text.count(old), old[:70]))
        text = text.replace(old, new)
    with open(path, "w", newline="") as handle:
        handle.write(text)
    print("patched: " + rel)


apply("src/Node.h", [
    ("#include <fstream>\n\nusing namespace std;",
     "#include <fstream>\n#include <utility>\n#include <vector>\n\nusing namespace std;"),
    ("\tint RangeQuery(char* queryString, int qStrLen, float threshold, \n"
     "\t\tchar* sString, char* eString);",
     "\t// HLP-grep bench hook: optional per-match result collector.\n"
     "\tint RangeQuery(char* queryString, int qStrLen, float threshold,\n"
     "\t\t char* sString, char* eString,\n"
     "\t\t std::vector<std::pair<int, float> >* out = 0);"),
])

apply("src/Node.cpp", [
    ("int Node::RangeQuery(char* queryString, int qStrLen, float threshold, \n"
     "\t\t\t\t\t char* sString, char* eString)",
     "// HLP-grep bench hook: see Node.h.\n"
     "int Node::RangeQuery(char* queryString, int qStrLen, float threshold,\n"
     "\t\t\t\t\t char* sString, char* eString,\n"
     "\t\t\t\t\t std::vector<std::pair<int, float> >* out)"),
    ("\t\t\tif( distance <= threshold)\n"
     "\t\t\t{\n"
     "\t\t\t\t//cout << strcomp->PrintString(tempString) << endl;\n"
     "\t\t\t\tcounter++;\n"
     "\t\t\t}",
     "\t\t\tif( distance <= threshold)\n"
     "\t\t\t{\n"
     "\t\t\t\t//cout << strcomp->PrintString(tempString) << endl;\n"
     "\t\t\t\tcounter++;\n"
     "\t\t\t\t// HLP-grep bench hook: record (id, distance); the\n"
     "\t\t\t\t// entry head holds the numeric dictionary id.\n"
     "\t\t\t\tif (out)\n"
     "\t\t\t\t\tout->push_back(std::make_pair(atoi(tempString), distance));\n"
     "\t\t\t}"),
    ("\t\t\tcounter += child->RangeQuery(queryString, qStrLen, threshold,\n"
     "\t\t\t\t\ttStartString, tEndString);",
     "\t\t\tcounter += child->RangeQuery(queryString, qStrLen, threshold,\n"
     "\t\t\t\t\ttStartString, tEndString, out);"),
])

apply("src/SBPTree.h", [
    ("#include <time.h>",
     "#include <time.h>\n#include <utility>\n#include <vector>"),
    ("\tint RangeQuery(char* queryString, int qStrLen, float threshold);",
     "\t// HLP-grep bench hook: optional per-match result collector.\n"
     "\tint RangeQuery(char* queryString, int qStrLen, float threshold,\n"
     "\t\t\t\t   std::vector<std::pair<int, float> >* out = 0);"),
])

apply("src/SBPTree.cpp", [
    ("int SBPTree::RangeQuery(char* queryString, int qStrLen, float threshold)",
     "// HLP-grep bench hook: see SBPTree.h.\n"
     "int SBPTree::RangeQuery(char* queryString, int qStrLen, float threshold,\n"
     "\t\t\t\t   std::vector<std::pair<int, float> >* out)"),
    ("\t\tresSize = root->RangeQuery(queryString, qStrLen, threshold, sString, eString);",
     "\t\tresSize = root->RangeQuery(queryString, qStrLen, threshold, sString, eString, out);"),
    ("\t\tcout << \"Node \" << root->GetID() << \" is new root\" << endl;",
     "\t\t// HLP-grep bench hook: silenced (bench keeps stdout clean).\n"
     "\t\t//cout << \"Node \" << root->GetID() << \" is new root\" << endl;"),
])

apply("src/BufferManager.h", [
    ("#include <fstream>\n",
     "#include <fstream>\n#include <string>\n"),
    ("\tBufferManager(int size, int constraint, int pSize);\n\t~BufferManager();",
     "\tBufferManager(int size, int constraint, int pSize);\n"
     "\t~BufferManager();\n"
     "\n"
     "\t// HLP-grep bench hook: redirect the backing file (default \"disk.dat\").\n"
     "\tstatic std::string diskPath;\n"
     "\tstatic void setDiskPath(const char* path) { diskPath = path; }"),
])

apply("src/BufferManager.cpp", [
    ("using namespace std;\n",
     "using namespace std;\n\n// HLP-grep bench hook: see BufferManager.h.\n"
     "std::string BufferManager::diskPath = \"disk.dat\";\n"),
    ("diskfile.open(\"disk.dat\", ios::out);",
     "diskfile.open(diskPath.c_str(), ios::out);"),
    ("diskfile.open(\"disk.dat\", ios::in | ios::out | ios::binary);",
     "diskfile.open(diskPath.c_str(), ios::in | ios::out | ios::binary);"),
    ("\tcout << \"Buffer manager created : \" << maxSize << \" \" << pageConstraint << endl;",
     "\t// HLP-grep bench hook: silenced (bench keeps stdout clean).\n"
     "\t//cout << \"Buffer manager created : \" << maxSize << \" \" << pageConstraint << endl;"),
])
EOF

# Build the reference binary with upstream's own Makefile (needs g++).
# A clean rebuild is required: the Makefile tracks no header
# dependencies, so objects predating the patch above would keep calling
# the unpatched RangeQuery signature.
make -C "$DEST/src" clean
make -C "$DEST/src"

echo "BED-tree ready at $DEST (reference binary: $DEST/bin/bedtree)"
