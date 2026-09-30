#!/usr/bin/env bash
# Fetch (and patch) the vendored hstree clone used by
# `bench --method hstree`.
#
# hstree lives inside
# https://github.com/TsinghuaDatabaseGroup/Similarity-Search-and-Join
# (hstree/src, exact threshold edit-distance search with the `-r` driver:
# hierarchical-substring filtering plus exact verification). The embedded
# Python below applies the small bench patch: hsearch() gains an optional
# per-match (id, distance) result collector (upstream only counts hits in
# the simres global), and two chatty stdout prints on the bench path are
# silenced. main.cc still compiles, so the reference `hstree` binary keeps
# building for manual cross-checks.
#
# Usage: tests/scripts/setup_hstree.sh
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
DEST="$ROOT/third_party/Similarity-Search-and-Join"

if [ ! -f "$DEST/hstree/src/hstree.h" ]; then
	mkdir -p "$ROOT/third_party"
	git clone --depth 1 \
	    https://github.com/TsinghuaDatabaseGroup/Similarity-Search-and-Join.git \
	    "$DEST"
else
	echo "third_party/Similarity-Search-and-Join already present; leaving it alone"
fi

python3 - "$DEST/hstree/src" <<'EOF'
import sys

dest = sys.argv[1]


def apply(rel, pairs):
    path = dest + "/" + rel
    with open(path, "rb") as handle:
        text = handle.read().decode("utf-8")
    if all(new in text for _, new in pairs):
        print("already patched: " + rel)
        return
    for old, new in pairs:
        if new in text:
            continue  # this hunk is already applied
        if text.count(old) != 1:
            sys.exit("patch failed in %s: pattern found %d times (want 1): %r"
                     % (rel, text.count(old), old[:70]))
        text = text.replace(old, new)
    with open(path, "w", newline="") as handle:
        handle.write(text)
    print("patched: " + rel)


apply("hstree.h", [
    ("void hsearch(const string& query, int tau);",
     "// HLP-grep bench hook: optional per-match result collector.\n"
     "void hsearch(const string& query, int tau,\n"
     "             vector<pair<int, int> >* out = 0);"),
])

apply("hsimSearch.cc", [
    ("void hsearch(const string& query, int tau)",
     "// HLP-grep bench hook: see hstree.h.\n"
     "void hsearch(const string& query, int tau,\n"
     "             vector<pair<int, int> >* out)"),
    ("                        int veres = caled(query,dict[id],tau);\n"
     "                        if(veres <= tau)\n"
     "                            simres++;\n",
     "                        int veres = caled(query,dict[id],tau);\n"
     "                        if(veres <= tau)\n"
     "                        {\n"
     "                            simres++;\n"
     "                            // HLP-grep bench hook: record (id, distance).\n"
     "                            if (out)\n"
     "                                out->push_back(make_pair(id, veres));\n"
     "                        }\n"),
    ("                    \t          int veres = caled(query,dict[id],tau);\n"
     "                        \t           if(veres <= tau)\n"
     "                        \t            simres++;\n",
     "                    \t          int veres = caled(query,dict[id],tau);\n"
     "                        \t           if(veres <= tau)\n"
     "                        \t           {\n"
     "                        \t               simres++;\n"
     "                        \t               // HLP-grep bench hook: record (id, distance).\n"
     "                        \t               if (out)\n"
     "                        \t                   out->push_back(make_pair(id, veres));\n"
     "                        \t           }\n"),
])

apply("prepare.cc", [
    ('  cout<<"initial finished!"<<endl;',
     '  // HLP-grep bench hook: silenced (bench keeps stdout clean).\n'
     '  //cout<<"initial finished!"<<endl;'),
    ('   cout<<"Index size: "<<(totalsize/(1024*1024))<<" MB"<<endl;',
     '   // HLP-grep bench hook: silenced (bench keeps stdout clean).\n'
     '   //cout<<"Index size: "<<(totalsize/(1024*1024))<<" MB"<<endl;'),
])

apply("free.cc", [
    ('   cout<<"delete finished!"<<endl;',
     '   // HLP-grep bench hook: silenced (bench keeps stdout clean).\n'
     '   //cout<<"delete finished!"<<endl;'),
])
EOF

# Build the reference binary with upstream's own compile recipe (needs g++).
# compile.sh removes its .o files, so a plain re-run is a full rebuild.
(cd "$DEST/hstree/src" && sh compile.sh)

echo "hstree ready at $DEST/hstree/src (reference binary: $DEST/hstree/src/hstree)"
