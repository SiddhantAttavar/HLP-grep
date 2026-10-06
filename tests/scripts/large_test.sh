#!/usr/bin/env bash
# Generate a testcase/index if needed, then benchmark every method.
# Accepts the same arguments as testcase_gen.py.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
GEN="$ROOT/tests/scripts/testcase_gen.py"
BUILD="${BUILD_DIR:-$ROOT/build}"

if [[ $# -eq 1 && ( "$1" == "-h" || "$1" == "--help" ) ]]; then
    exec python3 "$GEN" --help
fi

# Read generator arguments to locate its output using the generator's own
# default-name logic, keeping this wrapper in sync with testcase_gen.py.
readarray -t OUTPUT_INFO < <(python3 - "$ROOT" "$@" <<'PY'
import argparse
import os
import sys

root, *argv = sys.argv[1:]
sys.path.insert(0, os.path.join(root, "tests", "scripts"))
import testcase_gen

parser = argparse.ArgumentParser()
parser.add_argument("--dict-size", type=int, required=True)
parser.add_argument("--num-queries", type=int, required=True)
parser.add_argument("--avg-seq-len", type=int, required=True)
parser.add_argument("--ks", type=testcase_gen.parse_ks, required=True)
parser.add_argument("--sub-rate", type=float, required=True)
parser.add_argument("--indel-rate", type=float, required=True)
parser.add_argument("--name")
parser.add_argument("--out-dir", default=testcase_gen.DEFAULT_OUT_DIR)
parser.add_argument("--seed", type=int, default=42)
args = parser.parse_args(argv)
print(os.path.abspath(args.out_dir))
print(args.name or testcase_gen.default_name(args))
PY
)
OUT_DIR="${OUTPUT_INFO[0]}"
NAME="${OUTPUT_INFO[1]}"
TESTCASE="$OUT_DIR/$NAME.txt"
INDEX="$OUT_DIR/$NAME.hlpg"

if [[ ! -f "$TESTCASE" ]]; then
    python3 "$GEN" "$@"
fi

if [[ ! -f "$INDEX" ]]; then
    cmake --build "$BUILD" -j --target save bench
    "$BUILD/tests/save" "$TESTCASE" --out "$INDEX"
fi

if [[ ! -x "$BUILD/tests/bench" ]]; then
    cmake --build "$BUILD" -j --target bench
fi

exec "$BUILD/tests/bench" \
    --method hlp_grep \
    --method wfa \
    --method dt_patricia \
    "$TESTCASE" --index "$INDEX"
