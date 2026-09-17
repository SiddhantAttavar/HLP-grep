#!/usr/bin/env bash
# Test HLP-grep on real-world highly similar DNA datasets.
#
#   tests/scripts/test_datasets.sh [--quick] [--dataset NAME]...
#   tests/scripts/test_datasets.sh --full [--validate] [--dataset NAME]...
#
# --quick (default): small profiles, solve + validate every derived case.
# --full: larger profiles (still naive-solver bounded); without --validate it
#   only generates. --validate runs solve + test on full cases too.
# Logs go to tests/testcases/datasets/logs/ (git-ignored).
set -u

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
DATA="$ROOT/tests/testcases/datasets"
DERIVED="$DATA/derived"
LOGS="$DATA/logs"
BUILD="$ROOT/build"
MODE=quick
VALIDATE=yes
DATASETS=""

while [ $# -gt 0 ]; do
    case "$1" in
        --quick) MODE=quick; shift ;;
        --full) MODE=full; VALIDATE=no; shift ;;
        --validate) VALIDATE=yes; shift ;;
        --dataset) DATASETS="$DATASETS --dataset $2"; shift 2 ;;
        -h|--help)
            sed -n '2,12p' "$0"; exit 0 ;;
        *) echo "unknown arg: $1" >&2; exit 1 ;;
    esac
done

mkdir -p "$DERIVED" "$LOGS"
STAMP="$(date +%Y%m%d-%H%M%S)"
LOG="$LOGS/${MODE}-${STAMP}.log"

# Larger but still naive-solver bounded generation profile.
if [ "$MODE" = full ]; then
    GEN_ARGS="--dict-size 20 --num-queries 8"
else
    GEN_ARGS=""
fi
# shellcheck disable=SC2086
{
echo "== test_datasets.sh mode=$MODE $(date -u +%FT%TZ)"
echo "-- generate"
# shellcheck disable=SC2086
python3 "$ROOT/tests/scripts/generate_dataset_tests.py" $DATASETS $GEN_ARGS
echo "-- build"
cmake -S "$ROOT" -B "$BUILD" -DCMAKE_BUILD_TYPE=Release > /dev/null
cmake --build "$BUILD" -j
if [ "$VALIDATE" = yes ]; then
echo "-- solve + validate"
for d in "$DERIVED"/*/; do
    [ -d "$d" ] || continue
    "$BUILD/tests/solve" "$d"
    "$BUILD/tests/test_hlp_grep" "$d"
done
else
echo "-- validate skipped (pass --validate to run solve + test)"
fi
echo "== done"
} 2>&1 | tee "$LOG"
