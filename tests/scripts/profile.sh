#!/usr/bin/env bash
# CPU/memory profile of HLP-grep construction + query phases.
#
#   tests/scripts/profile.sh [--build-dir DIR] [--reps N] [--max-seqs N]
#                            [--tool all|time|stat|record|callgrind|massif]
#                            [testcase-file-or-dir...]
#
# Builds RelWithDebInfo with frame pointers (readable stacks, realistic -O2
# timings) into --build-dir, runs bench_poa for wall-clock phase timing, then
# runs the selected profiler tools when installed. Missing tools are skipped
# with a note. Reports go to <build-dir>/profile-<timestamp>/.
set -u

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
BUILD="$ROOT/build-prof"
REPS=3
MAX_SEQS=0
TOOL=all
ARGS=""

while [ $# -gt 0 ]; do
    case "$1" in
        --build-dir) BUILD="$2"; shift 2 ;;
        --reps) REPS="$2"; shift 2 ;;
        --max-seqs) MAX_SEQS="$2"; shift 2 ;;
        --tool) TOOL="$2"; shift 2 ;;
        -h|--help)
            sed -n '2,12p' "$0"; exit 0 ;;
        *) ARGS="$ARGS $1"; shift ;;
    esac
done

if [ -z "$ARGS" ]; then
    ARGS=" $ROOT/tests/testcases/manual"
fi
if [ "$MAX_SEQS" != 0 ]; then
    # shellcheck disable=SC2086
    BENCH_FLAGS="--max-seqs $MAX_SEQS --reps $REPS $ARGS"
else
    # shellcheck disable=SC2086
    BENCH_FLAGS="--reps $REPS $ARGS"
fi

STAMP="$(date +%Y%m%d-%H%M%S)"
OUT="$BUILD/profile-$STAMP"
mkdir -p "$OUT"

echo "== configure + build ($BUILD)"
cmake -S "$ROOT" -B "$BUILD" \
    -DCMAKE_BUILD_TYPE=RelWithDebInfo \
    -DCMAKE_CXX_FLAGS="-fno-omit-frame-pointer" > "$OUT/cmake.log" 2>&1
cmake --build "$BUILD" -j --target bench_poa >> "$OUT/cmake.log" 2>&1
echo "-- build log: $OUT/cmake.log"

BIN="$BUILD/tests/bench_poa"

want() { [ "$TOOL" = all ] || [ "$TOOL" = "$1" ]; }

if want time; then
    echo "== wall-clock timing (bench_poa separates build_s per config)"
    # shellcheck disable=SC2086
    { time $BIN $BENCH_FLAGS; } 2>&1 | tee "$OUT/time.log"
fi

if want stat; then
    if command -v perf > /dev/null; then
        echo "== perf stat (cycles, cache misses, page faults)"
        # shellcheck disable=SC2086
        perf stat -d $BIN $BENCH_FLAGS 2>&1 | tee "$OUT/perf-stat.log"
    else
        echo "-- perf not found, skipping stat"
    fi
fi

if want record; then
    if command -v perf > /dev/null; then
        echo "== perf record (sampling call stacks)"
        # shellcheck disable=SC2086
        perf record -g --call-graph dwarf -o "$OUT/perf.data" \
            $BIN $BENCH_FLAGS > "$OUT/perf-bench.log" 2>&1
        perf report --stdio -i "$OUT/perf.data" > "$OUT/perf-report.txt" 2>&1
        echo "-- report: $OUT/perf-report.txt"
        echo "-- hint: perf report -i $OUT/perf.data (interactive)"
    else
        echo "-- perf not found, skipping record"
    fi
fi

if want callgrind; then
    if command -v valgrind > /dev/null; then
        echo "== callgrind (call-path inclusive costs; slow)"
        # shellcheck disable=SC2086
        valgrind --tool=callgrind --callgrind-out-file="$OUT/callgrind.out" \
            $BIN --reps 1 $ARGS > "$OUT/callgrind-bench.log" 2>&1
        if command -v callgrind_annotate > /dev/null; then
            callgrind_annotate "$OUT/callgrind.out" > "$OUT/callgrind.txt" 2>&1
            echo "-- report: $OUT/callgrind.txt"
        else
            echo "-- report: $OUT/callgrind.out (view with kcachegrind)"
        fi
    else
        echo "-- valgrind not found, skipping callgrind"
    fi
fi

if want massif; then
    if command -v valgrind > /dev/null; then
        echo "== massif (heap over time)"
        # shellcheck disable=SC2086
        valgrind --tool=massif --massif-out-file="$OUT/massif.out" \
            $BIN --reps 1 $ARGS > "$OUT/massif-bench.log" 2>&1
        if command -v ms_print > /dev/null; then
            ms_print "$OUT/massif.out" > "$OUT/massif.txt" 2>&1
            echo "-- report: $OUT/massif.txt"
        fi
    else
        echo "-- valgrind not found, skipping massif"
    fi
fi

echo "== done: $OUT"
