#!/usr/bin/env python3
"""Run parameter experiments against HLP-grep and save benchmark results.

Usage:
    python3 tests/scripts/run_parameter_experiments.py DATASET_DIR RESULTS_DIR
    python3 tests/scripts/run_parameter_experiments.py DATASET_DIR RESULTS_DIR --quick

Each testcase compares HLP-grep with the selected baseline methods. When a
sibling .hlpg exists, both fresh-build and prebuilt-index query runs are
measured. Parameter runs use one worker. A separate scorer/thread experiment
compares HLP-grep's default and A* scorers over increasing worker counts on a
representative prebuilt index. Results are checked for identical query output.
"""

import argparse
import csv
import datetime
import json
import os
import subprocess
import sys


def method_signatures(data):
    """Return each benchmark method's result signature by query index."""
    return {
        method["method"]: [
            (query["index"], query["k"], query["num_matches"],
             query["results_hash"])
            for query in method["queries"]
        ]
        for method in data["methods"]
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("dataset_dir", help="folder containing generated testcase .txt files")
    parser.add_argument("results_dir", help="folder in which run JSON and CSV results are saved")
    parser.add_argument("--bench-bin", default=None,
                        help="path to bench executable (default: <repo>/build/tests/bench)")
    parser.add_argument("--methods", default="hlp_grep,wfa,hstree,dt_patricia,bed_tree,edlib,parasail",
                        help="comma-separated methods to compare (default: all supported methods)")
    parser.add_argument("--repeats", type=int, default=3,
                        help="runs per testcase/method/mode (default: 3)")
    parser.add_argument("--quick", action="store_true",
                        help="run one repeat; same datasets and scorer/mode coverage")
    args = parser.parse_args()
    if args.repeats < 1:
        parser.error("--repeats must be positive")
    methods = tuple(dict.fromkeys(name.strip() for name in args.methods.split(",")
                                  if name.strip()))
    if "hlp_grep" not in methods:
        parser.error("--methods must include hlp_grep for exact-result cross-checks")

    script_dir = os.path.dirname(os.path.abspath(__file__))
    repo_root = os.path.abspath(os.path.join(script_dir, "..", ".."))
    bench_bin = args.bench_bin or os.path.join(repo_root, "build", "tests", "bench")
    dataset_dir = os.path.abspath(args.dataset_dir)
    if not os.path.isdir(dataset_dir):
        parser.error(f"dataset folder does not exist: {dataset_dir}")
    if not os.path.isfile(bench_bin) or not os.access(bench_bin, os.X_OK):
        parser.error(f"bench executable not found or not executable: {bench_bin}; "
                     "build it with `cmake --build build --target bench` or pass --bench-bin")

    testcase_paths = sorted(
        os.path.join(dataset_dir, name)
        for name in os.listdir(dataset_dir)
        if name.endswith(".txt") and os.path.isfile(os.path.join(dataset_dir, name)))
    if not testcase_paths:
        parser.error(f"no testcase .txt files found in {dataset_dir}")

    stamp = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
    run_dir = os.path.join(os.path.abspath(args.results_dir), f"run-{stamp}")
    os.makedirs(run_dir, exist_ok=False)
    repeats = 1 if args.quick else args.repeats
    scorers = ("default", "a-star")
    summary_rows = []
    run_metadata = []
    baseline_signatures = {}

    def execute(testcase, name, mode, repeat, workers, selected_methods,
                index_path=None, scorer=None):
        worker_tag = f".t{workers}" if mode == "threads" else ""
        scorer_tag = f".{scorer}" if scorer else ""
        out_json = os.path.join(run_dir,
                                f"{name}.{mode}{worker_tag}{scorer_tag}.rep{repeat}.json")
        command = [bench_bin]
        for method_name in selected_methods:
            command.extend(("--method", method_name))
        command.append(testcase)
        if scorer:
            command.extend(("--scorer", scorer))
        command.extend(("--out", out_json))
        if index_path:
            command.extend(("--index", index_path))
        env = os.environ.copy()
        # Solver's automatic worker setting uses omp_get_max_threads()/2.
        # Scale the environment value so the requested worker count is exact.
        env["OMP_NUM_THREADS"] = str(1 if workers == 1 else workers * 2)
        env["OMP_DYNAMIC"] = "FALSE"
        print(f"RUN workers={workers} " + " ".join(command), flush=True)
        subprocess.run(command, check=True, env=env)
        with open(out_json, encoding="utf-8") as source:
            data = json.load(source)
        signatures = method_signatures(data)
        for method in data["methods"]:
            method_name = method["method"]
            for query in method["queries"]:
                summary_rows.append({
                    "testcase": name,
                    "mode": mode,
                    "method": method_name,
                    "scorer": scorer if method_name == "hlp_grep" and scorer else
                              ("default" if method_name == "hlp_grep" else "n/a"),
                    "repeat": repeat,
                    "thread_count": workers,
                    "build_time_ms": method["build_time_ms"],
                    "peak_memory_bytes": method["peak_memory_bytes"],
                    "query_index": query["index"],
                    "k": query["k"],
                    "query_length": query["len"],
                    "query_runtime_ms": query["runtime_ms"],
                    "num_matches": query["num_matches"],
                    "results_hash": query["results_hash"],
                    "json_file": os.path.basename(out_json),
                })
        run_metadata.append({
            "testcase": name,
            "mode": mode,
            "methods": selected_methods,
            "scorer": scorer,
            "repeat": repeat,
            "thread_count": workers,
            "json_file": os.path.basename(out_json),
        })
        return signatures

    for testcase in testcase_paths:
        name = os.path.splitext(os.path.basename(testcase))[0]
        index = os.path.splitext(testcase)[0] + ".hlpg"
        modes = [("build", None)]
        if os.path.isfile(index):
            modes.append(("prebuilt", index))
        else:
            print(f"note: no prebuilt index for {name}; skipping prebuilt mode",
                  file=sys.stderr)

        by_mode_repeat = {}
        for mode, index_path in modes:
            for repeat in range(1, repeats + 1):
                signatures = execute(testcase, name, mode, repeat, 1,
                                     methods, index_path)
                hlp_signature = signatures["hlp_grep"]
                mismatches = [method for method, signature in signatures.items()
                              if signature != hlp_signature]
                if mismatches:
                    raise RuntimeError(
                        f"result mismatch for {name}, {mode}, repeat {repeat}: "
                        f"{', '.join(mismatches)} differ from hlp_grep")
                by_mode_repeat[(mode, repeat)] = hlp_signature
                if mode == "prebuilt":
                    baseline_signatures[(name, repeat)] = hlp_signature

        for mode, _index_path in modes:
            for repeat in range(1, repeats + 1):
                if by_mode_repeat[(mode, repeat)] != by_mode_repeat[("build", repeat)]:
                    raise RuntimeError(f"result mismatch between build/prebuilt modes for {name}")

    # Explicit thread-scaling experiment on a representative indexed workload.
    indexed = {
        os.path.splitext(os.path.basename(path))[0]: path
        for path in testcase_paths
        if os.path.isfile(os.path.splitext(path)[0] + ".hlpg")
    }
    preferred = ("size-10000", "length-5000")
    thread_case = next((candidate for candidate in preferred if candidate in indexed), None)
    if thread_case is None:
        size_cases = sorted((name for name in indexed if name.startswith("size-")),
                            key=lambda name: int(name.split("-")[1]))
        thread_case = size_cases[-1] if size_cases else None
    if thread_case is None:
        thread_case = next((name for name in indexed if name.startswith("length-")), None)
    thread_counts = []
    if thread_case is not None:
        thread_testcase = os.path.join(dataset_dir, thread_case + ".txt")
        thread_index = os.path.splitext(thread_testcase)[0] + ".hlpg"
        max_workers = max(1, (os.cpu_count() or 1) // 2)
        thread_counts = [1]
        workers = 2
        while workers <= max_workers:
            thread_counts.append(workers)
            workers *= 2
        for workers in thread_counts:
            for repeat in range(1, repeats + 1):
                for scorer in scorers:
                    signatures = execute(
                        thread_testcase, thread_case, "threads", repeat,
                        workers, ("hlp_grep",), thread_index, scorer)
                    if signatures["hlp_grep"] != baseline_signatures.get(
                            (thread_case, repeat)):
                        raise RuntimeError(
                            f"thread/scorer result mismatch for {thread_case}, "
                            f"scorer={scorer}, workers={workers}, repeat={repeat}")
    else:
        print("note: no indexed size/length testcase; skipping thread sweep",
              file=sys.stderr)

    csv_path = os.path.join(run_dir, "summary.csv")
    columns = ("testcase", "mode", "method", "scorer", "repeat", "thread_count", "build_time_ms",
               "peak_memory_bytes", "query_index", "k", "query_length",
               "query_runtime_ms", "num_matches", "results_hash", "json_file")
    with open(csv_path, "w", newline="", encoding="utf-8") as output:
        writer = csv.DictWriter(output, fieldnames=columns)
        writer.writeheader()
        writer.writerows(summary_rows)

    manifest_path = os.path.join(run_dir, "RUN.json")
    with open(manifest_path, "w", encoding="utf-8") as output:
        json.dump({
            "dataset_dir": dataset_dir,
            "bench_bin": bench_bin,
            "repeats": repeats,
            "methods": methods,
            "scorers": scorers,
            "single_thread_workers": 1,
            "thread_sweep_testcase": thread_case,
            "thread_sweep_workers": thread_counts if thread_case else [],
            "testcases": [os.path.basename(path) for path in testcase_paths],
            "runs": run_metadata,
        }, output, indent=2)
        output.write("\n")
    print(f"Results saved under {run_dir}")
    print(f"Summary: {csv_path}")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except subprocess.CalledProcessError as error:
        print(f"benchmark command failed with exit code {error.returncode}",
              file=sys.stderr)
        sys.exit(error.returncode or 1)
