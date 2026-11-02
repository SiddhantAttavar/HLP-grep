#!/usr/bin/env python3
"""Run HLP-grep parameter/scorer experiments and save raw and summary results.

Usage:
    python3 tests/scripts/run_parameter_experiments.py DATASET_DIR RESULTS_DIR
    python3 tests/scripts/run_parameter_experiments.py DATASET_DIR RESULTS_DIR --quick

Each testcase is benchmarked with both the default (shared-prefix/DnC) and A*
scorers. When a sibling .hlpg exists, both fresh-build and prebuilt-index query
runs are measured. Raw bench JSON files and a flattened CSV are saved beneath
RESULTS_DIR. Results from both scorers are checked for identical query output.
"""

import argparse
import csv
import datetime
import json
import os
import subprocess
import sys


def signatures(data):
    """Return query result signatures by index for cross-scorer validation."""
    method = next(item for item in data["methods"]
                  if item["method"] == "hlp_grep")
    return [(query["index"], query["k"], query["num_matches"],
             query["results_hash"]) for query in method["queries"]]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("dataset_dir", help="folder containing generated testcase .txt files")
    parser.add_argument("results_dir", help="folder in which run JSON and CSV results are saved")
    parser.add_argument("--bench-bin", default=None,
                        help="path to bench executable (default: <repo>/build/tests/bench)")
    parser.add_argument("--repeats", type=int, default=3,
                        help="runs per testcase/scorer/mode (default: 3)")
    parser.add_argument("--quick", action="store_true",
                        help="run one repeat; same datasets and scorer/mode coverage")
    args = parser.parse_args()
    if args.repeats < 1:
        parser.error("--repeats must be positive")

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

    for testcase in testcase_paths:
        name = os.path.splitext(os.path.basename(testcase))[0]
        index = os.path.splitext(testcase)[0] + ".hlpg"
        modes = [("build", None)]
        if os.path.isfile(index):
            modes.append(("prebuilt", index))
        else:
            print(f"note: no prebuilt index for {name}; skipping prebuilt mode",
                  file=sys.stderr)

        by_scorer_mode = {}
        for mode, index_path in modes:
            for scorer in scorers:
                for repeat in range(1, repeats + 1):
                    out_json = os.path.join(
                        run_dir, f"{name}.{mode}.{scorer}.rep{repeat}.json")
                    command = [bench_bin, "--method", "hlp_grep", testcase,
                               "--scorer", scorer, "--out", out_json]
                    if index_path:
                        command.extend(("--index", index_path))
                    print("RUN", " ".join(command), flush=True)
                    subprocess.run(command, check=True)
                    with open(out_json, encoding="utf-8") as source:
                        data = json.load(source)
                    by_scorer_mode[(mode, repeat, scorer)] = signatures(data)
                    method = next(item for item in data["methods"]
                                  if item["method"] == "hlp_grep")
                    for query in method["queries"]:
                        summary_rows.append({
                            "testcase": name,
                            "mode": mode,
                            "scorer": scorer,
                            "repeat": repeat,
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
                        "scorer": scorer,
                        "repeat": repeat,
                        "json_file": os.path.basename(out_json),
                    })

        for mode, _index_path in modes:
            for repeat in range(1, repeats + 1):
                if (by_scorer_mode[(mode, repeat, "default")] !=
                        by_scorer_mode[(mode, repeat, "a-star")]):
                    raise RuntimeError(
                        f"scorer result mismatch for {name}, {mode}, repeat {repeat}")

    csv_path = os.path.join(run_dir, "summary.csv")
    columns = ("testcase", "mode", "scorer", "repeat", "build_time_ms",
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
            "scorers": scorers,
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
