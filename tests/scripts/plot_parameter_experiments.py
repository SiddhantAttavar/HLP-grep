#!/usr/bin/env python3
"""Plot summaries produced by run_parameter_experiments.py.

Usage:
    python3 tests/scripts/plot_parameter_experiments.py RESULTS_DIR
    python3 tests/scripts/plot_parameter_experiments.py RESULTS_DIR --out-dir plots

RESULTS_DIR may be a single timestamped run folder or its parent. PNG and SVG
figures are written to RESULTS_DIR/plots unless --out-dir is supplied.
"""

import argparse
import csv
import os
import statistics
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


SCORER_LABELS = {"default": "Default / DnC", "a-star": "A*"}
SCORER_COLORS = {"default": "#1f77b4", "a-star": "#ff7f0e"}
METHOD_LABELS = {
    "hlp_grep": "HLP-Grep",
    "wfa": "WFA2",
    "hstree": "HSTree",
    "dt_patricia": "DT-Patricia",
    "bed_tree": "BED-tree",
    "edlib": "Edlib",
    "parasail": "Parasail",
}
METHOD_COLORS = {
    "hlp_grep": "#1f77b4",
    "wfa": "#ff7f0e",
    "hstree": "#2ca02c",
    "dt_patricia": "#d62728",
    "bed_tree": "#9467bd",
    "edlib": "#8c564b",
    "parasail": "#e377c2",
}


def load_rows(results_dir):
    paths = sorted(os.path.join(root, name)
                   for root, _dirs, files in os.walk(results_dir)
                   for name in files if name == "summary.csv")
    rows = []
    for path in paths:
        with open(path, newline="", encoding="utf-8") as source:
            for row in csv.DictReader(source):
                row["k"] = int(row["k"])
                row["query_runtime_ms"] = float(row["query_runtime_ms"])
                row["num_matches"] = int(row["num_matches"])
                row["build_time_ms"] = float(row["build_time_ms"])
                row["peak_memory_bytes"] = int(row["peak_memory_bytes"])
                row["query_length"] = int(row["query_length"])
                row["repeat"] = int(row["repeat"])
                row["thread_count"] = int(row.get("thread_count", 1))
                row.setdefault("method", "hlp_grep")
                row.setdefault("scorer", "default")
                rows.append(row)
    return rows, paths


def save(fig, path):
    fig.tight_layout()
    fig.savefig(path + ".png", dpi=180)
    fig.savefig(path + ".svg")
    plt.close(fig)


def grouped_runtime(rows, mode="prebuilt"):
    grouped = {}
    for row in rows:
        if row["mode"] != mode:
            continue
        key = (row["testcase"], row["method"], row["scorer"])
        grouped.setdefault(key, []).append(row["query_runtime_ms"])
    return {key: statistics.median(values) for key, values in grouped.items()}


def line_plot(rows, output, title, xlabel, prefix, order_key, mode="prebuilt"):
    values = grouped_runtime(rows, mode)
    names = sorted({name for name, _method, _scorer in values if name.startswith(prefix)},
                   key=order_key)
    if not names:
        return False
    fig, ax = plt.subplots(figsize=(8.5, 5.2))
    methods = sorted({method for _name, method, _scorer in values
                      if method != "hlp_grep" or "default" in _scorer},
                     key=lambda method: (method != "hlp_grep", method))
    for method in methods:
        scorer = "default" if method == "hlp_grep" else "n/a"
        y = [values.get((name, method, scorer)) for name in names]
        xs = [order_key(name) for name in names]
        ax.plot(xs, y, marker="o", linewidth=1.8,
                label=METHOD_LABELS.get(method, method),
                color=METHOD_COLORS.get(method))
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Median query runtime (ms)")
    ax.set_yscale("log")
    ax.grid(True, which="both", linestyle=":", alpha=0.45)
    ax.legend()
    save(fig, output)
    return True


def plot_threshold(rows, output, mode="prebuilt"):
    grouped = {}
    for row in rows:
        if row["mode"] != mode or row["testcase"] != "threshold-sweep":
            continue
        grouped.setdefault((row["method"], row["scorer"], row["k"]), []).append(
            row["query_runtime_ms"])
    if not grouped:
        return False
    fig, ax = plt.subplots(figsize=(8.5, 5.2))
    ks = sorted({key[2] for key in grouped})
    for method in sorted({key[0] for key in grouped}):
        scorer = "default" if method == "hlp_grep" else "n/a"
        ax.plot(ks,
                [statistics.median(grouped[(method, scorer, k)]) for k in ks],
                marker="o", linewidth=1.8,
                label=METHOD_LABELS.get(method, method),
                color=METHOD_COLORS.get(method))
    ax.set_title("Query latency by edit-distance threshold")
    ax.set_xlabel("Threshold k")
    ax.set_ylabel("Median query runtime (ms)")
    ax.set_yscale("log")
    ax.grid(True, which="both", linestyle=":", alpha=0.45)
    ax.legend()
    save(fig, output)
    return True


def plot_matches(rows, output, mode="prebuilt"):
    grouped = {}
    for row in rows:
        if row["mode"] != mode or not row["testcase"].startswith("matches-"):
            continue
        pct = int(row["testcase"].split("-")[1].removesuffix("pct"))
        grouped.setdefault((row["method"], row["scorer"], pct), []).append(
            row["query_runtime_ms"])
    if not grouped:
        return False
    fig, ax = plt.subplots(figsize=(8.5, 5.2))
    pcts = sorted({key[2] for key in grouped})
    for method in sorted({key[0] for key in grouped}):
        scorer = "default" if method == "hlp_grep" else "n/a"
        ax.plot(pcts,
                [statistics.median(grouped[(method, scorer, pct)]) for pct in pcts],
                marker="o", linewidth=1.8,
                label=METHOD_LABELS.get(method, method),
                color=METHOD_COLORS.get(method))
    ax.set_title("Query latency by target match fraction")
    ax.set_xlabel("Target matching dictionary entries (%)")
    ax.set_ylabel("Query runtime (ms)")
    ax.set_yscale("log")
    ax.grid(True, which="both", linestyle=":", alpha=0.45)
    ax.legend()
    save(fig, output)
    return True


def plot_method_scatter(rows, output, mode="prebuilt"):
    values = grouped_runtime(rows, mode)
    names = sorted({name for name, _method, _scorer in values
                    if (name, "hlp_grep", "default") in values})
    methods = sorted({method for _name, method, _scorer in values
                      if method != "hlp_grep"})
    if not names or not methods:
        return False
    points = {
        method: ([values[(name, "hlp_grep", "default")] for name in names],
                 [values[(name, method, "n/a")] for name in names
                  if (name, method, "n/a") in values])
        for method in methods
    }
    methods = [method for method in methods
               if len(points[method][1]) == len(points[method][0])]
    if not methods:
        return False
    all_values = [value for method in methods for series in points[method]
                  for value in series]
    low = min(all_values)
    high = max(all_values)
    fig, ax = plt.subplots(figsize=(6.4, 6.0))
    for method in methods:
        x, y = points[method]
        ax.scatter(x, y, alpha=0.8, color=METHOD_COLORS.get(method),
                   label=METHOD_LABELS.get(method, method))
    ax.plot([low, high], [low, high], linestyle="--", color="#555555",
            label="Equal latency")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("HLP-Grep median latency (ms)")
    ax.set_ylabel("Baseline median latency (ms)")
    ax.set_title("Baseline method latency vs HLP-Grep")
    ax.grid(True, which="both", linestyle=":", alpha=0.45)
    ax.legend()
    save(fig, output)
    return True


def plot_build_resources(rows, output_dir):
    # Each build measurement is repeated once per query in summary.csv; reduce
    # duplicate rows to one value per testcase/method/repeat before aggregation.
    builds = {}
    memories = {}
    for row in rows:
        if row["mode"] != "build" or not row["testcase"].startswith("size-"):
            continue
        key = (row["testcase"], row["method"], row["repeat"])
        builds[key] = row["build_time_ms"]
        memories[key] = row["peak_memory_bytes"] / (1024 * 1024)
    names = sorted({key[0] for key in builds},
                   key=lambda name: int(name.split("-")[1]))
    if not names:
        return False
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.8))
    methods = sorted({key[1] for key in builds},
                     key=lambda method: (method != "hlp_grep", method))
    for method in methods:
        build_medians = []
        memory_medians = []
        for name in names:
            build_medians.append(statistics.median(
                value for (case, method_name, _repeat), value in builds.items()
                if case == name and method_name == method))
            memory_medians.append(statistics.median(
                value for (case, method_name, _repeat), value in memories.items()
                if case == name and method_name == method))
        sizes = [int(name.split("-")[1]) for name in names]
        axes[0].plot(sizes, build_medians, marker="o",
                     label=METHOD_LABELS.get(method, method),
                     color=METHOD_COLORS.get(method))
        axes[1].plot(sizes, memory_medians, marker="o",
                     label=METHOD_LABELS.get(method, method),
                     color=METHOD_COLORS.get(method))
    axes[0].set_xlabel("Dictionary size")
    axes[0].set_ylabel("Index build time (ms)")
    axes[0].set_title("Index construction")
    axes[1].set_xlabel("Dictionary size")
    axes[1].set_ylabel("Peak RSS (MiB)")
    axes[1].set_title("Peak memory")
    for ax in axes:
        ax.set_xscale("log")
        ax.grid(True, which="both", linestyle=":", alpha=0.45)
        ax.legend()
    save(fig, os.path.join(output_dir, "build-memory-by-size"))
    return True


def plot_thread_scaling(rows, output):
    grouped = {}
    for row in rows:
        if row["mode"] != "threads":
            continue
        grouped.setdefault((row["scorer"], row["thread_count"]), []).append(
            row["query_runtime_ms"])
    if not grouped:
        return False
    workers = sorted({count for _scorer, count in grouped})
    fig, ax = plt.subplots(figsize=(8.5, 5.2))
    for scorer in ("default", "a-star"):
        latency = [statistics.median(grouped[(scorer, count)]) for count in workers]
        baseline = latency[workers.index(1)]
        ax.plot(workers, latency, marker="o", color=SCORER_COLORS[scorer],
                label=f"{SCORER_LABELS[scorer]} latency")
    ax.set_xlabel("Worker threads")
    ax.set_ylabel("Median query latency (ms)")
    ax.set_yscale("log")
    ax.grid(True, which="both", linestyle=":", alpha=0.45)
    ax2 = ax.twinx()
    for scorer in ("default", "a-star"):
        latency = [statistics.median(grouped[(scorer, count)]) for count in workers]
        baseline = latency[workers.index(1)]
        speedup = [baseline / value for value in latency]
        ax2.plot(workers, speedup, marker="s", linestyle="--",
                 color=SCORER_COLORS[scorer],
                 label=f"{SCORER_LABELS[scorer]} speedup")
    ax2.set_ylabel("Speedup", color="#d62728")
    ax2.axhline(1.0, color="#777777", linewidth=0.8, linestyle=":")
    ax.set_title("A* vs default scorer across thread counts")
    lines, labels = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(lines + lines2, labels + labels2, loc="best")
    save(fig, output)
    return True


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("results_dir", help="run folder or parent containing summary.csv files")
    parser.add_argument("--out-dir", default=None,
                        help="plot output folder (default: RESULTS_DIR/plots)")
    parser.add_argument("--mode", choices=("prebuilt", "build"), default="prebuilt",
                        help="use query-only prebuilt results or fresh-build runs (default: prebuilt)")
    args = parser.parse_args()
    results_dir = os.path.abspath(args.results_dir)
    if not os.path.isdir(results_dir):
        parser.error(f"results folder does not exist: {results_dir}")
    rows, summaries = load_rows(results_dir)
    if not rows:
        parser.error(f"no non-empty summary.csv files found under {results_dir}")
    output_dir = os.path.abspath(args.out_dir or os.path.join(results_dir, "plots"))
    os.makedirs(output_dir, exist_ok=True)
    mode = args.mode
    figures = []
    definitions = [
        ("query-latency-by-size", "Dictionary-size sweep",
         "Dictionary size", "size-",
         lambda name: int(name.split("-")[1])),
        ("query-latency-by-mutation", "Mutation-rate sweep",
         "Mutation rate (%)", "mutation-",
         lambda name: float(name.removeprefix("mutation-").removesuffix("pct").replace("p", "."))),
        ("query-latency-by-length", "Sequence-length sweep",
         "Sequence length (bp)", "length-",
         lambda name: int(name.split("-")[1])),
    ]
    for filename, title, xlabel, prefix, key in definitions:
        if line_plot(rows, os.path.join(output_dir, filename), title, xlabel,
                     prefix, key, mode):
            figures.append(filename)
    if plot_threshold(rows, os.path.join(output_dir, "query-latency-by-threshold"), mode):
        figures.append("query-latency-by-threshold")
    if plot_matches(rows, os.path.join(output_dir, "query-latency-by-matches"), mode):
        figures.append("query-latency-by-matches")
    if plot_method_scatter(rows, os.path.join(output_dir, "method-comparison"), mode):
        figures.append("method-comparison")
    if mode == "build" and plot_build_resources(rows, output_dir):
        figures.append("build-memory-by-size")
    if plot_thread_scaling(rows, os.path.join(output_dir, "thread-scaling")):
        figures.append("thread-scaling")
    if not figures:
        parser.error(f"no plot-compatible rows found for mode={mode}")
    with open(os.path.join(output_dir, "PLOTS.txt"), "w", encoding="utf-8") as output:
        output.write(f"Input summaries ({len(summaries)}):\n")
        output.writelines(f"- {path}\n" for path in summaries)
        output.write(f"\nMode: {mode}\nFigures (PNG and SVG):\n")
        output.writelines(f"- {figure}\n" for figure in figures)
    print(f"Wrote {len(figures)} plots (PNG + SVG) to {output_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
