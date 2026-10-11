#!/usr/bin/env python3
"""Plot generated benchmark results, preferring a fresh HLP-grep run.

For matching benchmark JSONs, query timings are grouped by threshold and
summarized across the ten query sequences. Build time is not included because
the JSON's runtime_ms records query time only.
"""

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_BASELINE = ROOT / (
    "tests/testcases/generated/"
    "gen_d100000_q10_l1000_k1-2-3-5-7-10-15-20-40_"
    "s0.01_i0.01_seed42.bench.json"
)
DEFAULT_HLP = ROOT / (
    "tests/testcases/generated/"
    "gen_d100000_q10_l1000_k1-2-3-5-7-10-15-20-40_"
    "s0.01_i0.01_seed42.bench3.json"
)
DEFAULT_OUTPUT = ROOT / "notes/images"

COLORS = {
    "bed_tree": "#9467bd",
    "hstree": "#2ca02c",
    "wfa": "#ff7f0e",
    "dt_patricia": "#d62728",
    "hlp_grep": "#1f77b4",
}
LABELS = {
    "bed_tree": "BED-tree",
    "hstree": "HSTree",
    "wfa": "WFA",
    "dt_patricia": "DT-Patricia",
    "hlp_grep": "HLP-grep (bench3)",
}


def load_json(path):
    with path.open() as stream:
        return json.load(stream)


def index_methods(data):
    return {method["method"]: method for method in data["methods"]}


def query_groups(method):
    grouped = {}
    for query in method["queries"]:
        grouped.setdefault(int(query["k"]), []).append(query)
    return grouped


def validate_matching_runs(baseline, fresh_hlp):
    if baseline["dict_size"] != fresh_hlp["dict_size"]:
        raise ValueError("Benchmark runs have different dictionary sizes")
    if baseline["num_queries"] != fresh_hlp["num_queries"]:
        raise ValueError("Benchmark runs have different query counts")

    old_hlp = index_methods(baseline).get("hlp_grep")
    new_hlp = index_methods(fresh_hlp).get("hlp_grep")
    if old_hlp is None or new_hlp is None:
        raise ValueError("Both benchmark files must contain HLP-grep")

    def signatures(method):
        return [
            (q["index"], q["k"], q["len"], q["num_matches"], q["results_hash"])
            for q in method["queries"]
        ]

    if signatures(old_hlp) != signatures(new_hlp):
        raise ValueError("The fresh HLP-grep run does not match baseline queries/results")


def plot_latency(methods, ks, output, summary):
    latency_ks = [k for k in ks if k <= 20]
    fig = plt.figure(figsize=(9.2, 5.7))
    ax = fig.add_subplot(1, 1, 1)
    for name, method in methods.items():
        groups = query_groups(method)
        samples = [
            np.array([q["runtime_ms"] for q in groups[k]]) for k in latency_ks
        ]
        means = np.array([np.mean(values) for values in samples])
        style = {
            "label": LABELS.get(name, name),
            "color": COLORS.get(name),
            "marker": "o",
            "linewidth": 2.0,
            "markersize": 5,
        }

        if summary == "mean":
            ax.plot(latency_ks, means, **style)
        elif summary == "mean-std":
            stds = np.array([np.std(values, ddof=1) for values in samples])
            ax.errorbar(latency_ks, means, yerr=stds, capsize=2.5, **style)
        else:
            lower = np.array([np.percentile(values, 25) for values in samples])
            medians = np.array([np.percentile(values, 50) for values in samples])
            upper = np.array([np.percentile(values, 75) for values in samples])
            ax.errorbar(
                latency_ks,
                medians,
                yerr=np.vstack((medians - lower, upper - medians)),
                capsize=2.5,
                **style,
            )

    ax.set_yscale("log")
    ax.set_ylim(bottom=0.01)
    ax.set_xticks(latency_ks)
    ax.set_xlim(min(latency_ks) - 0.5, max(latency_ks) + 0.5)
    ax.set_xlabel("Edit-distance threshold k")
    ax.set_ylabel("Query runtime (ms, log scale)")
    ax.set_title("Query latency (k ≤ 20) — 100,000 sequences, length ≈ 1,000")
    ax.grid(True, which="both", linestyle=":", alpha=0.45)
    ax.legend(frameon=True, ncol=2)
    footnotes = {
        "mean": "Mean over 10 queries; build excluded.",
        "mean-std": (
            "Mean ± sample SD over 10 queries; lower bounds may be clipped by "
            "the log scale. Build excluded."
        ),
        "quartiles": (
            "Median markers; whiskers show 25th–75th percentiles over 10 "
            "queries. Build excluded."
        ),
    }
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    fig.text(
        0.12,
        0.015,
        footnotes[summary],
        fontsize=8,
        color="#444444",
    )
    save_figure(fig, output)


def plot_speedup(methods, ks, output):
    groups = {name: query_groups(method) for name, method in methods.items()}
    wfa_means = np.array(
        [np.mean([q["runtime_ms"] for q in groups["wfa"][k]]) for k in ks]
    )

    fig = plt.figure(figsize=(9.2, 5.7))
    ax = fig.add_subplot(1, 1, 1)
    for name, method_groups in groups.items():
        means = np.array(
            [np.mean([q["runtime_ms"] for q in method_groups[k]]) for k in ks]
        )
        ax.plot(
            ks,
            wfa_means / means,
            label=LABELS.get(name, name),
            color=COLORS.get(name),
            marker="o",
            linewidth=2.0,
            markersize=5,
        )

    ax.axhline(1.0, color="#333333", linewidth=1, linestyle="--")
    ax.set_yscale("log")
    ax.set_xticks(ks)
    ax.set_xlabel("Edit-distance threshold k")
    ax.set_ylabel("Speedup vs WFA (WFA mean / method mean)")
    ax.set_title("Relative query latency — 100,000 sequences, length ≈ 1,000")
    ax.grid(True, which="both", linestyle=":", alpha=0.45)
    ax.legend(frameon=True, ncol=2)
    ax.text(
        0.01,
        0.015,
        "Above 1× is faster than WFA. Query time only; build excluded.",
        transform=ax.transAxes,
        fontsize=8,
        color="#444444",
    )
    fig.tight_layout()
    save_figure(fig, output)


def save_figure(fig, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=180)
    fig.savefig(path.with_suffix(".svg"))
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path, default=DEFAULT_BASELINE)
    parser.add_argument("--hlp", type=Path, default=DEFAULT_HLP)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument(
        "--latency-summary",
        choices=("mean", "mean-std", "quartiles"),
        default="quartiles",
        help=(
            "absolute latency plot summary: mean, mean with sample standard "
            "deviation, or median with 25th/75th percentile whiskers "
            "(default: quartiles)"
        ),
    )
    args = parser.parse_args()

    baseline = load_json(args.baseline)
    fresh_hlp = load_json(args.hlp)
    validate_matching_runs(baseline, fresh_hlp)

    methods = index_methods(baseline)
    methods["hlp_grep"] = index_methods(fresh_hlp)["hlp_grep"]
    required = ["bed_tree", "hstree", "wfa", "dt_patricia", "hlp_grep"]
    missing = [name for name in required if name not in methods]
    if missing:
        raise ValueError(f"Missing methods in benchmark data: {', '.join(missing)}")
    methods = {name: methods[name] for name in required}

    groups = {name: query_groups(method) for name, method in methods.items()}
    ks = sorted(groups["hlp_grep"])
    for name in required:
        if sorted(groups[name]) != ks or any(len(groups[name][k]) != 10 for k in ks):
            raise ValueError(f"Unexpected threshold/query grouping for {name}")

    plot_latency(
        methods,
        ks,
        args.output_dir / "query_latency_d100000_l1000.png",
        args.latency_summary,
    )
    plot_speedup(
        methods,
        ks,
        args.output_dir / "speedup_vs_wfa_d100000_l1000.png",
    )
    print(f"Plots written to {args.output_dir}")
    print("HLP-grep timings sourced from:", args.hlp)


if __name__ == "__main__":
    main()
