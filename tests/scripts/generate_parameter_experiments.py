#!/usr/bin/env python3
"""Generate testcase files for the parameter/scorer experiment matrix.

Usage:
    python3 tests/scripts/generate_parameter_experiments.py OUTPUT_DIR
    python3 tests/scripts/generate_parameter_experiments.py OUTPUT_DIR --quick
    python3 tests/scripts/generate_parameter_experiments.py OUTPUT_DIR --overwrite

The full profile follows the planned 1,000 bp synthetic experiments. ``--quick``
reduces sizes for smoke testing. Existing HLA FASTA records under
tests/testcases/datasets/raw/hla/ are also used to create a real-data validation
case when available. Thread and scorer comparisons reuse the generated inputs.
"""

import argparse
import glob
import os
import random
import subprocess
import sys

from generate_mitomap_mtdna_test import (
    DEFAULT_CACHE_DIR as MITOMAP_CACHE_DIR,
    load_accessions as load_mitomap_accessions,
    select_records as select_mitomap_records,
)

DNA = "ACGT"
SIZE_LEVELS = (100, 500, 1_000, 5_000, 10_000, 50_000, 100_000)
MUTATION_RATES = (0.001, 0.005, 0.01, 0.02, 0.05)
THRESHOLDS = (1, 3, 5, 10, 15, 20)


def random_sequence(rng, length):
    return "".join(rng.choices(DNA, k=length))


def mutate(rng, sequence, rate):
    """Apply substitutions at an exact rounded per-base mutation rate."""
    result = list(sequence)
    edits = min(len(result), max(1, round(len(result) * rate))) if rate > 0 else 0
    for position in rng.sample(range(len(result)), edits):
        old = result[position]
        result[position] = rng.choice([base for base in DNA if base != old])
    return "".join(result)


def write_case(path, dictionary, query_pairs):
    with open(path, "w", encoding="ascii") as output:
        output.write("AGCT\n1 1 0 1\n")
        output.write(f"{len(dictionary)}\n")
        output.writelines(sequence + "\n" for sequence in dictionary)
        output.write(f"{len(query_pairs)}\n")
        output.writelines(f"{k} {sequence}\n" for k, sequence in query_pairs)


def fresh_queries(rng, base, count, rate):
    return [mutate(rng, base, rate) for _ in range(count)]


def add_case(out_dir, manifest, name, dictionary, queries, details,
             overwrite=False):
    path = os.path.join(out_dir, name + ".txt")
    if os.path.exists(path) and not overwrite:
        print(f"skip existing {path}")
    else:
        write_case(path, dictionary, queries)
        print(f"wrote {path} (dict={len(dictionary)} queries={len(queries)})")
    manifest.append((name, details))


def read_fasta(path):
    records = []
    chunks = []
    with open(path, encoding="ascii", errors="replace") as source:
        for line in source:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if chunks:
                    records.append("".join(chunks).upper())
                chunks = []
            else:
                chunks.append(line)
    if chunks:
        records.append("".join(chunks).upper())
    return [sequence for sequence in records
            if sequence and all(base in DNA for base in sequence)]


def add_match_cases(out_dir, manifest, rng, base, n, length, query_count,
                    overwrite=False):
    """Create five one-query cases with controlled exact-match cardinality."""
    for fraction in (0.0, 0.01, 0.1, 0.5, 1.0):
        target = round(n * fraction)
        dictionary = [random_sequence(rng, length) for _ in range(n)]
        query = random_sequence(rng, length)
        if target:
            # One query per testcase lets each requested cardinality be exact
            # without assigning more than n dictionary slots across queries.
            dictionary[:target] = [query] * target
        queries = [(0, query)]
        name = f"matches-{int(fraction * 100):03d}pct"
        add_case(out_dir, manifest, name, dictionary, queries,
                 f"n={n}; L={length}; target_matches={target}; k=0",
                 overwrite)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output_dir", help="folder in which testcase files are written")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--quick", action="store_true",
                        help="small profile for smoke tests; not the full planned matrix")
    parser.add_argument("--overwrite", action="store_true",
                        help="replace existing testcase, index, and manifest files")
    parser.add_argument("--length", type=int, default=1_000,
                        help="synthetic sequence length (default: 1000)")
    parser.add_argument("--queries", type=int, default=100,
                        help="query count per case (default: 100)")
    parser.add_argument("--save-bin", default=None,
                        help="path to save executable (default: <repo>/build/tests/save)")
    args = parser.parse_args()
    if args.length < 20 or args.queries < 1:
        parser.error("--length must be at least 20 and --queries must be positive")

    script_dir = os.path.dirname(os.path.abspath(__file__))
    repo_root = os.path.abspath(os.path.join(script_dir, "..", ".."))
    save_bin = args.save_bin or os.path.join(repo_root, "build", "tests", "save")
    if not os.path.isfile(save_bin) or not os.access(save_bin, os.X_OK):
        parser.error(f"save executable not found or not executable: {save_bin}; "
                     "build it with `cmake --build build --target save` or pass --save-bin")

    sizes = (10, 25, 50) if args.quick else SIZE_LEVELS
    max_size = max(sizes)
    query_count = min(args.queries, 5) if args.quick else args.queries
    length = min(args.length, 100) if args.quick else args.length
    rng = random.Random(args.seed)
    os.makedirs(args.output_dir, exist_ok=True)
    manifest = []

    # Nested subsets make dictionary-size comparisons directly comparable.
    size_base = random_sequence(rng, length)
    size_master = [size_base] + [mutate(rng, size_base, 0.01)
                                 for _ in range(max_size - 1)]
    size_queries = [(10, sequence) for sequence in
                    fresh_queries(rng, size_base, query_count, 0.01)]
    for size in sizes:
        add_case(args.output_dir, manifest, f"size-{size}",
                 size_master[:size], size_queries,
                 f"n={size}; L={length}; mutation=1%; k=10", args.overwrite)

    # Mutation-rate sweep holds n, L, k, and query count constant.
    mutation_n = 1_000 if not args.quick else 50
    mutation_base = random_sequence(rng, length)
    for rate in MUTATION_RATES:
        dictionary = [mutation_base] + [mutate(rng, mutation_base, rate)
                                        for _ in range(mutation_n - 1)]
        queries = [(10, sequence) for sequence in
                   fresh_queries(rng, mutation_base, query_count, rate)]
        label = f"{rate * 100:g}pct".replace(".", "p")
        add_case(args.output_dir, manifest, f"mutation-{label}", dictionary,
                 queries, f"n={mutation_n}; L={length}; mutation={rate:.1%}; k=10",
                 args.overwrite)

    # Threshold sweep uses generated synthetic data only (never real-validation
    # records), with the same dictionary and queries at every threshold.
    threshold_n = 100_000 if not args.quick else 50
    threshold_base = random_sequence(rng, length)
    threshold_dict = [threshold_base] + [mutate(rng, threshold_base, 0.01)
                                         for _ in range(threshold_n - 1)]
    threshold_queries = fresh_queries(rng, threshold_base, query_count, 0.01)
    threshold_pairs = [(k, sequence) for k in THRESHOLDS
                       for sequence in threshold_queries]
    add_case(args.output_dir, manifest, "threshold-sweep", threshold_dict,
             threshold_pairs,
             f"n={threshold_n}; L={length}; mutation=1%; k={THRESHOLDS}",
             args.overwrite)

    # Sequence-length sweep holds n, mutation rate, k, and query count fixed.
    length_levels = (40, 75, 100) if args.quick else (100, 500, 1_000, 5_000)
    length_n = 50 if args.quick else 500
    for seq_length in length_levels:
        base = random_sequence(rng, seq_length)
        dictionary = [base] + [mutate(rng, base, 0.01)
                               for _ in range(length_n - 1)]
        queries = [(10, sequence) for sequence in
                   fresh_queries(rng, base, query_count, 0.01)]
        add_case(args.output_dir, manifest, f"length-{seq_length}",
                 dictionary, queries,
                 f"n={length_n}; L={seq_length}; mutation=1%; k=10",
                 args.overwrite)

    # 2x2x2 interaction matrix.
    combined_n = (25, 50) if args.quick else (500, 10_000)
    for size_name, n in zip(("small", "large"), combined_n):
        for mutation_name, rate in (("low", 0.001), ("high", 0.05)):
            base = random_sequence(rng, length)
            dictionary = [base] + [mutate(rng, base, rate) for _ in range(n - 1)]
            for threshold_name, k in (("low", 1), ("high", 20)):
                queries = [(k, sequence) for sequence in
                           fresh_queries(rng, base, query_count, rate)]
                name = f"combined-{size_name}-{mutation_name}-{threshold_name}"
                add_case(args.output_dir, manifest, name, dictionary, queries,
                         f"n={n}; L={length}; mutation={rate:.1%}; k={k}",
                         args.overwrite)

    # Match-density cases keep n and L fixed; each query has the target number
    # of exact dictionary copies. Distinct random query sequences prevent
    # results from one query being altered while constructing the next.
    match_n = 100_000 if not args.quick else 50
    match_base = random_sequence(rng, length)
    add_match_cases(args.output_dir, manifest, rng, match_base, match_n,
                    length, query_count, args.overwrite)

    # Real HLA validation is optional on installations without the raw FASTA.
    hla_files = sorted(glob.glob(os.path.join(
        repo_root, "tests", "testcases", "datasets", "raw", "hla", "*.fasta")))
    hla_records = []
    for fasta in hla_files:
        hla_records.extend(read_fasta(fasta))
    if len(hla_records) > 1:
        rng.shuffle(hla_records)
        dictionary = hla_records[:min(500, len(hla_records) - 1)]
        heldout = hla_records[len(dictionary):len(dictionary) + query_count]
        if heldout:
            queries = [(k, sequence) for k in (0, 1, 5, 10, 20)
                       for sequence in heldout]
            add_case(args.output_dir, manifest, "real-validation-hla",
                     dictionary, queries,
                     f"source=raw/hla; n={len(dictionary)}; q={len(heldout)}; k=0,1,5,10,20",
                     args.overwrite)
    else:
        manifest.append(("real-validation-hla", "NOT GENERATED: add HLA FASTA under repository raw/hla/"))
        print("note: skipped real HLA validation (no usable raw HLA FASTA found)",
              file=sys.stderr)

    # Large real-data workload sampled from the locally cached MITOMAP set.
    mito_name = (f"mitomap_mtdna_d50000_q100_k1-3-5-10-15-20_"
                 f"seed{args.seed}")
    mito_accession_file = os.path.join(
        MITOMAP_CACHE_DIR, "mitomap_full_length_accessions.txt")
    mito_sequence_file = os.path.join(
        MITOMAP_CACHE_DIR, "full_length_sequences.fasta")
    if args.quick:
        manifest.append((mito_name, "NOT GENERATED in --quick mode"))
    elif not (os.path.isfile(mito_accession_file)
              and os.path.isfile(mito_sequence_file)):
        manifest.append((mito_name,
                         "NOT GENERATED: run tests/scripts/download_mitomap_mtdna.py first"))
        print("note: skipped 50k MITOMAP mtDNA experiment; full MITOMAP sequence "
              "cache is missing", file=sys.stderr)
    else:
        try:
            mito_accessions, _accession_path = load_mitomap_accessions(
                MITOMAP_CACHE_DIR)
            mito_records, _fetched = select_mitomap_records(
                mito_accessions, 50_100, random.Random(args.seed),
                MITOMAP_CACHE_DIR)
        except (OSError, RuntimeError) as error:
            print(f"cannot prepare MITOMAP mtDNA experiment: {error}",
                  file=sys.stderr)
            return 1
        mito_dict = [sequence for _accession, sequence in mito_records[:50_000]]
        mito_heldout = [sequence for _accession, sequence in mito_records[50_000:]]
        mito_queries = [(k, sequence) for sequence in mito_heldout
                        for k in (1, 3, 5, 10, 15, 20)]
        add_case(args.output_dir, manifest, mito_name, mito_dict, mito_queries,
                 "source=MITOMAP full-length human mtDNA; n=50000; "
                 "q=100; k=1,3,5,10,15,20; natural sequence lengths",
                 args.overwrite)

    # Every generated testcase receives its adjacent .hlpg index through save.
    # The real-data manifest entry is skipped when its FASTA source is absent.
    for name, _details in manifest:
        testcase = os.path.join(args.output_dir, name + ".txt")
        if not os.path.isfile(testcase):
            continue
        index = os.path.splitext(testcase)[0] + ".hlpg"
        if os.path.exists(index) and not args.overwrite:
            print(f"skip existing {index}")
            continue
        try:
            save_env = os.environ.copy()
            save_env["OMP_NUM_THREADS"] = "1"
            save_env["OMP_DYNAMIC"] = "FALSE"
            subprocess.run([save_bin, testcase], check=True, env=save_env)
        except subprocess.CalledProcessError as error:
            print(f"failed to save index for {testcase} (exit {error.returncode})",
                  file=sys.stderr)
            return error.returncode or 1

    # Keep the human-readable summary separate from testcase generation.
    summary_path = os.path.join(args.output_dir, "MANIFEST.md")
    if os.path.exists(summary_path) and not args.overwrite:
        print(f"skip existing {summary_path}")
        print(f"datasets are in {args.output_dir}")
        return 0
    with open(summary_path, "w", encoding="utf-8") as output:
        output.write("# Parameter experiment datasets\n\n")
        output.write(f"Seed: `{args.seed}`. Sequence length: `{length}`. "
                     f"Queries per workload: `{query_count}`.\n\n")
        output.write("Algorithm comparisons use the same inputs; the direct A* vs default scorer comparison is isolated to the thread-count experiment.\n\n")
        output.write("Each generated `.txt` testcase has a sibling `.hlpg` index created by `save`.\n\n")
        output.write("Match-density cases contain one query each so the target cardinality is exact.\n\n")
        output.write("| Dataset | Parameters |\n|---|---|\n")
        for name, details in manifest:
            if isinstance(name, str) and name not in ("", "# Parameter experiment datasets"):
                output.write(f"| `{name}` | {details} |\n")
    print(f"wrote {summary_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
