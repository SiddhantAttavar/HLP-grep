#!/usr/bin/env python3
"""Generate synthetic HLP-grep testcase files from a random reference.

Builds one random DNA reference sequence of length --avg-seq-len, then
derives the dictionary and the queries as per-base Bernoulli variations
of it (substitutions at --sub-rate, insertions/deletions splitting
--indel-rate 50/50). Writes a single testcase file in the format
described in tests/testcases/README.md, crossing every query with every
threshold in --ks (dict[0] and query 0 are the exact reference, so k=0
has a guaranteed hit).

    python3 tests/scripts/testcase_gen.py --dict-size 100 --num-queries 5 \\
        --avg-seq-len 200 --ks 0,2,5 --sub-rate 0.02 --indel-rate 0.01

Only the standard library is used. Solution (.sol) files are NOT written
here; use build/tests/solve for that, then build/tests/test_hlp_grep.
"""

import argparse
import os
import random
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_OUT_DIR = os.path.normpath(
    os.path.join(HERE, "..", "testcases", "generated"))

ALPHABET = "AGCT"
DNA = "ACGT"
COST_LINE = "1 1 0 1"


def random_seq(rng, length):
    return "".join(rng.choice(DNA) for _ in range(length))


def vary(rng, seq, sub_rate, indel_rate):
    """Per-base variation: substitute with prob sub_rate, insert or delete
    with prob indel_rate/2 each (inserted base is uniform over ACGT)."""
    out = []
    half_indel = indel_rate / 2
    for base in seq:
        r = rng.random()
        if r < sub_rate:
            out.append(rng.choice([c for c in DNA if c != base]))
        elif r < sub_rate + half_indel:
            out.append(rng.choice(DNA))
            out.append(base)
        elif r < sub_rate + indel_rate:
            pass  # deletion
        else:
            out.append(base)
    return "".join(out)


def write_testcase(path, dictionary, queries):
    with open(path, "w") as out:
        out.write(f"{ALPHABET}\n{COST_LINE}\n")
        out.write(f"{len(dictionary)}\n")
        for seq in dictionary:
            out.write(f"{seq}\n")
        out.write(f"{len(queries)}\n")
        for k, seq in queries:
            out.write(f"{k} {seq}\n")


def default_name(args):
    return (f"gen_d{args.dict_size}_q{args.num_queries}_l{args.avg_seq_len}"
            f"_k{'-'.join(str(k) for k in args.ks)}"
            f"_s{args.sub_rate}_i{args.indel_rate}_seed{args.seed}")


def parse_ks(text):
    try:
        ks = tuple(int(x) for x in text.split(","))
    except ValueError:
        raise argparse.ArgumentTypeError(
            f"ks must be comma-separated integers, got {text!r}")
    if not ks or any(k < 0 for k in ks):
        raise argparse.ArgumentTypeError("ks must be non-empty, non-negative")
    return ks


def main():
    parser = argparse.ArgumentParser(description="Generate synthetic testcases")
    parser.add_argument("--dict-size", type=int, required=True,
                        help="number of dictionary sequences (>= 1)")
    parser.add_argument("--num-queries", type=int, required=True,
                        help="number of query sequences (>= 1)")
    parser.add_argument("--avg-seq-len", type=int, required=True,
                        help="reference sequence length (>= 1)")
    parser.add_argument("--ks", type=parse_ks, required=True,
                        help="comma-separated thresholds, e.g. 0,2,5")
    parser.add_argument("--sub-rate", type=float, required=True,
                        help="per-base substitution probability in [0, 1]")
    parser.add_argument("--indel-rate", type=float, required=True,
                        help="per-base indel probability in [0, 1], "
                        "split 50/50 between insertions and deletions")
    parser.add_argument("--name", default=None,
                        help="output basename (default: derived from params)")
    parser.add_argument("--out-dir", default=DEFAULT_OUT_DIR)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    if args.dict_size < 1:
        parser.error("--dict-size must be >= 1")
    if args.num_queries < 1:
        parser.error("--num-queries must be >= 1")
    if args.avg_seq_len < 1:
        parser.error("--avg-seq-len must be >= 1")
    if not 0.0 <= args.sub_rate <= 1.0:
        parser.error("--sub-rate must be in [0, 1]")
    if not 0.0 <= args.indel_rate <= 1.0:
        parser.error("--indel-rate must be in [0, 1]")

    rng = random.Random(args.seed)
    reference = random_seq(rng, args.avg_seq_len)
    dictionary = [reference] + [vary(rng, reference, args.sub_rate,
                                     args.indel_rate)
                                for _ in range(args.dict_size - 1)]
    variants = [reference] + [vary(rng, reference, args.sub_rate,
                                   args.indel_rate)
                              for _ in range(args.num_queries - 1)]
    queries = [(k, seq) for seq in variants for k in args.ks]

    os.makedirs(args.out_dir, exist_ok=True)
    path = os.path.join(args.out_dir, f"{args.name or default_name(args)}.txt")
    if os.path.exists(path):
        print(f"refusing to overwrite existing {path} "
              "(pass --name to pick another)", file=sys.stderr)
        return 1
    write_testcase(path, dictionary, queries)
    print(f"wrote {path} (dict={len(dictionary)} "
          f"queries={len(queries)} len~{len(reference)})")
    with open(os.path.join(args.out_dir, "MANIFEST.md"), "w") as out:
        out.write("# Generated testcase manifest\n"
                  f"seed={args.seed}\n\n"
                  f"- {os.path.basename(path)} dict={len(dictionary)} "
                  f"queries={len(queries)} len~{len(reference)}\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
