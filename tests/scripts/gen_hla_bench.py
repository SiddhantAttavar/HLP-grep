#!/usr/bin/env python3
"""Generate a large benchmark testcase from a raw dataset folder.

Builds a testcase in the format described in tests/testcases/README.md:
a dictionary of 500 (by default) full-length sequences taken from the
FASTA files in a folder under tests/testcases/datasets/raw/ and 500
(by default) held-out sequences from the same source, all queried at a
single edit-distance threshold given on the command line.

    python3 tests/scripts/gen_hla_bench.py --source hla/A_nuc.fasta --k 3
    python3 tests/scripts/gen_hla_bench.py --source sarscov2 --k 5 --seed 7 \
        --dict-size 500 --num-queries 500 --seq-len 350 --offset 100 --expand 20

The source argument is either a FASTA file or a folder (inside the raw
directory or given as a path); folders contribute the records of every
*.fasta file they contain. Sequence selection is the same regardless of
the source: fixed-length windows (--seq-len/--offset), ACGT-only records.

Only the standard library is used. Solution (.sol) files are NOT written
here; use build/tests/solve for that, then build/tests/test_hlp_grep.
"""

import argparse
import glob
import os
import random
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DATASETS_DIR = os.path.normpath(os.path.join(HERE, "..", "testcases", "datasets"))
RAW_DIR = os.path.join(DATASETS_DIR, "raw")
DERIVED_DIR = os.path.join(DATASETS_DIR, "derived")

DNA = "ACGT"

MIN_LEN = 20


def resolve_source(src):
    """Return the FASTA files backing the --source argument.

    A file is used directly; a folder contributes every *.fasta below it
    (first level). The argument may be an absolute/relative path or a name
    relative to the raw directory (e.g. 'hla', 'sarscov2',
    'hla/A_nuc.fasta')."""
    path = src if os.path.isabs(src) else os.path.join(RAW_DIR, src)
    if not os.path.exists(path) and os.path.exists(src):
        path = src
    if os.path.isdir(path):
        found = sorted(
            p for p in glob.glob(os.path.join(path, "*.fasta"))
            if os.path.isfile(p))
        if not found:
            raise ValueError(f"{path}: no *.fasta files")
        return found
    if os.path.isfile(path):
        if path.endswith(".fasta"):
            return [path]
        raise ValueError(f"{path}: expected a *.fasta file or a folder")
    raise FileNotFoundError(path)


def read_fasta(path):
    """Return [(header, sequence)] with whitespace stripped, uppercased."""
    records, header, chunks = [], None, []
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(chunks).upper()))
                header, chunks = line[1:], []
            elif header is None:
                raise ValueError(f"{path}: sequence before header")
            else:
                chunks.append("".join(line.split()))
    if header is not None:
        records.append((header, "".join(chunks).upper()))
    return records


def clip(seq, seq_len, offset):
    """Fixed window [offset, offset + seq_len) of the sequence."""
    start = max(0, min(offset, len(seq) - 1))
    if seq_len is None:
        return seq[start:]
    return seq[start:start + seq_len]


def longest_clean_run(seq):
    """Longest ACGT-only stretch (references may carry N spacers)."""
    best, cur = "", ""
    for char in seq:
        if char in DNA:
            cur += char
            if len(cur) > len(best):
                best = cur
        else:
            cur = ""
    return best


def mutate(rng, seq, edits):
    """Random edits: 70% substitutions, 15% deletions, 15% insertions."""
    seq = list(seq)
    for _ in range(edits):
        if not seq:
            break
        op = rng.random()
        pos = rng.randrange(len(seq))
        if op < 0.7:
            choices = [c for c in DNA if c != seq[pos]]
            seq[pos] = rng.choice(choices)
        elif op < 0.85:
            del seq[pos]
        else:
            seq.insert(pos, rng.choice(DNA))
    return "".join(seq)


def write_testcase(path, dictionary, queries, k):
    with open(path, "w") as out:
        out.write("AGCT\n1 1 0 1\n")
        out.write(f"{len(dictionary)}\n")
        for seq in dictionary:
            out.write(f"{seq}\n")
        out.write(f"{len(queries)}\n")
        for seq in queries:
            out.write(f"{k} {seq}\n")


def main():
    parser = argparse.ArgumentParser(
        description="Generate a benchmark testcase from a raw dataset "
        "folder or FASTA file")
    parser.add_argument("--source", required=True,
                        help="FASTA file or raw-dataset folder (e.g. "
                        "'sarscov2', 'hla/A_nuc.fasta'); folders "
                        "contribute every contained *.fasta")
    parser.add_argument("--k", required=True, type=int,
                        help="edit-distance threshold applied to all queries")
    parser.add_argument("--dict-size", type=int, default=500)
    parser.add_argument("--num-queries", type=int, default=500)
    parser.add_argument("--seq-len", type=int, default=None,
                        help="clip every sequence to this length (approximate: "
                        "clipped pieces may be shorter near the locus end)")
    parser.add_argument("--offset", type=int, default=0,
                        help="start clipping at this sequence position")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--expand", type=int, default=0,
                        help="when the source has too few usable sequences, "
                        "grow the pool with r mutants per existing sequence "
                        "(random edits), e.g. --expand 500 for "
                        "single-reference sources like sarscov2/mtdna")
    parser.add_argument("--mutant-edits", type=int, default=5)
    parser.add_argument("--out-dir", default=None)
    args = parser.parse_args()

    if args.k < 0:
        parser.error("--k must be non-negative")
    if args.dict_size < 1 or args.num_queries < 1:
        parser.error("--dict-size and --num-queries must be positive")
    if args.seq_len is not None and args.seq_len < 1:
        parser.error("--seq-len must be positive")
    if args.offset < 0:
        parser.error("--offset must be non-negative")
    if args.seq_len is not None and args.seq_len < MIN_LEN:
        parser.error(f"--seq-len must be at least {MIN_LEN}")
    if args.expand < 0:
        parser.error("--expand must be non-negative")
    if args.mutant_edits < 0:
        parser.error("--mutant-edits must be non-negative")

    try:
        files = resolve_source(args.source)
    except (FileNotFoundError, ValueError) as exc:
        print(f"bad --source: {exc}\n"
              "  sources live in tests/testcases/datasets/raw/; fetch "
              "new ones with tests/scripts/download_datasets.py",
              file=sys.stderr)
        return 1
    print("source files:")
    for path in files:
        print(f"  {os.path.relpath(path, RAW_DIR)} ({os.path.getsize(path)} B)")

    if len(files) == 1:
        name = os.path.splitext(os.path.basename(files[0]))[0]
    else:
        name = os.path.basename(os.path.normpath(args.source))
    seqs = []
    dropped_non_acgt, dropped_short = 0, 0
    total_records = 0
    for path in files:
        records = read_fasta(path)
        total_records += len(records)
        for header, seq in records:
            if (args.seq_len is not None
                    and len(seq) < args.offset + args.seq_len):
                dropped_short += 1
                continue
            piece = clip(seq, args.seq_len, args.offset)
            if any(c not in DNA for c in piece):
                dropped_non_acgt += 1
                continue
            if len(piece) < MIN_LEN:
                dropped_short += 1
                continue
            seqs.append(piece)
    print(f"{total_records} record(s), usable={len(seqs)} "
          f"(dropped {dropped_non_acgt} non-ACGT, {dropped_short} too short)")

    if not seqs:
        # References with non-DNA spacers (e.g. rCRS): keep the longest
        # clean run of each record, full length only.
        for path in files:
            for header, seq in read_fasta(path):
                run = longest_clean_run(seq)
                if len(run) >= MIN_LEN:
                    seqs.append(run)
        if seqs:
            print(f"fallback: using longest ACGT-only runs "
                  f"({len(seqs)} sequence(s))")

    expanded = 0
    if len(seqs) < args.dict_size + args.num_queries:
        if args.expand > 0:
            rng = random.Random(args.seed)
            need = args.dict_size + args.num_queries - len(seqs)
            base_pool = list(seqs)
            for i in range(need):
                expanded += 1
                seqs.append(mutate(rng, base_pool[i % len(base_pool)],
                                   args.mutant_edits))
            print(f"expanded pool with {expanded} mutants "
                  f"(--mutant-edits {args.mutant_edits})")
        else:
            print(f"only {len(seqs)} usable sequences, need "
                  f"{args.dict_size + args.num_queries} "
                  f"('--expand N' can grow the pool with random mutants)",
                  file=sys.stderr)
            return 1

    rng = random.Random(args.seed)
    rng.shuffle(seqs)
    dictionary = seqs[:args.dict_size]
    queries = seqs[args.dict_size:args.dict_size + args.num_queries]

    out_dir = args.out_dir or os.path.join(DERIVED_DIR, "bench")
    os.makedirs(out_dir, exist_ok=True)
    path = os.path.join(out_dir, f"{name}_k{args.k}.txt")
    write_testcase(path, dictionary, queries, args.k)

    lengths = sorted({len(s) for s in dictionary})
    med = lengths[len(lengths) // 2]
    outliers = sum(1 for s in queries if abs(len(s) - med) > args.k)
    print(f"wrote {path} (dict={len(dictionary)} queries={len(queries)} "
          f"k={args.k} len~{med}, {outliers} query length outliers)")
    max_len = max((len(s) for s in dictionary), default=0)
    print(f"length range in dictionary: {min(lengths)}-{max_len}")

    if outliers > len(queries) / 2:
        print("warning: most queries are farther than k in length from the "
              "typical dictionary sequence; they will match few/no paths",
              file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
