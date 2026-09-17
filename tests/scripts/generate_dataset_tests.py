#!/usr/bin/env python3
"""Generate HLP-grep testcase files from real-world DNA sources.

Reads FASTA records from tests/testcases/datasets/raw/ (or --from-fasta),
builds highly similar dictionaries (per-locus alleles, genome windows,
mutant cohorts) and writes testcase files in the format described in
tests/testcases/README.md to tests/testcases/datasets/derived/<dataset>/.

    python3 tests/scripts/generate_dataset_tests.py --list
    python3 tests/scripts/generate_dataset_tests.py --dataset synthetic
    python3 tests/scripts/generate_dataset_tests.py --dataset hla --seed 7

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

# Per-dataset defaults: (window, dict_size, num_queries, mutant_edits, ks).
# Kept small enough for the O(n*q*L^2) naive reference solver.
PROFILES = {
    "hla": {"window": 0, "dict_size": 10, "num_queries": 5,
            "mutant_edits": 4, "ks": (0, 2, 5)},
    "mhc": {"window": 1000, "dict_size": 8, "num_queries": 4,
            "mutant_edits": 12, "ks": (0, 5, 20)},
    "sarscov2": {"window": 1500, "dict_size": 6, "num_queries": 4,
                 "mutant_edits": 8, "ks": (0, 5, 10)},
    "mtdna": {"window": 1500, "dict_size": 6, "num_queries": 4,
              "mutant_edits": 6, "ks": (0, 5, 10)},
    "markers": {"window": 0, "dict_size": 10, "num_queries": 5,
                "mutant_edits": 3, "ks": (0, 1, 3)},
    "theseus": {"window": 2000, "dict_size": 6, "num_queries": 4,
                "mutant_edits": 15, "ks": (0, 10, 30)},
    "synthetic": {"window": 0, "dict_size": 12, "num_queries": 5,
                  "mutant_edits": 6, "ks": (0, 2, 6)},
    "uniform": {"window": 0, "dict_size": 12, "num_queries": 3,
                "mutant_edits": 0, "ks": (0, 1, 2)},
}

# Raw files backing each dataset in curated mode (see download_datasets.py).
# For hla, every locus FASTA present in raw/hla/ is used, so --full
# downloads are picked up automatically.
RAW_SOURCES = {
    "hla": ["hla/DMA_nuc.fasta", "hla/DPA2_nuc.fasta", "hla/DPB2_nuc.fasta",
            "hla/DQA2_nuc.fasta", "hla/DQB2_nuc.fasta"],
    "mhc": ["mhc/OK649231_w0_10kb.fasta", "mhc/OK649231_w2500k_10kb.fasta"],
    "sarscov2": ["sarscov2/NC_045512.fasta"],
    "mtdna": ["mtdna/NC_012920.fasta"],
}


def resolve_sources(name):
    """Raw files backing a dataset: every locus FASTA present in raw/hla/
    for hla (so --full downloads are picked up), else the curated list."""
    if name == "hla":
        found = sorted(glob.glob(os.path.join(RAW_DIR, "hla", "*.fasta")))
        if found:
            return [os.path.relpath(p, RAW_DIR) for p in found]
    return RAW_SOURCES.get(name, [])


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


def clean(seq, min_len=20):
    """Keep ACGT-only sequences; return None otherwise."""
    if len(seq) < min_len or any(c not in DNA for c in seq):
        return None
    return seq


def mutate(rng, seq, edits):
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


def random_seq(rng, length):
    return "".join(rng.choice(DNA) for _ in range(length))


def tile(seq, window):
    """Non-overlapping windows; whole sequence if window <= 0."""
    if window <= 0 or len(seq) <= window:
        return [seq]
    return [seq[i:i + window] for i in range(0, len(seq) - window + 1, window)]


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


def write_testcase(path, dictionary, queries):
    with open(path, "w") as out:
        out.write("AGCT\n1 1 1 1\n1 1 1 1\n"
                  "0 1 1 1\n1 0 1 1\n1 1 0 1\n1 1 1 0\n")
        out.write(f"{len(dictionary)}\n")
        for seq in dictionary:
            out.write(f"{seq}\n")
        out.write(f"{len(queries)}\n")
        for k, seq in queries:
            out.write(f"{k} {seq}\n")


def build_queries(rng, dictionary, heldout, num_queries, mutant_edits, ks):
    """Exact member (k=0) + held-out/mutant sequences across thresholds."""
    queries = [(0, dictionary[0])]
    pool = ([("heldout", s) for s in heldout]
            + [("mutant", mutate(rng, rng.choice(dictionary), mutant_edits))
               for _ in range(num_queries)])
    i = 0
    while len(queries) < num_queries and pool:
        kind, seq = pool[i % len(pool)]
        k = ks[min(1 + (i // 2), len(ks) - 1)] if kind == "heldout" else ks[-1]
        if kind == "heldout" and len(ks) > 1:
            k = ks[1 + (i % (len(ks) - 1))]
        queries.append((k, seq))
        i += 1
        if i > len(pool) * len(ks):
            break
    return queries[:num_queries]


def from_split(name, records, out_dir, seed, window, split, expand,
               mutant_edits, ks, min_len=20, max_queries=None):
    """Full-length, use-all-records testcase: every usable record is kept
    (no subsampling, no windowing unless --window is set); single-molecule
    sources are expanded into a full-length mutant cohort of size `expand`.
    The shuffled records are divided into (split, 1 - split) dict/query,
    and each held-out query is tested at every threshold in ks."""
    rng = random.Random(seed)
    seqs = []
    for header, seq in records:
        for piece in tile(seq, window):
            cleaned = clean(piece, min_len)
            if cleaned is not None:
                seqs.append(cleaned)
    if not seqs:
        # e.g. rCRS carries an N spacer: keep the longest clean run.
        for header, seq in records:
            run = longest_clean_run(seq)
            if len(run) >= min_len:
                seqs.append(run)
    if len(seqs) == 1 and expand > 1:
        seqs += [mutate(rng, seqs[0], mutant_edits) for _ in range(expand - 1)]
    if len(seqs) < 2:
        raise ValueError(f"{name}: only {len(seqs)} usable sequence(s)")
    rng.shuffle(seqs)
    nq = min(max(1, int(round(len(seqs) * (1 - split)))), len(seqs) - 1)
    dictionary, heldout = seqs[:len(seqs) - nq], seqs[len(seqs) - nq:]
    if max_queries is not None:
        heldout = heldout[:max_queries]
    queries = [(k, h) for h in heldout for k in ks]
    os.makedirs(out_dir, exist_ok=True)
    path = os.path.join(out_dir, f"{name}.txt")
    write_testcase(path, dictionary, queries)
    return path, len(dictionary), len(queries), len(dictionary[0])


def from_records(name, records, out_dir, seed, window, dict_size,
                 num_queries, mutant_edits, ks, min_len=20):
    rng = random.Random(seed)
    seqs = []
    for header, seq in records:
        for piece in tile(seq, window):
            cleaned = clean(piece, min_len)
            if cleaned is not None:
                seqs.append(cleaned)
    # De-duplicate, keep order; drop sequences seen verbatim twice is NOT
    # wanted (duplicates are a valid case) — shuffle then split instead.
    rng.shuffle(seqs)
    if len(seqs) < 2:
        raise ValueError(f"{name}: only {len(seqs)} usable sequence(s)")
    dictionary = seqs[:dict_size]
    heldout = seqs[dict_size:dict_size + max(num_queries, 1)]
    queries = build_queries(rng, dictionary, heldout or dictionary[1:],
                            num_queries, mutant_edits, ks)
    os.makedirs(out_dir, exist_ok=True)
    path = os.path.join(out_dir, f"{name}.txt")
    write_testcase(path, dictionary, queries)
    return path, len(dictionary), len(queries), len(dictionary[0])


def split_dataset(name, prof, out_dir, args):
    """Split-mode driver: full-length records, no subsampling; see from_split."""
    made = []
    if name == "synthetic":
        rng = random.Random(args.seed)
        base = random_seq(rng, args.synthetic_length)
        records = [("base", base)] + [
            (f"mut{i}", mutate(rng, base, prof["mutant_edits"]))
            for i in range(args.expand - 1)]
        made.append(from_split("synthetic_mut", records, out_dir, args.seed,
                               0, args.split, 0, prof["mutant_edits"],
                               prof["ks"], max_queries=args.max_queries))
    elif name == "uniform":
        rng = random.Random(args.seed)
        records = [(f"rand{i}", random_seq(rng, 20))
                   for i in range(args.expand)]
        made.append(from_split("uniform_ctrl", records, out_dir, args.seed,
                               0, args.split, 0, 0, prof["ks"],
                               max_queries=args.max_queries))
    else:
        sources = ([args.from_fasta] if args.from_fasta
                   else [os.path.join(RAW_DIR, s)
                         for s in resolve_sources(name)])
        if args.from_fasta is None and name in RAW_SOURCES:
            missing = [s for s in sources if not os.path.exists(s)]
            if missing:
                print(f"skip {name}: missing raw file(s): "
                      f"{', '.join(missing)}\n  run "
                      "tests/scripts/download_datasets.py first "
                      "or pass --from-fasta")
                return made
        if args.from_fasta is None and name not in RAW_SOURCES:
            print(f"skip {name}: needs --from-fasta "
                  "(see tests/testcases/README.md)")
            return made
        for src in sources:
            stem = os.path.splitext(os.path.basename(src))[0]
            records = read_fasta(src)
            made.append(from_split(f"{name}_{stem}", records, out_dir,
                                   args.seed, prof["window"], args.split,
                                   args.expand, prof["mutant_edits"],
                                   prof["ks"],
                                   max_queries=args.max_queries))
    return made


def main():
    parser = argparse.ArgumentParser(description="Generate dataset testcases")
    parser.add_argument("--dataset", action="append", default=[],
                        choices=sorted(PROFILES),
                        help="repeatable; default: all with raw present + synthetic/uniform")
    parser.add_argument("--from-fasta", default=None,
                        help="override input FASTA (for markers/theseus/manual)")
    parser.add_argument("--out-dir", default=None)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--dict-size", type=int, default=None)
    parser.add_argument("--num-queries", type=int, default=None)
    parser.add_argument("--window", type=int, default=None)
    parser.add_argument("--mutant-edits", type=int, default=None)
    parser.add_argument("--ks", default=None,
                        help="comma-separated thresholds, e.g. 0,2,5")
    parser.add_argument("--synthetic-length", type=int, default=1000)
    parser.add_argument("--split", type=float, default=None,
                        help="dict fraction, e.g. 0.9: keep every record at "
                        "full length, split dict/query, cross all ks; "
                        "single-molecule sources expand to --expand mutants")
    parser.add_argument("--expand", type=int, default=20,
                        help="cohort size for single-molecule sources in "
                        "--split mode")
    parser.add_argument("--max-queries", type=int, default=None,
                        help="cap held-out queries per file in --split mode "
                        "(dict keeps all records)")
    parser.add_argument("--list", action="store_true",
                        help="list dataset profiles and exit")
    args = parser.parse_args()

    if args.list:
        for name in sorted(PROFILES):
            srcs = RAW_SOURCES.get(name, ["--from-fasta / synthetic"])
            print(f"{name}: {PROFILES[name]} <- {', '.join(srcs)}")
        return 0

    selected = args.dataset or ["hla", "mhc", "sarscov2", "mtdna",
                                "synthetic", "uniform"]
    ks = tuple(int(x) for x in args.ks.split(",")) if args.ks else None
    made = []
    for name in selected:
        prof = dict(PROFILES[name])
        if args.dict_size is not None:
            prof["dict_size"] = args.dict_size
        if args.num_queries is not None:
            prof["num_queries"] = args.num_queries
        if args.window is not None:
            prof["window"] = args.window
        if args.mutant_edits is not None:
            prof["mutant_edits"] = args.mutant_edits
        if ks is not None:
            prof["ks"] = ks
        out_dir = args.out_dir or os.path.join(DERIVED_DIR, name)

        if args.split is not None:
            made.extend(split_dataset(name, prof, out_dir, args))
            continue

        if name == "synthetic":
            rng = random.Random(args.seed)
            base = random_seq(rng, args.synthetic_length)
            records = [(f"mut{i}",
                        mutate(rng, base, rng.randint(0, 2 * prof["mutant_edits"])))
                       for i in range(prof["dict_size"] + prof["num_queries"])]
            made.append(from_records("synthetic_mut", records, out_dir,
                                     args.seed, 0, prof["dict_size"],
                                     prof["num_queries"], prof["mutant_edits"],
                                     prof["ks"]))
        elif name == "uniform":
            rng = random.Random(args.seed)
            records = [(f"rand{i}", random_seq(rng, 20))
                       for i in range(prof["dict_size"] + prof["num_queries"])]
            made.append(from_records("uniform_ctrl", records, out_dir,
                                     args.seed, 0, prof["dict_size"],
                                     prof["num_queries"], 0, prof["ks"]))
        else:
            sources = ([args.from_fasta] if args.from_fasta
                       else [os.path.join(RAW_DIR, s)
                             for s in resolve_sources(name)])
            if args.from_fasta is None and name in RAW_SOURCES:
                missing = [s for s in sources if not os.path.exists(s)]
                if missing:
                    print(f"skip {name}: missing raw file(s): "
                          f"{', '.join(missing)}\n  run "
                          "tests/scripts/download_datasets.py first "
                          "or pass --from-fasta")
                    continue
            if args.from_fasta is None and name not in RAW_SOURCES:
                print(f"skip {name}: needs --from-fasta "
                      "(see tests/testcases/README.md)")
                continue
            for src in sources:
                stem = os.path.splitext(os.path.basename(src))[0]
                records = read_fasta(src)
                if name == "hla":
                    # One testcase per locus file; alleles are the dictionary.
                    made.append(from_records(f"hla_{stem}", records, out_dir,
                                             args.seed, 0, prof["dict_size"],
                                             prof["num_queries"],
                                             prof["mutant_edits"], prof["ks"]))
                elif name in ("mhc", "sarscov2", "mtdna", "theseus"):
                    # Long molecule(s): tile into windows, or cohort of
                    # mutants when a single reference is given.
                    if len(records) == 1 and name in ("sarscov2", "mtdna"):
                        rng = random.Random(args.seed)
                        ref = clean(records[0][1])
                        if ref is None:
                            # e.g. rCRS carries an N spacer at position 3107.
                            ref = longest_clean_run(records[0][1])
                        if len(ref) < 20:
                            raise ValueError(f"{src}: no usable ACGT run")
                        wins = tile(ref, prof["window"])
                        cohort = []
                        for w in wins[:2]:
                            cohort.append(w)
                            for _ in range(prof["dict_size"]):
                                cohort.append(mutate(rng, w, prof["mutant_edits"]))
                        recs = [(f"v{i}", s) for i, s in enumerate(cohort)]
                        made.append(from_records(f"{name}_{stem}", recs, out_dir,
                                                 args.seed, 0, prof["dict_size"],
                                                 prof["num_queries"],
                                                 prof["mutant_edits"], prof["ks"]))
                    else:
                        recs = [(h, s) for h, s in records]
                        made.append(from_records(f"{name}_{stem}", recs, out_dir,
                                                 args.seed, prof["window"],
                                                 prof["dict_size"],
                                                 prof["num_queries"],
                                                 prof["mutant_edits"], prof["ks"]))
                else:  # markers / generic fasta: sample records directly
                    made.append(from_records(f"{name}_{stem}", records, out_dir,
                                             args.seed, prof["window"],
                                             prof["dict_size"], prof["num_queries"],
                                             prof["mutant_edits"], prof["ks"]))

    manifest = ["# Derived testcase manifest", f"seed={args.seed}", ""]
    for path, n, q, length in made:
        print(f"wrote {path} (dict={n} queries={q} len~{length})")
        manifest.append(f"- {os.path.relpath(path, DATASETS_DIR)} "
                        f"dict={n} queries={q} len~{length}")
    if made:
        out_base = args.out_dir or DERIVED_DIR
        os.makedirs(out_base, exist_ok=True)
        with open(os.path.join(out_base, "MANIFEST.md"), "w") as out:
            out.write("\n".join(manifest) + "\n")
    if not made:
        print("nothing generated")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
