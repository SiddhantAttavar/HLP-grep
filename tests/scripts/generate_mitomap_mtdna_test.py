#!/usr/bin/env python3
"""Create an HLP-grep testcase from MITOMAP's full-length mtDNA set.

Example:
    python3 tests/scripts/generate_mitomap_mtdna_test.py \
        --ks 1,3,5,10,15,20 --dict-size 100 --num-queries 10

MITOMAP supplies the accession list; sequences are fetched from NCBI E-utilities
in batches and cached below tests/testcases/datasets/raw/mtdna/mitomap/. Only
full-length A/C/G/T records (>=15.4 kbp) are used. The testcase and its .hlpg
index are written beneath tests/testcases/datasets/derived/mitomap_mtdna/.
"""

import argparse
import os
import random
import re
import subprocess
import sys
import time
import urllib.error
import urllib.parse
import urllib.request


SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "..", ".."))
DATASETS_DIR = os.path.join(REPO_ROOT, "tests", "testcases", "datasets")
DEFAULT_CACHE_DIR = os.path.join(DATASETS_DIR, "raw", "mtdna", "mitomap")
DEFAULT_OUT_DIR = os.path.join(DATASETS_DIR, "derived", "mitomap_mtdna")
MITOMAP_IDS_URL = "https://mitomap.org/cgi-bin/genbank_ids.cgi"
NCBI_EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
ACCESSION_RE = re.compile(r"^[A-Z]{1,4}_?\d+(?:\.\d+)?$")
MIN_FULL_LENGTH = 15_400
ALPHABET = "AGCT"
DNA = set("ACGT")
USER_AGENT = "HLP-grep-MITOMAP-testdata/1.0 (academic research)"


def parse_ks(text):
    try:
        values = tuple(int(item.strip()) for item in text.split(","))
    except ValueError as error:
        raise argparse.ArgumentTypeError(
            "--ks must be comma-separated non-negative integers, e.g. 1,3,5"
        ) from error
    if not values or any(value < 0 for value in values):
        raise argparse.ArgumentTypeError("--ks must contain non-negative integers")
    if len(set(values)) != len(values):
        raise argparse.ArgumentTypeError("--ks values must not be repeated")
    return values


def request_text(url, data=None, timeout=90):
    request = urllib.request.Request(
        url,
        data=data,
        headers={"User-Agent": USER_AGENT, "Accept": "text/plain, */*"},
        method="POST" if data is not None else "GET",
    )
    with urllib.request.urlopen(request, timeout=timeout) as response:
        return response.read().decode("utf-8", errors="replace")


def load_accessions(cache_dir, refresh=False, accession_file=None):
    os.makedirs(cache_dir, exist_ok=True)
    path = os.path.join(cache_dir, "mitomap_full_length_accessions.txt")
    if accession_file:
        accessions = []
        with open(accession_file, encoding="utf-8") as source:
            for line in source:
                for token in line.strip().split():
                    token = token.strip(",;")
                    if ACCESSION_RE.fullmatch(token):
                        accessions.append(token)
        accessions = list(dict.fromkeys(accessions))
        if not accessions:
            raise RuntimeError(f"no GenBank accessions found in {accession_file}")
        with open(path, "w", encoding="ascii") as output:
            output.write("\n".join(accessions) + "\n")
        return accessions, path
    if os.path.isfile(path) and not refresh:
        with open(path, encoding="utf-8") as source:
            accessions = [line.strip() for line in source if line.strip()]
        if accessions:
            return accessions, path

    text = request_text(MITOMAP_IDS_URL)
    accessions = []
    for line in text.splitlines():
        value = line.strip().split()[0] if line.strip() else ""
        if ACCESSION_RE.fullmatch(value) and value.lower() != "genbank_id":
            accessions.append(value)
    # Preserve source order but discard any duplicate accession entries.
    accessions = list(dict.fromkeys(accessions))
    if not accessions:
        raise RuntimeError(
            "MITOMAP returned no GenBank accessions. Check network access or "
            f"the accession-list endpoint: {MITOMAP_IDS_URL}"
        )
    with open(path, "w", encoding="ascii") as output:
        output.write("\n".join(accessions) + "\n")
    return accessions, path


def parse_fasta(text):
    records = {}
    header = None
    chunks = []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records[header] = "".join(chunks).upper()
            # NCBI FASTA headers start with the accession.version token.
            header = line[1:].split()[0].split("|")[-1]
            chunks = []
        elif header is not None:
            chunks.append("".join(line.split()))
    if header is not None:
        records[header] = "".join(chunks).upper()
    return records


def iter_fasta(path):
    header = None
    chunks = []
    with open(path, encoding="ascii") as source:
        for line in source:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks).upper()
                header = line[1:].split()[0].split("|")[-1]
                chunks = []
            elif header is not None:
                chunks.append("".join(line.split()))
    if header is not None:
        yield header, "".join(chunks).upper()


def append_sequence_cache(path, records):
    with open(path, "a", encoding="ascii") as output:
        for accession, sequence in records.items():
            output.write(f">{accession}\n")
            for start in range(0, len(sequence), 80):
                output.write(sequence[start:start + 80] + "\n")


def fetch_batch(accessions):
    payload = urllib.parse.urlencode({
        "db": "nuccore",
        "id": ",".join(accessions),
        "rettype": "fasta",
        "retmode": "text",
    }).encode("ascii")
    return parse_fasta(request_text(NCBI_EFETCH_URL, payload))


def is_usable(sequence):
    return (len(sequence) >= MIN_FULL_LENGTH
            and all(base in DNA for base in sequence))


def select_records(accessions, wanted, rng, cache_dir, batch_size=40):
    cache_path = os.path.join(cache_dir, "full_length_sequences.fasta")
    accession_set = set(accessions)
    cached_accessions = set()
    selected = []
    usable_seen = 0

    def consider(accession, sequence):
        nonlocal usable_seen
        if accession not in accession_set or not is_usable(sequence):
            return
        usable_seen += 1
        if len(selected) < wanted:
            selected.append((accession, sequence))
        else:
            position = rng.randrange(usable_seen)
            if position < wanted:
                selected[position] = (accession, sequence)

    if os.path.isfile(cache_path):
        for accession, sequence in iter_fasta(cache_path):
            cached_accessions.add(accession)
            consider(accession, sequence)
        if len(selected) >= wanted:
            rng.shuffle(selected)
            return selected, 0

    candidates = list(accessions)
    rng.shuffle(candidates)
    fetched = 0

    for start in range(0, len(candidates), batch_size):
        if len(selected) >= wanted:
            break
        batch = [accession for accession in candidates[start:start + batch_size]
                 if accession not in cached_accessions]
        if batch:
            try:
                downloaded = fetch_batch(batch)
            except (OSError, urllib.error.URLError, TimeoutError) as error:
                raise RuntimeError(
                    f"failed fetching NCBI sequences for batch beginning "
                    f"{batch[0]}: {error}"
                ) from error
            append_sequence_cache(cache_path, downloaded)
            cached_accessions.update(downloaded)
            fetched += len(downloaded)
            for accession, sequence in downloaded.items():
                consider(accession, sequence)
            print(f"fetched {len(downloaded)} record(s); "
                  f"usable cache now has {usable_seen}", flush=True)
            # NCBI recommends spaced requests; one batch typically covers many IDs.
            time.sleep(0.4)

    if len(selected) < wanted:
        raise RuntimeError(
            f"only {len(selected)} usable full-length ACGT records found; "
            f"need {wanted}. MITOMAP accession list has {len(accessions)} IDs."
        )
    rng.shuffle(selected)
    return selected, fetched


def write_testcase(path, dictionary, queries):
    with open(path, "w", encoding="ascii") as output:
        output.write(f"{ALPHABET}\n1 1 0 1\n")
        output.write(f"{len(dictionary)}\n")
        output.writelines(sequence + "\n" for _accession, sequence in dictionary)
        output.write(f"{len(queries)}\n")
        output.writelines(f"{k} {sequence}\n" for k, sequence in queries)


def save_index(save_bin, testcase):
    environment = os.environ.copy()
    environment["OMP_NUM_THREADS"] = "1"
    environment["OMP_DYNAMIC"] = "FALSE"
    subprocess.run([save_bin, testcase], check=True, env=environment)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ks", required=True, type=parse_ks,
                        help="comma-separated thresholds, e.g. 1,3,5,10,15,20")
    parser.add_argument("--dict-size", required=True, type=int,
                        help="number of dictionary sequences")
    parser.add_argument("--num-queries", required=True, type=int,
                        help="number of distinct held-out sequences; each is tested at every k")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--out-dir", default=DEFAULT_OUT_DIR)
    parser.add_argument("--cache-dir", default=DEFAULT_CACHE_DIR)
    parser.add_argument("--accession-file", default=None,
                        help="use an already-downloaded MITOMAP GenBank ID list")
    parser.add_argument("--save-bin", default=os.path.join(REPO_ROOT, "build", "tests", "save"))
    parser.add_argument("--refresh-accessions", action="store_true",
                        help="refresh MITOMAP's accession-list cache")
    parser.add_argument("--overwrite", action="store_true",
                        help="replace an existing testcase, index, and manifest")
    args = parser.parse_args()
    if args.dict_size < 1 or args.num_queries < 1:
        parser.error("--dict-size and --num-queries must be positive")
    if not os.path.isfile(args.save_bin) or not os.access(args.save_bin, os.X_OK):
        parser.error(f"save executable not found or not executable: {args.save_bin}; "
                     "build it with `cmake --build build --target save` or pass --save-bin")

    os.makedirs(args.out_dir, exist_ok=True)
    basename = (f"mitomap_mtdna_d{args.dict_size}_q{args.num_queries}_"
                f"k{'-'.join(map(str, args.ks))}_seed{args.seed}")
    testcase = os.path.join(args.out_dir, basename + ".txt")
    index = os.path.splitext(testcase)[0] + ".hlpg"
    manifest = os.path.join(args.out_dir, "MANIFEST.md")
    if not args.overwrite and any(os.path.exists(path)
                                  for path in (testcase, index)):
        parser.error(f"output already exists; choose another --out-dir or use --overwrite: {testcase}")

    try:
        accessions, accession_path = load_accessions(
            args.cache_dir, refresh=args.refresh_accessions,
            accession_file=args.accession_file)
        rng = random.Random(args.seed)
        records, fetched = select_records(
            accessions, args.dict_size + args.num_queries, rng, args.cache_dir)
    except (OSError, RuntimeError, urllib.error.URLError) as error:
        parser.error(str(error))

    dictionary_records = records[:args.dict_size]
    heldout = records[args.dict_size:]
    query_pairs = [(k, sequence) for _accession, sequence in heldout
                   for k in args.ks]
    write_testcase(testcase, dictionary_records, query_pairs)
    save_index(args.save_bin, testcase)

    with open(manifest, "w", encoding="utf-8") as output:
        output.write("# MITOMAP human mtDNA testcase\n\n")
        output.write(f"- Seed: `{args.seed}`\n")
        output.write(f"- Dictionary sequences: `{len(dictionary_records)}`\n")
        output.write(f"- Held-out sequences: `{len(heldout)}`\n")
        output.write(f"- Thresholds: `{','.join(map(str, args.ks))}`\n")
        output.write(f"- Testcase queries: `{len(query_pairs)}`\n")
        output.write(f"- Newly fetched records: `{fetched}`\n")
        output.write(f"- MITOMAP accession list: `{accession_path}`\n")
        output.write("\n## Accessions\n\n")
        output.write("### Dictionary\n\n")
        output.write("\n".join(f"- `{accession}`" for accession, _ in dictionary_records))
        output.write("\n\n### Queries\n\n")
        output.write("\n".join(f"- `{accession}`" for accession, _ in heldout))
        output.write("\n")

    print(f"wrote {testcase} (dict={len(dictionary_records)} "
          f"heldout={len(heldout)} query-lines={len(query_pairs)})")
    print(f"saved index {index}")
    print(f"manifest: {manifest}")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except subprocess.CalledProcessError as error:
        print(f"save failed with exit code {error.returncode}", file=sys.stderr)
        sys.exit(error.returncode or 1)
