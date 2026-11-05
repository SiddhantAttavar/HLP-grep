#!/usr/bin/env python3
"""Download MITOMAP's full-length human mtDNA GenBank set as FASTA.

The MITOMAP accession-list endpoint is used by default. If it is inaccessible,
download the "GenBank FL ID set" file from https://www.mitomap.org/MITOMAP/Mitobank
in a browser and pass its path with --accession-file.

    python3 tests/scripts/download_mitomap_mtdna.py
    python3 tests/scripts/download_mitomap_mtdna.py --accession-file ~/Downloads/genbank_ids.txt

Sequences are fetched from NCBI E-utilities in batches and appended to a
restartable FASTA cache under tests/testcases/datasets/raw/mtdna/mitomap/.
"""

import argparse
import os
import sys
import time
import urllib.error

from generate_mitomap_mtdna_test import (
    DEFAULT_CACHE_DIR,
    is_usable,
    load_accessions,
    parse_fasta,
    request_text,
)


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
                chunks.append(line)
    if header is not None:
        yield header, "".join(chunks).upper()


def fetch_batch(accessions):
    import urllib.parse
    import urllib.request

    payload = urllib.parse.urlencode({
        "db": "nuccore",
        "id": ",".join(accessions),
        "rettype": "fasta",
        "retmode": "text",
    }).encode("ascii")
    return parse_fasta(request_text(
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi", payload,
        timeout=180))


def append_records(path, records):
    with open(path, "a", encoding="ascii") as output:
        for accession, sequence in records.items():
            output.write(f">{accession}\n")
            for start in range(0, len(sequence), 80):
                output.write(sequence[start:start + 80] + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--accession-file", default=None,
                        help="locally downloaded MITOMAP full-length ID list")
    parser.add_argument("--cache-dir", default=DEFAULT_CACHE_DIR)
    parser.add_argument("--batch-size", type=int, default=100,
                        help="accessions per NCBI EFetch request (default: 100)")
    parser.add_argument("--delay", type=float, default=0.4,
                        help="seconds between NCBI requests (default: 0.4)")
    parser.add_argument("--refresh-accessions", action="store_true",
                        help="refresh the accession list from MITOMAP")
    args = parser.parse_args()
    if args.batch_size < 1 or args.batch_size > 500:
        parser.error("--batch-size must be between 1 and 500")
    if args.delay < 0:
        parser.error("--delay must be non-negative")

    os.makedirs(args.cache_dir, exist_ok=True)
    sequence_path = os.path.join(args.cache_dir, "full_length_sequences.fasta")
    try:
        accessions, accession_path = load_accessions(
            args.cache_dir, refresh=args.refresh_accessions,
            accession_file=args.accession_file)
    except (OSError, RuntimeError, urllib.error.URLError) as error:
        parser.error(str(error))

    cached = set()
    if os.path.isfile(sequence_path):
        for accession, _sequence in iter_fasta(sequence_path):
            cached.add(accession)
    pending = [accession for accession in accessions if accession not in cached]
    print(f"MITOMAP accessions: {len(accessions)}; cached: {len(cached)}; "
          f"to fetch: {len(pending)}")
    if pending:
        for start in range(0, len(pending), args.batch_size):
            batch = pending[start:start + args.batch_size]
            error = None
            records = {}
            for attempt in range(4):
                try:
                    records = fetch_batch(batch)
                    error = None
                    break
                except (OSError, urllib.error.URLError, TimeoutError) as caught:
                    error = caught
                    time.sleep(min(2 ** attempt, 8))
            if error is not None:
                parser.error(f"NCBI EFetch failed for batch starting at "
                             f"{batch[0]}: {error}")
            append_records(sequence_path, records)
            cached.update(records)
            if (start // args.batch_size + 1) % 10 == 0 or start + len(batch) >= len(pending):
                print(f"fetched {min(start + len(batch), len(pending))}/"
                      f"{len(pending)} accessions; received {len(cached)} total records",
                      flush=True)
            if args.delay:
                time.sleep(args.delay)

    # Count coverage/usability by streaming rather than loading ~1 GB into RAM.
    total_records = 0
    usable_records = 0
    total_bases = 0
    for _accession, sequence in iter_fasta(sequence_path):
        total_records += 1
        total_bases += len(sequence)
        if is_usable(sequence):
            usable_records += 1

    missing = len(accessions) - len(cached.intersection(accessions))
    print(f"FASTA: {sequence_path}")
    print(f"accession list: {accession_path}")
    print(f"records saved: {total_records}; usable ACGT full-length: {usable_records}; "
          f"mean length: {total_bases / max(total_records, 1):.1f} bp; "
          f"accessions not returned: {missing}")
    if missing:
        print("Some accessions were not returned by NCBI; rerun to retry the missing IDs.",
              file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
