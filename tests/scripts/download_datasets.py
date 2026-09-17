#!/usr/bin/env python3
"""Download real-world DNA datasets for HLP-grep threshold-search testing.

Default mode downloads small curated subsets (KBs, seconds). Full datasets
(MBs/GBs, minutes to hours) require the explicit `--full` flag, and the
largest ones additionally require `--yes`.

    python3 tests/scripts/download_datasets.py --dry-run
    python3 tests/scripts/download_datasets.py
    python3 tests/scripts/download_datasets.py --dataset hla
    python3 tests/scripts/download_datasets.py --full --dataset hla --yes

Layout (git-ignored except .gitkeep files):
    tests/testcases/datasets/raw/      downloaded sources
    tests/testcases/datasets/derived/  generated .txt/.sol testcases
    tests/testcases/datasets/logs/     runner logs

Sources and licences are documented in tests/testcases/README.md.
Only the standard library is used.
"""

import argparse
import os
import sys
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
DATASETS_DIR = os.path.normpath(os.path.join(HERE, "..", "testcases", "datasets"))
RAW_DIR = os.path.join(DATASETS_DIR, "raw")

IMGT_BASE = "https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/fasta"
NCBI_EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"


def ncbi_url(accession, start=None, stop=None):
    url = (
        f"{NCBI_EFETCH}?db=nucleotide&id={accession}"
        f"&rettype=fasta&retmode=text"
    )
    if start is not None and stop is not None:
        url += f"&seq_start={start}&seq_stop={stop}"
    return url


# Each entry: (relative raw path, url, approx bytes, needs_full, needs_yes)
FILES = {
    # Curated HLA: small per-locus cDNA alignments as FASTA (real alleles,
    # KBs). DPA2 has 6 alleles, DPB2 7, DQA2/DQB2/DMA are larger panels.
    "hla": [
        ("hla/DMA_nuc.fasta", f"{IMGT_BASE}/DMA_nuc.fasta", 88665, False, False),
        ("hla/DPA2_nuc.fasta", f"{IMGT_BASE}/DPA2_nuc.fasta", 4701, False, False),
        ("hla/DPB2_nuc.fasta", f"{IMGT_BASE}/DPB2_nuc.fasta", 5771, False, False),
        ("hla/DQA2_nuc.fasta", f"{IMGT_BASE}/DQA2_nuc.fasta", 34341, False, False),
        ("hla/DQB2_nuc.fasta", f"{IMGT_BASE}/DQB2_nuc.fasta", 34620, False, False),
        # Full: highly polymorphic loci (MBs each).
        ("hla/A_nuc.fasta", f"{IMGT_BASE}/A_nuc.fasta", 9291105, True, False),
        ("hla/B_nuc.fasta", f"{IMGT_BASE}/B_nuc.fasta", 11171138, True, False),
        ("hla/C_nuc.fasta", f"{IMGT_BASE}/C_nuc.fasta", 9635306, True, False),
        ("hla/DRB1_nuc.fasta", f"{IMGT_BASE}/DRB1_nuc.fasta", 2215640, True, False),
        ("hla/DQB1_nuc.fasta", f"{IMGT_BASE}/DQB1_nuc.fasta", 2110526, True, False),
        ("hla/DPB1_nuc.fasta", f"{IMGT_BASE}/DPB1_nuc.fasta", 2201676, True, False),
        ("hla/DQA1_nuc.fasta", f"{IMGT_BASE}/DQA1_nuc.fasta", 823390, True, False),
        ("hla/DPA1_nuc.fasta", f"{IMGT_BASE}/DPA1_nuc.fasta", 749625, True, False),
    ],
    # Curated MHC: two 10 kb windows of finished haplotype APD (OK649231,
    # 4.9 Mbp). Full: six finished haplotypes, ~30 MB total.
    "mhc": [
        ("mhc/OK649231_w0_10kb.fasta", ncbi_url("OK649231", 1, 10000), 10200, False, False),
        ("mhc/OK649231_w2500k_10kb.fasta", ncbi_url("OK649231", 2500001, 2510000), 10200, False, False),
        ("mhc/OK649231.fasta", ncbi_url("OK649231"), 5000000, True, True),
        ("mhc/OK649232.fasta", ncbi_url("OK649232"), 5100000, True, True),
        ("mhc/OK649233.fasta", ncbi_url("OK649233"), 5000000, True, True),
        ("mhc/OK649234.fasta", ncbi_url("OK649234"), 5100000, True, True),
        ("mhc/OK649235.fasta", ncbi_url("OK649235"), 5000000, True, True),
        ("mhc/OK649236.fasta", ncbi_url("OK649236"), 5100000, True, True),
    ],
    # Curated SARS-CoV-2: Wuhan-Hu-1 reference (29.9 kb). Full mode additionally
    # documents the Theseus 2732-genome set (manual download, see README).
    "sarscov2": [
        ("sarscov2/NC_045512.fasta", ncbi_url("NC_045512.2"), 31000, False, False),
    ],
    # Curated mtDNA: rCRS reference (16.6 kb).
    "mtdna": [
        ("mtdna/NC_012920.fasta", ncbi_url("NC_012920.1"), 17500, False, False),
    ],
}

DATASETS = sorted(FILES)


def plan(dataset, full):
    items = []
    for name in ([dataset] if dataset != "all" else DATASETS):
        for rel, url, size, needs_full, needs_yes in FILES[name]:
            if needs_full and not full:
                continue
            items.append((name, rel, url, size, needs_yes))
    return items


def download(url, dest):
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    req = urllib.request.Request(url, headers={"User-Agent": "HLP-grep-testdata/0.1"})
    with urllib.request.urlopen(req, timeout=120) as resp, open(dest, "wb") as out:
        while True:
            chunk = resp.read(1 << 20)
            if not chunk:
                break
            out.write(chunk)
    return os.path.getsize(dest)


def main():
    parser = argparse.ArgumentParser(description="Download HLP-grep test datasets")
    parser.add_argument("--dataset", default="all", choices=["all"] + DATASETS)
    parser.add_argument("--full", action="store_true",
                        help="include full (large) dataset files")
    parser.add_argument("--yes", action="store_true",
                        help="confirm the largest full downloads")
    parser.add_argument("--dry-run", action="store_true",
                        help="print the download manifest without fetching")
    parser.add_argument("--force", action="store_true",
                        help="re-download files that already exist")
    args = parser.parse_args()

    items = plan(args.dataset, args.full)
    need_yes = [rel for _, rel, _, _, needs_yes in items if needs_yes]
    if need_yes and args.full and not args.yes:
        print("Refusing: the following full downloads need --yes "
              f"(~30 MB MHC haplotypes): {', '.join(need_yes)}")
        return 1

    total = sum(size for _, _, _, size, _ in items)
    print(f"mode={'full' if args.full else 'curated'} "
          f"files={len(items)} approx={total / 1e6:.1f} MB -> {RAW_DIR}")
    if args.dry_run:
        for name, rel, url, size, _ in items:
            print(f"  [{name}] {rel} (~{size} B)\n    {url}")
        return 0

    failed = 0
    for name, rel, url, size, _ in items:
        dest = os.path.join(RAW_DIR, rel)
        if os.path.exists(dest) and not args.force:
            print(f"  skip (exists): {rel}")
            continue
        try:
            got = download(url, dest)
            print(f"  ok: {rel} ({got} B)")
        except Exception as exc:  # noqa: BLE001 - report URL, keep going
            print(f"  FAILED: {rel}\n    {url}\n    {exc}")
            failed += 1
    if failed:
        print(f"{failed} download(s) failed; see tests/testcases/README.md "
              "for manual alternatives")
        return 1
    print("done")
    return 0


if __name__ == "__main__":
    sys.exit(main())
