# Testcase format

Each testcase is a single text file containing a single dictionary of strings followed by a set of queries. Additionally, the alphabet and cost model are specified.

## Structure
1. **Alphabet**: A single string containing the characters in the alphabet $\Sigma$.
2. **Cost model**: A single line containing 4 space-separated integers: insertion cost, deletion cost, match cost, mismatch cost.
3. **Dictionary header**: A single line containing the number of strings $n$ in the dictionary.
4. **Dictionary sequences**: one string per line $dict_i$.
5. **Query header**: A single integer $q$ giving the number of queries.
6. **Queries**: one query per line; the value of $k$ (the edit distance threshold) followed by the query sequence $query_i$, separated by a space.

## Example

```text
AGCT
1 1 0 1
4
ACGT
ACGA
TTTT
ACG
3
1 ACGT
0 TTTT
2 AAAA
```

- The alphabet used is ${A, G, C, T}$
- The cost of insertion, deletion and mismatch are all $1$, and match costs $0$
- The dictionary contains 4 sequences.
- There are 3 queries:
  - `ACGT` with $k = 1$,
  - `TTTT` with $k = 0$,
  - `AAAA` with $k = 2$.

## Testcase structure
The `testcases` folder contains 3 subfolders:
1. `manual`: Short manually-created testcases intended for debugging
2. `generated`: Larger testscases generated using `testcase_gen.py`
3. `datasets`: Real-world highly similar DNA datasets (see below)

## Real-world datasets

`datasets/` holds downloaded sources and derived testcases for threshold
search over highly similar dictionaries — the regime where POA heavy chains
are shared. Only scripts and docs are tracked; `raw/`, `derived/` and
`logs/` contents are git-ignored (see root `.gitignore`):

```text
tests/testcases/datasets/
  raw/      downloaded sources (per-dataset subdirs)
  derived/  generated .txt testcases + .sol solutions + MANIFEST.md
  logs/     test_datasets.sh run logs
tests/scripts/
  download_datasets.py    fetch sources (curated by default)
  generate_dataset_tests.py  FASTA -> testcase .txt files
  test_datasets.sh        generate + solve + validate end to end
```

### Download

Curated subsets (KBs, seconds) are the default; full collections need
`--full`, the largest additionally need `--yes`:

```bash
python3 tests/scripts/download_datasets.py --dry-run
python3 tests/scripts/download_datasets.py
python3 tests/scripts/download_datasets.py --dataset hla
python3 tests/scripts/download_datasets.py --full --dataset hla --yes
```

| Dataset | Curated subset | Full (`--full`) | Source / licence |
|---|---|---|---|
| `hla` | `DMA/DPA2/DPB2/DQA2/DQB2_nuc.fasta` (~170 KB real alleles) | `A/B/C/DRB1/DQB1/DPB1/DQA1/DPA1_nuc.fasta` (~40 MB) | IPD-IMGT/HLA, GitHub `ANHIG/IMGTHLA` (`Latest`), cite Robinson et al.; free for academic use |
| `mhc` | Two 10 kb windows of finished haplotype APD (`OK649231`) | Six finished haplotypes `OK649231-OK649236` (~30 MB) | GenBank (INSDC, open) |
| `sarscov2` | Wuhan-Hu-1 `NC_045512.2` (29.9 kb) | Theseus `covid_19_complete.fasta` (2732 genomes, manual, Zenodo `18482097`) | GenBank (open); Theseus per original licences |
| `mtdna` | rCRS `NC_012920.1` (16.6 kb) | MITOMAP / GenBank bulk (manual, 65k+ sequences) | GenBank (open); MITOMAP cite Brandon et al. |
| `markers`/`theseus` | Bring your own FASTA via `--from-fasta` | Same | 16S/COX1 GenBank loci; Theseus Zenodo sets (MTB/HIV/monkeypox) |

### Generate

```bash
python3 tests/scripts/generate_dataset_tests.py --list
python3 tests/scripts/generate_dataset_tests.py --dataset synthetic
python3 tests/scripts/generate_dataset_tests.py --dataset hla --seed 7
python3 tests/scripts/generate_dataset_tests.py --dataset markers \
  --from-fasta path/to/16s.fasta --dict-size 10 --ks 0,1,3
```

Per-dataset profiles (window, dict size, thresholds) live in
`PROFILES` in `generate_dataset_tests.py`. Conventions: one testcase per
HLA locus file (alleles are the dictionary); long molecules
(MHC/SARS-CoV-2/mtDNA) are tiled into windows or expanded into mutant
cohorts so the `O(n·q·L²)` naive solver finishes in seconds; queries mix
exact (`k=0`), held-out and synthetic-mutant sequences. Records with
non-`ACGT` bases are dropped and reported.

### Test

```bash
tests/scripts/test_datasets.sh            # quick: generate + solve + validate
tests/scripts/test_datasets.sh --full --validate  # larger profiles
build/tests/solve tests/testcases/datasets/derived/hla
build/tests/test_hlp_grep tests/testcases/datasets/derived/hla
```

Why these sets: alleles of one HLA gene differ by 1–tens of edits
(maximal chain sharing); MHC haplotype windows stress long chains;
SARS-CoV-2/mtDNA give 16–30 kb near-identical cohorts; `synthetic`
isolates sharing vs length vs `k`; `uniform` is the low-sharing control.

## *Optional*: testcase solution
For each testcase file `{testcase_name}.txt`, solutions may be stored in `{testcase_name}.sol`. The solution file contains 3 lines for each query
1. Number of sequences $l$ in the dictionary satisfying $edit\_dist(query, dict_j) \le k$
2. $l$ space separated integers $a_1, a_2, \cdots a_l$ ($1 \le a_1 \le a_2 \cdots a_l \le |dict|$) representing the $1$-based index of strings $edit\_dist(query, dict_j) \le k$
3. $l$ space separated integers representing the edit distances of the strings given in the previous line

Solution files are generated from testcases with the `solve` executable
(`build/tests/solve <testcase-file-or-dir>...`), which uses its built-in
`NaiveSolver` reference implementation. They are validated with the
`test_hlp_grep` executable (`build/tests/test_hlp_grep
<testcase-file-or-dir>...`), which runs the library `Solver`. Both accept a
mix of testcase files and folders containing `*.txt` files.
