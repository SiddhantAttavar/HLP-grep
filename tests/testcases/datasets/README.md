# Datasets

Real-world highly similar DNA for threshold-search testing. Documented in
`tests/testcases/README.md` ("Real-world datasets" section).

- `raw/` — downloaded sources (git-ignored).
- `derived/` — generated `.txt`/`.sol` testcases + `MANIFEST.md` (git-ignored).
- `logs/` — `test_datasets.sh` logs (git-ignored).

Quickstart:

```bash
python3 tests/scripts/download_datasets.py
python3 tests/scripts/generate_dataset_tests.py --dataset synthetic
tests/scripts/test_datasets.sh
```
