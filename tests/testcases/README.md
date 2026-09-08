# Testcase format

Each testcase is a single text file containing a single dictionary of strings followed by a set of queries. Additionally, the alphabet and cost model are specified.

## Structure
1. **Alphabet**: A single string containing the characters in the alphabet $\Sigma$.
2. **Cost model**: 
   1. The first line contains two integers: the cost of insertion and deletion respectively.
   2. The next $|\Sigma|$ lines containing $|\Sigma|$ integers each. The $j^{th}$ integer in the $i^{th}$ line is the cost of match/mismatch between $\Sigma_i$ and $\Sigma_j$. This matrix must be symmetric, and all elements along the diagonal must be $0$.
3. **Dictionary header**: A single line containing the number of strings $n$ in the dictionary.
4. **Dictionary sequences**: one string per line $dict_i$.
5. **Query header**: A single integer $q$ giving the number of queries.
6. **Queries**: one query per line; the value of $k$ (the edit distance threshold) followed by the query sequence $query_i$, separated by a space.

## Example

```text
AGCT
1 1
0 1 1 1
1 0 1 1
1 1 0 1
1 1 1 0
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
- The cost of insertion, deletion and mismatch are all $1$
- The dictionary contains 4 sequences.
- There are 3 queries:
  - `ACGT` with $k = 1$,
  - `TTTT` with $k = 0$,
  - `AAAA` with $k = 2$.

## Testcase structure
The `testcases` folder contains 3 subfolders:
1. `manual`: Short manually-created testcases intended for debugging
2. `generated`: Larger testscases generated using `testcase_gen.py`

## *Optional*: testcase solution
For each testcase file `{testcase_name}.txt`, solutions may be stored in `{testcase_name}.sol`. The solution file contains 3 lines for each query
1. Number of sequences $l$ in the dictionary satisfying $edit\_dist(query, dict_j) \le k$
2. $l$ space separated integers $a_1, a_2, \cdots a_l$ ($1 \le a_1 \le a_2 \cdots a_l \le |dict|$) representing the $1$-based index of strings $edit\_dist(query, dict_j) \le k$
3. $l$ space separated integers representing the edit distances of the strings given in the previous line

Solution files are generated from testcases with the `solve` executable (`build/tests/solve <solver-name> <testcase-file-or-dir>...`) and are validated against with the `test_hlp_grep` executable (`build/tests/test_hlp_grep <solver-name> <testcase-file-or-dir>...`). Both accept a mix of testcase files and folders containing `*.txt` files. Currently supported solvers: `NaiveSolver`.
