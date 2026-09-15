# HLP-grep
**H**eavy **L**ight **P**artial-order-alignment **grep** (**HLP-grep**): Edit distance search queries on DNA sequence dictionary using heavy-light decomposition on partial order alignment (POA) graph paths

## Build

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build
```

The build also produces two test executables, `build/tests/solve` and
`build/tests/test_hlp_grep` (see [Testing](#testing)).

## Install

```sh
cmake --install build --prefix <install-prefix>
```

## Usage

Header-only library:

```cmake
find_package(hlp_grep REQUIRED)
target_link_libraries(your_target PRIVATE hlp_grep::hlp_grep)
```

All types live in namespace `hlp_grep` and are available through the
single public API header `hlp_grep/hlp_grep.hpp`. The query answer type
carries the match:

```cpp
struct hlp_grep::Result {
	std::size_t id; ///< index of the matched sequence in the dictionary (0-based)
	int dist;       ///< edit distance between the match and the query string
};
```

### Basic query

Construct a `hlp_grep::Solver` over the dictionary, then query with a
string and a threshold `k`. The result vector enumerates every
dictionary sequence whose edit distance to the query is `<= k`, sorted
by id.

```cpp
#include <hlp_grep/hlp_grep.hpp>
#include <iostream>
#include <string>
#include <vector>

int main() {
	const std::vector<std::string> dict{
	    "ACGT", // id 0
	    "ACGA", // id 1
	    "TTTT", // id 2
	};
	hlp_grep::Solver solver(dict);

	for (const hlp_grep::Result &r : solver.query("ACGT", 1))
		std::cout << "dict[" << r.id << "] dist " << r.dist << '\n';
	// dist("ACGT", "ACGT") = 0
	// dist("ACGA", "ACGT") = 1 -> replaces the final A with T

	// Raise the threshold and re-query freely; the dictionary graph is
	// built once per Solver, only the per-query precomputation repeats.
}
```

Output:

```text
dict[0] dist 0
dict[1] dist 1
```

`Result::id` is a 0-based index into the `dict` vector as passed to the
constructor. On disk, testcase solution files store the same ids as
1-based values and the test harness converts them.

### Cost models

`CostModel(ins, del, alphabet, matrix)` defines the edit operation
costs: `ins` and `del` default to 1, and the substitution matrix
defaults to the unit model over the DNA alphabet `AGCT` (0 on the
diagonal, 1 elsewhere). The matrix is looked up as
`CostModel::match(alphabet[i], alphabet[j])`; characters absent from
the alphabet fall back to index 0.

```cpp
// Heavier noise operations: alignments prefer matching bases and
// thresholds stay strict (ins = 2, del = 3, unit substitutions).
hlp_grep::CostModel strict(2, 3, "AGCT", {});

// Minimal-cost transitions only between purines (A<->G) and
// pyrimidines (C<->T); rows/columns follow the alphabet "AGCT".
const std::vector<std::vector<int>> purine_pyrimidine{
    {0, 4, 1, 1}, //      A  G  C  T
    {4, 0, 1, 1},
    {1, 1, 0, 4},
    {1, 1, 4, 0},
};
hlp_grep::CostModel transitions(1, 1, "AGCT", purine_pyrimidine);

hlp_grep::Solver solver(dict, transitions);
```

hlp_grep::Solver solver(dict, transitions);
```

## Testing

Testcases are plain text files (alphabet, cost model, dictionary, queries);
the format is described in [tests/testcases/README.md](tests/testcases/README.md).

- `solve <testcase-file-or-dir>...` generates `<name>.sol` solution files
  with its built-in `NaiveSolver` reference implementation.
- `test_hlp_grep <testcase-file-or-dir>...` validates the library `Solver`
  against those solution files.

Both accept a mix of testcase files and folders containing `*.txt` files.
The ctest suite also runs unit tests for `DistMatrix` and the POA graph.
The `tests/testcases/manual` suite covers the implemented `Solver::query()`.
