# HLP-grep
**H**eavy **L**ight **P**artial-order-alignment **grep** (**HLP-grep**): Edit distance search queries on DNA sequence dictionary using heavy-light decomposition on partial order alignment (POA) graph paths

## Build

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build
```

The build also produces the test executables under `build/tests/`
(`solve`, `test_hlp_grep`, plus `POAGraph` and `DistMatrix` unit tests;
see [Testing](#testing)).

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

### Save and load the graph

`Solver::save_graph()` writes the finalized POA graph, including node labels,
edges, paths, topological order, and heavy-chain metadata. `Solver::load_graph()`
reconstructs a solver from that file without rebuilding the dictionary:

```cpp
hlp_grep::Solver first(dict);
first.save_graph("index.hlpg");

hlp_grep::Solver cached = hlp_grep::Solver::load_graph("index.hlpg");
auto matches = cached.query(query, k);
```

The cache is a versioned binary format and includes its cost model; loading
with a different model fails. Solver query helpers are reconstructed from
the saved graph, while query-specific tables are still built for each query.
`POAGraph` also provides `save_file()`, `load_file()`, and stream-based
`save()`/`load()` methods for graph-only use.

The `save` and `load` helpers apply this to testcase files: `save` builds
the graph for a testcase dictionary and writes the `.hlpg` file, while
`load` answers the testcase queries from the saved graph (skipping the
dictionary build) and writes the `.sol` solution file:

```sh
build/tests/save <testcase>.txt --out index.hlpg
build/tests/load <testcase>.txt index.hlpg --out <testcase>.sol
```

`bench` can also run the `hlp_grep` method from a prebuilt index: with
`--index`, the dictionary build is skipped and the reported build time is
zero, isolating query-time comparisons:

```sh
build/tests/bench --method hlp_grep <testcase>.txt --index index.hlpg
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

`CostModel(ins, del, match, mismatch)` holds the costs of the basic edit
operations: consuming equal characters costs `match` (default 0), any
distinct pair costs `mismatch` (default 1), and `ins`/`del` default to 1.
`Solver` and the other consumers (`POAGraph`, `BinaryLifter`,
`NaiveSolver`) hold a `CostModel&` reference, so construct the model as a
named object and keep it alive for as long as the consumer is used
(`DEFAULT_COST_MODEL` is a permanent unit-cost object used as the default
argument when no explicit model is supplied).

The constructor enforces the structural precondition behind the
argmin-staircase composition used by `DistMatrix::min_plus_product`
(non-crossing shortest paths in the edit-distance DAG): `match == 0`,
positive `ins`/`del`, nonnegative `mismatch`, and
`ins + del >= mismatch`; violating models
throw `std::invalid_argument`.

```cpp
// Heavier noise operations: alignments prefer matching bases and
// thresholds stay strict (ins = 2, del = 3, unit substitutions).
hlp_grep::CostModel strict(2, 3);

hlp_grep::Solver solver(dict, strict);
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

## Performance flags

All default to off so default builds stay portable; enable them explicitly
when benchmarking on a fixed machine:

- `-DHLP_GREP_MARCH=ON`: `-march=native` (GNU/Clang) for all `hlp_grep`
  consumers (unlocks AVX-512 etc. on this machine).
- `-DHLP_GREP_LTO=ON`: `-flto` at compile and link time.
- `-DHLP_GREP_PGO=GENERATE|USE`: profile-guided optimization. Two-phase
  workflow against a representative workload (the same flags must match in
  both phases; keep `HLP_GREP_MARCH/LTO` identical):
  ```bash
  cmake -S . -B build-pgo -DCMAKE_BUILD_TYPE=Release -DHLP_GREP_MARCH=ON \
    -DHLP_GREP_LTO=ON -DHLP_GREP_PGO=GENERATE
  cmake --build build-pgo
  build-pgo/tests/bench --method hlp_grep <representative-testcase>
  cmake -S . -B build-pgo -DHLP_GREP_PGO=USE   # same dir, reuse .gcda
  cmake --build build-pgo
  ```
