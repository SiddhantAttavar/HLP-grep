# HLP-grep
**H**eavy **L**ight **P**angenome **grep** (**HLP-grep**): Edit distance search queries on DNA sequence dictionary using heavy-light decomposition on pangenome paths

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

```cpp
#include <hlp_grep/hlp_grep.hpp>

Solver solver({"ACGT", "ACGA", "TTTT"});
for (const Result &r : solver.query("ACGT", 1))
	// r.id = index in the dictionary, r.dist = edit distance
	;

// custom insertion/deletion costs; the substitution matrix is
// configurable via CostModel as well
Solver custom({"ACGT", "ACGA", "TTTT"}, CostModel(2, 3));
```

Note: the heavy-light decomposition engine is still under development;
`Solver::query()` currently throws `std::logic_error`.

## Testing

Testcases are plain text files (alphabet, cost model, dictionary, queries);
the format is described in [tests/testcases/README.md](tests/testcases/README.md).

- `solve <testcase-file-or-dir>...` generates `<name>.sol` solution files
  with its built-in `NaiveSolver` reference implementation.
- `test_hlp_grep <testcase-file-or-dir>...` validates the library `Solver`
  against those solution files.

Both accept a mix of testcase files and folders containing `*.txt` files.
The ctest suite runs `test_hlp_grep` over `tests/testcases/manual` and
`tests/testcases/generated`; the manual suite is expected to fail until
`Solver::query()` is implemented.
