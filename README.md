# HLP-grep
**H**eavy **L**ight **P**angenome **grep** (**HLP-grep**): Edit distance search queries on DNA sequence dictionary using heavy-light decomposition on pangenome paths

## Build

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build
```

## Install

```sh
cmake --install build --prefix <install-prefix>
```

## Usage

```cmake
find_package(hlp_grep REQUIRED)
target_link_libraries(your_target PRIVATE hlp_grep::hlp_grep)
```

```cpp
#include <hlp_grep/hlp_grep.hpp>

// the reference solver is directly usable
Solver solver({"ACGT", "ACGA", "TTTT"});
for (const Result &r : solver.query("ACGT", 1))
	// r.id = index in the dictionary, r.dist = edit distance
	;
```
