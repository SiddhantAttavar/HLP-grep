add_executable(solve naive_solver.cpp)
target_link_libraries(solve PRIVATE hlp_grep::hlp_grep)

add_executable(save save.cpp)
target_link_libraries(save PRIVATE hlp_grep::hlp_grep)

add_executable(load load.cpp)
target_link_libraries(load PRIVATE hlp_grep::hlp_grep)

option(HLP_GREP_WITH_WFA2 "Enable the WFA2-lib method in bench (needs third_party/WFA2-lib)" ON)
option(HLP_GREP_WITH_DT_PATRICIA "Enable the DT-Patricia method in bench (needs third_party/DT-Patricia)" ON)
option(HLP_GREP_WITH_BED_TREE "Enable the BED-tree method in bench (needs third_party/bed-tree)" ON)
option(HLP_GREP_WITH_HSTREE "Enable the hstree method in bench (needs third_party/Similarity-Search-and-Join)" ON)
option(HLP_GREP_WITH_EDLIB "Enable the edlib method in bench (needs third_party/edlib)" ON)
option(HLP_GREP_WITH_PARASAIL "Enable the parasail method in bench (needs third_party/parasail)" ON)

add_executable(test_hlp_grep test.cpp)
target_link_libraries(test_hlp_grep PRIVATE hlp_grep::hlp_grep)

add_executable(poa_test poa_test.cpp)
target_link_libraries(poa_test PRIVATE hlp_grep::hlp_grep)

add_executable(dist_matrix_test dist_matrix_test.cpp)
target_link_libraries(dist_matrix_test PRIVATE hlp_grep::hlp_grep)

add_executable(bench_poa bench_poa.cpp)
target_link_libraries(bench_poa PRIVATE hlp_grep::hlp_grep)

add_executable(bench bench.cpp)
target_link_libraries(bench PRIVATE hlp_grep::hlp_grep)

# Optional WFA2-lib backend for the bench executable: run
# tests/scripts/setup_wfa2.sh once to clone and patch third_party/WFA2-lib
# (WFA2's own ctest registration would otherwise fail every `ctest` run).
set(HLP_GREP_WFA2_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../third_party/WFA2-lib)
if(HLP_GREP_WITH_WFA2 AND EXISTS ${HLP_GREP_WFA2_DIR}/CMakeLists.txt)
	add_subdirectory(${HLP_GREP_WFA2_DIR}
	                 ${CMAKE_BINARY_DIR}/third_party/WFA2-lib
	                 EXCLUDE_FROM_ALL)
	target_compile_definitions(bench PRIVATE HLP_GREP_WITH_WFA2)
	target_link_libraries(bench PRIVATE wfa2::wfa2_static)
endif()

# Optional DT-Patricia backend (header-only C++20): run
# tests/scripts/setup_dt_patricia.sh once to clone third_party/DT-Patricia.
set(HLP_GREP_DT_PATRICIA_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../third_party/DT-Patricia)
if(HLP_GREP_WITH_DT_PATRICIA AND EXISTS ${HLP_GREP_DT_PATRICIA_DIR}/CMakeLists.txt)
	add_subdirectory(${HLP_GREP_DT_PATRICIA_DIR}
	                 ${CMAKE_BINARY_DIR}/third_party/DT-Patricia
	                 EXCLUDE_FROM_ALL)
	target_compile_definitions(bench PRIVATE HLP_GREP_WITH_DT_PATRICIA)
	target_link_libraries(bench PRIVATE dt_patricia::dt_patricia)
endif()

# Optional BED-tree backend (legacy C++ without its own CMake build): run
# tests/scripts/setup_bed_tree.sh once to clone and patch
# third_party/bed-tree. main.cpp is excluded; the library sources are
# compiled straight into bench.
	set(HLP_GREP_BED_TREE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../third_party/bed-tree)
if(HLP_GREP_WITH_BED_TREE AND EXISTS ${HLP_GREP_BED_TREE_DIR}/src/SBPTree.cpp)
	target_compile_definitions(bench PRIVATE HLP_GREP_WITH_BED_TREE)
	target_include_directories(bench PRIVATE ${HLP_GREP_BED_TREE_DIR}/src)
	target_sources(bench PRIVATE
		${HLP_GREP_BED_TREE_DIR}/src/BufferManager.cpp
		${HLP_GREP_BED_TREE_DIR}/src/GramDictOrder.cpp
		${HLP_GREP_BED_TREE_DIR}/src/GramLocOrder.cpp
		${HLP_GREP_BED_TREE_DIR}/src/HashFunction.cpp
		${HLP_GREP_BED_TREE_DIR}/src/InvertedList.cpp
		${HLP_GREP_BED_TREE_DIR}/src/Node.cpp
		${HLP_GREP_BED_TREE_DIR}/src/NormGramDictOrder.cpp
		${HLP_GREP_BED_TREE_DIR}/src/Page.cpp
		${HLP_GREP_BED_TREE_DIR}/src/Queue.cpp
		${HLP_GREP_BED_TREE_DIR}/src/SBPTree.cpp
		${HLP_GREP_BED_TREE_DIR}/src/StrComp.cpp
		${HLP_GREP_BED_TREE_DIR}/src/StrDictOrder.cpp
	)
endif()

# Optional hstree backend (exact threshold search from
# TsinghuaDatabaseGroup/Similarity-Search-and-Join, no CMake build of its
# own): run tests/scripts/setup_hstree.sh once to clone and patch
# third_party/Similarity-Search-and-Join. main.cc is excluded; the library
# sources are compiled straight into bench.
set(HLP_GREP_HSTREE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../third_party/Similarity-Search-and-Join/hstree/src)
if(HLP_GREP_WITH_HSTREE AND EXISTS ${HLP_GREP_HSTREE_DIR}/hstree.h)
	target_compile_definitions(bench PRIVATE HLP_GREP_WITH_HSTREE)
	target_include_directories(bench PRIVATE ${HLP_GREP_HSTREE_DIR})
	target_sources(bench PRIVATE
		${HLP_GREP_HSTREE_DIR}/free.cc
		${HLP_GREP_HSTREE_DIR}/hsimSearch.cc
		${HLP_GREP_HSTREE_DIR}/hstopk.cc
		${HLP_GREP_HSTREE_DIR}/io.cc
		${HLP_GREP_HSTREE_DIR}/prepare.cc
	)
endif()

# Optional edlib backend (exact bit-vector alignment library): run
# tests/scripts/setup_edlib.sh once to clone third_party/edlib.
# Examples, utilities and the test tree are switched off; only the
# edlib::edlib library is consumed.
set(HLP_GREP_EDLIB_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../third_party/edlib)
if(HLP_GREP_WITH_EDLIB AND EXISTS ${HLP_GREP_EDLIB_DIR}/CMakeLists.txt)
	set(EDLIB_BUILD_EXAMPLES OFF)
	set(EDLIB_BUILD_UTILITIES OFF)
	set(EDLIB_ENABLE_INSTALL OFF)
	set(BUILD_TESTING OFF)
	add_subdirectory(${HLP_GREP_EDLIB_DIR}
	                 ${CMAKE_BINARY_DIR}/third_party/edlib
	                 EXCLUDE_FROM_ALL)
	target_compile_definitions(bench PRIVATE HLP_GREP_WITH_EDLIB)
	target_link_libraries(bench PRIVATE edlib::edlib)
endif()

# Optional parasail backend (SIMD pairwise alignment library): run
# tests/scripts/setup_parasail.sh once to clone and patch
# third_party/parasail.
set(HLP_GREP_PARASAIL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../third_party/parasail)
if(HLP_GREP_WITH_PARASAIL AND EXISTS ${HLP_GREP_PARASAIL_DIR}/CMakeLists.txt)
	add_subdirectory(${HLP_GREP_PARASAIL_DIR}
	                 ${CMAKE_BINARY_DIR}/third_party/parasail
	                 EXCLUDE_FROM_ALL)
	target_compile_definitions(bench PRIVATE HLP_GREP_WITH_PARASAIL)
	target_link_libraries(bench PRIVATE parasail)
endif()

add_executable(graph_stats graph_stats.cpp)
target_link_libraries(graph_stats PRIVATE hlp_grep::hlp_grep)

add_executable(lifter_stats lifter_stats.cpp)
target_link_libraries(lifter_stats PRIVATE hlp_grep::hlp_grep)

add_executable(visualize_graph visualize_graph.cpp)
target_link_libraries(visualize_graph PRIVATE hlp_grep::hlp_grep)

add_test(NAME hlp_grep_manual COMMAND test_hlp_grep ${CMAKE_CURRENT_SOURCE_DIR}/testcases/manual)
add_test(NAME hlp_grep_generated COMMAND test_hlp_grep ${CMAKE_CURRENT_SOURCE_DIR}/testcases/generated)
add_test(NAME hlp_grep_poa COMMAND poa_test)
add_test(NAME hlp_grep_dist_matrix COMMAND dist_matrix_test)
add_test(NAME hlp_grep_bench COMMAND bench_poa --reps 1 ${CMAKE_CURRENT_SOURCE_DIR}/testcases/manual)
