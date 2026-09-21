add_executable(solve naive_solver.cpp)
target_link_libraries(solve PRIVATE hlp_grep::hlp_grep)

option(HLP_GREP_WITH_WFA2 "Enable the WFA2-lib method in bench (needs third_party/WFA2-lib)" ON)
option(HLP_GREP_WITH_DT_PATRICIA "Enable the DT-Patricia method in bench (needs third_party/DT-Patricia)" ON)

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

add_executable(graph_stats graph_stats.cpp)
target_link_libraries(graph_stats PRIVATE hlp_grep::hlp_grep)

add_executable(lifter_stats lifter_stats.cpp)
target_link_libraries(lifter_stats PRIVATE hlp_grep::hlp_grep)

add_test(NAME hlp_grep_manual COMMAND test_hlp_grep ${CMAKE_CURRENT_SOURCE_DIR}/testcases/manual)
add_test(NAME hlp_grep_generated COMMAND test_hlp_grep ${CMAKE_CURRENT_SOURCE_DIR}/testcases/generated)
add_test(NAME hlp_grep_poa COMMAND poa_test)
add_test(NAME hlp_grep_dist_matrix COMMAND dist_matrix_test)
add_test(NAME hlp_grep_bench COMMAND bench_poa --reps 1 ${CMAKE_CURRENT_SOURCE_DIR}/testcases/manual)
