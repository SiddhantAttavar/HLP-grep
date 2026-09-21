add_executable(solve naive_solver.cpp)
target_link_libraries(solve PRIVATE hlp_grep::hlp_grep)

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

add_executable(graph_stats graph_stats.cpp)
target_link_libraries(graph_stats PRIVATE hlp_grep::hlp_grep)

add_executable(lifter_stats lifter_stats.cpp)
target_link_libraries(lifter_stats PRIVATE hlp_grep::hlp_grep)

add_test(NAME hlp_grep_manual COMMAND test_hlp_grep ${CMAKE_CURRENT_SOURCE_DIR}/testcases/manual)
add_test(NAME hlp_grep_generated COMMAND test_hlp_grep ${CMAKE_CURRENT_SOURCE_DIR}/testcases/generated)
add_test(NAME hlp_grep_poa COMMAND poa_test)
add_test(NAME hlp_grep_dist_matrix COMMAND dist_matrix_test)
add_test(NAME hlp_grep_bench COMMAND bench_poa --reps 1 ${CMAKE_CURRENT_SOURCE_DIR}/testcases/manual)
