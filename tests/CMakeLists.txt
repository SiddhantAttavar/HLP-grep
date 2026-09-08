add_executable(solve naive_solver.cpp)
target_link_libraries(solve PRIVATE hlp_grep::hlp_grep)

add_executable(test_hlp_grep test.cpp)
target_link_libraries(test_hlp_grep PRIVATE hlp_grep::hlp_grep)

add_test(NAME hlp_grep_manual COMMAND test_hlp_grep ${CMAKE_CURRENT_SOURCE_DIR}/testcases/manual)
add_test(NAME hlp_grep_generated COMMAND test_hlp_grep ${CMAKE_CURRENT_SOURCE_DIR}/testcases/generated)
