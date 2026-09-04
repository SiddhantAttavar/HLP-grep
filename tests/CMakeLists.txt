add_executable(test_hlp_grep test_hlp_grep.cpp)
target_link_libraries(test_hlp_grep PRIVATE hlp_grep::hlp_grep)
target_include_directories(test_hlp_grep PRIVATE ${CMAKE_SOURCE_DIR}/src/debug)
add_test(NAME hlp_grep COMMAND test_hlp_grep)
