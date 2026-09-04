cmake_minimum_required(VERSION 3.21)

project(hlp_grep
	VERSION 0.1.0
	DESCRIPTION "Edit distance search queries on DNA sequence dictionary using heavy-light decomposition on pangenome paths"
	LANGUAGES CXX
)

option(HLP_GREP_BUILD_TESTS "Build hlp_grep tests" ${PROJECT_IS_TOP_LEVEL})

if(PROJECT_IS_TOP_LEVEL)
	set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
endif()

add_library(hlp_grep INTERFACE)
add_library(hlp_grep::hlp_grep ALIAS hlp_grep)

target_compile_features(hlp_grep INTERFACE cxx_std_17)

include(GNUInstallDirs)

target_include_directories(hlp_grep INTERFACE
	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
	$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

install(TARGETS hlp_grep
	EXPORT hlp_grepTargets
	ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
	LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
	RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(EXPORT hlp_grepTargets
	FILE hlp_grepTargets.cmake
	NAMESPACE hlp_grep::
	DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/hlp_grep
)

include(CMakePackageConfigHelpers)

configure_package_config_file(
	${CMAKE_CURRENT_SOURCE_DIR}/cmake/hlp_grepConfig.cmake.in
	${CMAKE_CURRENT_BINARY_DIR}/hlp_grepConfig.cmake
	INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/hlp_grep
)

write_basic_package_version_file(
	${CMAKE_CURRENT_BINARY_DIR}/hlp_grepConfigVersion.cmake
	VERSION ${PROJECT_VERSION}
	COMPATIBILITY SameMajorVersion
)

install(FILES
	${CMAKE_CURRENT_BINARY_DIR}/hlp_grepConfig.cmake
	${CMAKE_CURRENT_BINARY_DIR}/hlp_grepConfigVersion.cmake
	DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/hlp_grep
)

if(HLP_GREP_BUILD_TESTS)
	enable_testing()
	add_subdirectory(tests)
endif()
