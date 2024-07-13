cmake_minimum_required(VERSION 3.22)

project(mhd)

add_subdirectory(src)

add_subdirectory(external/googletest)
enable_testing()
include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

add_executable(grid_tests test/grid_tests.cpp)
target_link_libraries(grid_tests GTest::gtest GTest::gtest_main)
add_test(GridTests grid_tests)
