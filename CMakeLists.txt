cmake_minimum_required(VERSION 3.22)

project(mhd)

set(CMAKE_CXX_FLAGS "-fPIC")

add_subdirectory(src)

add_subdirectory(external/googletest)
enable_testing()
include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

add_subdirectory(test)