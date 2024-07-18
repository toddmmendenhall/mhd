add_executable(grid_tests grid_tests.cpp)
target_link_libraries(grid_tests GTest::gtest GTest::gtest_main)
add_test(GridTests grid_tests)
