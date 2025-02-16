add_executable(mhd_tests api_tests.cpp)

target_link_libraries(mhd_tests GTest::gtest GTest::gtest_main)
target_link_libraries(mhd_tests api)
# target_link_libraries(mhd_tests grid)
# target_link_libraries(mhd_tests solver)
# target_link_libraries(mhd_tests utilities)

add_test(APITests mhd_tests)
# add_test(GridTests mhd_tests)
# add_test(SolverTests mhd_tests)
