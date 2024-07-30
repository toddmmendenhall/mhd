# Setup executable
add_executable(mhd main.cpp)

set_target_properties(mhd PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

target_link_libraries(mhd PUBLIC api)

add_subdirectory(api)
add_subdirectory(grid)
add_subdirectory(solver)
add_subdirectory(utilities)