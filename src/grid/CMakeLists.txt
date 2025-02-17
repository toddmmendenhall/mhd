# Accumulate sources
set(sources grid.cpp)

# Accumulate includes
set(includes grid.hpp)

# Setup library
add_library(grid ${sources} ${includes})

target_include_directories(grid PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(grid PRIVATE ${CMAKE_SOURCE_DIR}/include)

add_subdirectory(cartesian)

target_link_libraries(grid PRIVATE cartesian)