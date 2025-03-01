# Accumulate sources
set(sources grid.hpp grid.cpp)

set(1d_sources 1d.hpp 1d.cpp)

# Setup library
add_library(grid ${sources} ${1d_sources} ${includes})

target_include_directories(grid PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(grid PRIVATE ${CMAKE_SOURCE_DIR}/include)
