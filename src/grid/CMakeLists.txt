# Accumulate sources
set(sources grid_factory.cpp)

# Accumulate includes
set(includes grid_factory.hpp)

# Setup library
add_library(grid ${sources} ${includes})

target_include_directories(grid PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_subdirectory(cartesian)

target_link_libraries(grid PUBLIC api)
target_link_libraries(grid PRIVATE cartesian)
