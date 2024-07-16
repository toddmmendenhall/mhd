# Accumulate sources
set(sources cartesian.cpp
            cell.cpp
            grid.cpp
            grid_factory.cpp)

# Accumulate includes
set(includes cartesian.hpp
             cell.hpp
             grid.hpp
             grid_factory.hpp)

# Setup library
add_library(grid ${sources} ${includes})

target_include_directories(grid PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(grid PUBLIC api)
