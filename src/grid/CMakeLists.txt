# Accumulate sources
set(sources cartesian.cpp
            cell.cpp
            grid.cpp
            grid_factory.cpp)

# Accumulate includes
set(includes ${CMAKE_CURRENT_SOURCE_DIR}/include/cartesian.hpp
             ${CMAKE_CURRENT_SOURCE_DIR}/include/cell.hpp
             ${CMAKE_CURRENT_SOURCE_DIR}/include/grid.hpp
             ${CMAKE_CURRENT_SOURCE_DIR}/include/grid_factory.hpp)

# Setup library
add_library(grid ${sources} ${includes})

target_include_directories(grid PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(grid PUBLIC api)
