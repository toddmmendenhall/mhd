# Accumulate sources
set(sources 1d.cpp
            2d.cpp
            3d.cpp)

# Accumulate includes
set(includes 1d.hpp
             2d.hpp
             3d.hpp)

# Setup library
add_library(cartesian ${sources} ${includes})

target_include_directories(cartesian PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(cartesian PRIVATE grid)
target_link_libraries(cartesian PUBLIC utilities)