# Accumulate sources
set(sources 1d.cpp)

# Accumulate includes
set(includes 1d.hpp)

# Setup library
add_library(cartesian ${sources} ${includes})

target_include_directories(cartesian PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(cartesian PRIVATE ${CMAKE_SOURCE_DIR}/include)

target_link_libraries(cartesian PUBLIC grid)
target_link_libraries(cartesian PUBLIC utilities)