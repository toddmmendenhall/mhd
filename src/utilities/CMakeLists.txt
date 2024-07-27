# Accumulate sources
set(sources point.cpp
            vector.cpp)

# Accumulate includes
set(includes constants.hpp
             point.hpp
             vector.hpp)

# Setup library
add_library(utilities ${sources} ${includes})

target_include_directories(utilities PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
