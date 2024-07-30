# Accumulate sources
set(sources solver_impl.cpp)

# Accumulate includes
set(includes solver_impl.hpp)

# Setup library
add_library(solver ${sources} ${includes})

target_include_directories(solver PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(solver PRIVATE api)