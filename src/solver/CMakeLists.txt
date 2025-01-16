# Accumulate sources
set(sources solver.cpp
            thermo_data.cpp)

# Accumulate includes
set(includes kernels.hpp
             solver.hpp
             state.hpp
             thermo_data.hpp)

# Setup library
add_library(solver ${sources} ${includes})

target_include_directories(solver PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(solver PUBLIC ${CMAKE_SOURCE_DIR}/include)

target_link_libraries(solver PRIVATE grid)