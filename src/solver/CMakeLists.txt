# Accumulate sources
set(sources solver.cpp
            )

# Accumulate includes
set(includes execution_controller.hpp
             kernels.hpp
             solver.hpp
             variable_store.hpp)

# Setup library
add_library(solver ${sources} ${includes})

target_include_directories(solver PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(solver PUBLIC ${CMAKE_SOURCE_DIR}/include)
