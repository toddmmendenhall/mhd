# Accumulate sources
set(sources flux_scheme.cpp
            boundary_condition.cpp
            reconstruction.cpp
            solver.cpp)

# Accumulate includes
set(includes context.hpp
             electric_field.hpp
             boundary_condition.hpp
             execution_controller.hpp
             flux_scheme.hpp
             integration.hpp
             kernels.hpp
             reconstruction.hpp
             solver.hpp
             variable_store.hpp)

# Setup library
add_library(solver ${sources} ${includes})

target_link_libraries(solver PUBLIC grid)
target_link_libraries(solver PUBLIC utilities)

target_include_directories(solver PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(solver PUBLIC ${CMAKE_SOURCE_DIR}/include)
