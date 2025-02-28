# Accumulate sources
set(sources solver.cpp)

set(reconstruction_sources reconstruction/reconstruction.hpp
                           reconstruction/reconstruction.cpp)

set(flux_sources flux/flux_scheme.hpp
                 flux/flux_scheme.cpp)

set(boundary_condition_sources boundary_condition/boundary_condition.hpp
                               boundary_condition/boundary_condition.cpp)

# Accumulate includes
set(includes context.hpp
             execution_controller.hpp
             integration.hpp
             kernels.hpp
             residual.hpp
             solver.hpp
             variable_store.hpp)

# Setup library
add_library(solver ${sources} ${reconstruction_sources} ${flux_sources} ${boundary_condition_sources} ${includes})

target_link_libraries(solver PUBLIC grid)
target_link_libraries(solver PUBLIC utilities)

target_include_directories(solver PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(solver PUBLIC ${CMAKE_SOURCE_DIR}/include)
