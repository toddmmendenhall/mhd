set(library_name grid)

# Accumulate sources
set(sources grid.cpp)

# Accumulate includes
set(includes grid.hpp)

# Setup library
add_library(${library_name} STATIC ${sources} ${includes})

target_include_directories(${library_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(${library_name} PRIVATE ${CMAKE_SOURCE_DIR}/include)

add_subdirectory(cartesian)

target_link_libraries(${library_name} PRIVATE cartesian)