set(library_name cartesian)

# Accumulate sources
set(sources 1d.cpp
            2d.cpp
            3d.cpp)

# Accumulate includes
set(includes 1d.hpp
             2d.hpp
             3d.hpp)

# Setup library
add_library(${library_name} STATIC ${sources} ${includes})

target_include_directories(${library_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(${library_name} PRIVATE ${CMAKE_SOURCE_DIR}/include)

target_link_libraries(${library_name} PRIVATE grid)
target_link_libraries(${library_name} PRIVATE utilities)