# Accumulate sources
set(sources cartesian.cpp)

# Accumulate private includes
set(include_dir ${CMAKE_CURRENT_SOURCE_DIR}/include)

set(includes ${include_dir}/cartesian.hpp)

# Accumulate public includes
set(public_include_dir ${CMAKE_SOURCE_DIR}/include)

set(public_includes ${public_include_dir}/grid.hpp)

# Setup library
add_library(grid ${sources} ${includes} ${public_includes})

# Add include directories
#target_include_directories(api PRIVATE ${include_dir})
#target_include_directories(api PUBLIC ${public_include_dir})

target_link_libraries(grid ${CMAKE_PROJECT_NAME})