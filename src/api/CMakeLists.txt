# Accumulate sources
set(sources calc.cpp
            grid.cpp
            profile.cpp)

# Set public include dir
set(include_dir ${CMAKE_SOURCE_DIR}/include)

# Accumulate public includes
set(includes ${include_dir}/calc.hpp
             ${include_dir}/grid.hpp
             ${include_dir}/profile.hpp)

# Set private include dir
set(private_include_dir ${CMAKE_CURRENT_SOURCE_DIR}/include)

# Accumulate private includes
set(private_includes ${private_include_dir}/exception.hpp
                     ${private_include_dir}/logger.hpp
                     ${private_include_dir}/safe.hpp)

# Setup library
add_library(api ${sources} ${includes} ${private_includes})

# Add include directories
target_include_directories(api PUBLIC ${include_dir})
target_include_directories(api PRIVATE ${private_include_dir})
