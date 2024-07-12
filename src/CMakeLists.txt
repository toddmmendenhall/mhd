# Accumulate sources
set(sources calc.cpp
            error.cpp
            exception.cpp
            grid.cpp
            logger.cpp
            profile.cpp)

# Accumulate includes
set(include_dir ${CMAKE_CURRENT_SOURCE_DIR}/include)

set(includes ${include_dir}/exception.hpp
             ${include_dir}/logger.hpp
             ${include_dir}/calc.hpp
             ${include_dir}/error.hpp
             ${include_dir}/grid.hpp
             ${include_dir}/profile.hpp)

# Setup library
add_library(${CMAKE_PROJECT_NAME} ${sources} ${includes})

# Add include directories
target_include_directories(${CMAKE_PROJECT_NAME} PRIVATE ${include_dir})
