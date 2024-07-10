# Accumulate sources
set(sources calc.cpp
            error.cpp
            exception.cpp
            grid.cpp
            logger.cpp
            profile.cpp
            safe.cpp)

# Accumulate private includes
set(include_dir ${CMAKE_CURRENT_SOURCE_DIR}/include)

set(includes ${private_include_dir}/exception.hpp
             ${private_include_dir}/logger.hpp
             ${private_include_dir}/safe.hpp)

# Accumulate public includes
set(public_include_dir ${CMAKE_SOURCE_DIR}/include)

set(includes ${public_include_dir}/calc.hpp
             ${public_include_dir}/error.hpp
             ${public_include_dir}/grid.hpp
             ${public_include_dir}/profile.hpp)

# Setup library
add_library(${CMAKE_PROJECT_NAME} ${sources} ${includes} ${public_includes})

# Add include directories
target_include_directories(${CMAKE_PROJECT_NAME} PRIVATE ${include_dir})
target_include_directories(${CMAKE_PROJECT_NAME} PUBLIC ${public_include_dir})

