# Accumulate sources
set(sources calc.cpp
            error.cpp
            exception.cpp
            logger.cpp
            profile.cpp)

# Accumulate includes
set(includes ${CMAKE_CURRENT_SOURCE_DIR}/include/calc.hpp
             ${CMAKE_CURRENT_SOURCE_DIR}/include/error.hpp
             ${CMAKE_CURRENT_SOURCE_DIR}/include/exception.hpp
             ${CMAKE_CURRENT_SOURCE_DIR}/include/logger.hpp
             ${CMAKE_CURRENT_SOURCE_DIR}/include/profile_options.hpp
             ${CMAKE_CURRENT_SOURCE_DIR}/include/profile.hpp)

# Setup library
add_library(mhd ${sources} ${includes})

target_include_directories(mhd PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

add_subdirectory(grid)

target_link_libraries(mhd PUBLIC grid)