# Accumulate sources
set(sources calc.cpp
            grid.cpp
            profile.cpp)

# Accumulate includes
set(includes exception.hpp
             logger.hpp)

set(public_includes ${CMAKE_SOURCE_DIR}/include/calc.hpp
                    ${CMAKE_SOURCE_DIR}/include/error.hpp
                    ${CMAKE_SOURCE_DIR}/include/grid.hpp
                    ${CMAKE_SOURCE_DIR}/include/profile.hpp
                    ${CMAKE_SOURCE_DIR}/include/profile_options.hpp)

# Setup library
add_library(api ${sources} ${includes} ${public_includes})

target_include_directories(api PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(api PUBLIC ${CMAKE_SOURCE_DIR}/include)

target_link_libraries(api PUBLIC grid)
