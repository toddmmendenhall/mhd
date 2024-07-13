# Accumulate sources
set(sources calc.cpp
            error.cpp
            exception.cpp
            logger.cpp
            profile.cpp)

# Accumulate includes
set(includes ${CMAKE_CURRENT_SOURCE_DIR}/include/exception.hpp
             ${CMAKE_CURRENT_SOURCE_DIR}/include/logger.hpp)

set(public_includes ${CMAKE_SOURCE_DIR}/include/calc.hpp
                    ${CMAKE_SOURCE_DIR}/include/error.hpp
                    ${CMAKE_SOURCE_DIR}/include/profile_options.hpp
                    ${CMAKE_SOURCE_DIR}/include/profile.hpp)

# Setup library
add_library(api ${sources} ${includes} ${public_includes})

target_include_directories(api PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories(api PUBLIC ${CMAKE_SOURCE_DIR}/include)

target_link_libraries(api PUBLIC grid)