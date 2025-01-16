# Accumulate sources
set(sources calc.cpp
            profile.cpp)

# Accumulate includes
set(includes exception.hpp
             logger.hpp)

set(public_include_dir ${CMAKE_SOURCE_DIR}/include)
set(public_includes ${public_include_dir}/calc.hpp
                    ${public_include_dir}/error.hpp
                    ${public_include_dir}/profile.hpp
                    ${public_include_dir}/profile_options.hpp)

# Setup library
add_library(api SHARED ${sources} ${includes} ${public_includes})

target_include_directories(api PUBLIC ${public_include_dir})
target_include_directories(api PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(api PRIVATE grid)
target_link_libraries(api PRIVATE solver)