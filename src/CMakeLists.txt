set(sources api_safe.cpp
            calc.cpp
            grid.cpp
            profile.cpp)

set(include_dir ${CMAKE_CURRENT_SOURCE_DIR}/include)
set(includes ${include_dir}/api_safe.hpp
             ${include_dir}/error.hpp
             ${include_dir}/exception.hpp
             ${include_dir}/logger.hpp)

set(public_include_dir ${CMAKE_SOURCE_DIR}/include)
set(public_includes ${public_include_dir}/calc.hpp
                    ${public_include_dir}/grid.hpp
                    ${public_include_dir}/profile.hpp)

set(library_name mhd)
add_library(${library_name} STATIC ${sources} ${includes} ${public_includes})

target_include_directories(${library_name} PRIVATE ${include_dir})
target_include_directories(${library_name} PUBLIC ${public_include_dir})
