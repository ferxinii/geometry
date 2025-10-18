cmake_minimum_required(VERSION 3.20)
project(geometry LANGUAGE C)

add_subdirectory(predicates)

add_library(geometry STATIC
    geometry.c
)

target_include_directories(geometry
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(geometry
    PRIVATE 
        predicates
)


target_compile_options(geometry
    PRIVATE
        -Wall -Werror -Wextra -O3
)
