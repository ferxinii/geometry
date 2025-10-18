cmake_minimum_required(VERSION 3.20)
project(geometry LANGUAGES C)

add_library(geometry STATIC
    geometry.c

    # Include the predicates .c files directly
    predicates/src/constants.c
    predicates/src/predicates.c
    predicates/src/random.c
    predicates/src/printing.c
)

target_include_directories(geometry
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/predicates/include
        ${CMAKE_CURRENT_SOURCE_DIR}/predicates/src
)

target_compile_options(geometry
    PRIVATE
        -Wall -Werror -Wextra -O3
)
