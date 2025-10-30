cmake_minimum_required(VERSION 3.13)
project(predicates_shewchuck LANGUAGES C)

add_library(predicates STATIC
    predicates.c
)

set_target_properties(predicates PROPERTIES
    C_STANDARD 99
    C_STANDARD_REQUIRED ON
    POSITION_INDEPENDENT_CODE ON
)


target_compile_options(predicates
    PRIVATE
        -Wno-deprecated-non-prototype -O3 -fno-fast-math
)

