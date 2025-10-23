cmake_minimum_required(VERSION 3.20)
project(geometry LANGUAGES C)


add_library(geometry STATIC
    src/geometry.c
    src/convh.c

    # Include the predicates .c files directly
    external/predicates/src/constants.c
    external/predicates/src/predicates.c
    external/predicates/src/random.c
    external/predicates/src/printing.c
)

target_include_directories(geometry
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/external/predicates/include
)

target_compile_options(geometry
    PRIVATE
        -Wall -Werror -Wextra -O3
)

if (NOT TARGET convhull_3d)
    add_subdirectory(external/convhull_3d)
endif()

target_link_libraries(geometry
    PRIVATE
        convhull_3d
)

