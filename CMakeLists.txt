cmake_minimum_required(VERSION 3.20)
project(geometry LANGUAGES C)


add_library(geometry STATIC
    src/points.c
    src/gtests.c
    src/convh.c
    src/ch_quickhull3D.c
    src/ch_intersect.c
)

target_include_directories(geometry
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
    # PRIVATE
        # ${CMAKE_CURRENT_SOURCE_DIR}/external/predicates/include
)

target_compile_options(geometry
    PRIVATE
        -Wall -Werror -Wextra
)


if (NOT TARGET predicates)
    add_subdirectory(external/predicates_shewchuck)
endif()


target_link_libraries(geometry
    PRIVATE
        predicates
)

