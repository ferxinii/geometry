cmake_minimum_required(VERSION 3.20)
project(geometry LANGUAGES C)

option(GEOMETRY_COMPILE_TESTS_POWER "Enable the power predicates" OFF)
option(GEOMETRY_COMPILE_LP "Enable the Linear Programming functions" OFF)
option(GEOMETRY_COMPILE_TESTS_IMPLICIT "Enable the implicit (indirect) predicates" OFF)

add_library(geometry STATIC
    src/points.c
    src/gtests.c
    src/convh.c
    src/ch_quickhull3D.c
    src/ch_intersect.c
    src/trimesh.c
    src/trimesh_repair.c
)

target_include_directories(geometry
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_compile_options(geometry
    PRIVATE
        -Wall -Werror -Wextra
)

if (NOT TARGET robust_predicates)
    add_subdirectory(external/robust_predicates)
endif()

target_link_libraries(geometry
    PRIVATE
        predicates_basic
)


if(GEOMETRY_COMPILE_TESTS_POWER)
    target_link_libraries(geometry PRIVATE predicates_power)
    target_compile_definitions(geometry PUBLIC GEOMETRY_COMPILE_TESTS_POWER)
endif()

if(GEOMETRY_COMPILE_LP)
    target_link_libraries(geometry PRIVATE predicates_lp)
    target_compile_definitions(geometry PUBLIC GEOMETRY_COMPILE_LP)
    target_sources(geometry PRIVATE src/lp.c)
endif()

if(GEOMETRY_COMPILE_TESTS_IMPLICIT)
    if (NOT TARGET implicit_predicates)
        add_subdirectory(external/implicit_predicates)
    endif()
    target_link_libraries(geometry PRIVATE implicit_predicates)
    target_compile_definitions(geometry PUBLIC GEOMETRY_COMPILE_TESTS_IMPLICIT)
endif()

