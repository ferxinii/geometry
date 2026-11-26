cmake_minimum_required(VERSION 3.20)
project(geometry LANGUAGES C)

# Default build (no CMAKE_BUILD_TYPE): -O3 + link-time optimization, asserts kept
# (they are load-bearing). An explicit CMAKE_BUILD_TYPE (e.g. Debug) overrides it.
if (NOT CMAKE_BUILD_TYPE)
    add_compile_options(-O3)
    include(CheckIPOSupported)
    check_ipo_supported(RESULT _ipo_ok OUTPUT _ipo_msg LANGUAGES C)
    if (_ipo_ok)
        set(CMAKE_INTERPROCEDURAL_OPTIMIZATION ON)
    else()
        message(STATUS "LTO not supported, building without it: ${_ipo_msg}")
    endif()
endif()

option(GEOMETRY_COMPILE_TESTS_POWER "Enable the power predicates" OFF)
option(GEOMETRY_COMPILE_LP "Enable the Linear Programming functions" OFF)
option(GEOMETRY_COMPILE_TESTS_IMPLICIT "Enable the implicit (indirect) predicates" OFF)

# Weighted LP (lp3_* with non-NULL weights) uses test_orthosphere, which only
# exists with the power predicates: LP therefore turns them on.
if(GEOMETRY_COMPILE_LP)
    set(GEOMETRY_COMPILE_TESTS_POWER ON CACHE BOOL "" FORCE)
endif()

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

