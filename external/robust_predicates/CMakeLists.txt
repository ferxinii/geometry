cmake_minimum_required(VERSION 3.13)
project(robust_predicates LANGUAGES CXX)

add_library(robust_predicates STATIC
    robust_predicates_n1.cpp
    robust_predicates_n2.cpp
    robust_predicates_n3.cpp
)

add_library(voronoi_predicates EXCLUDE_FROM_ALL STATIC
    voronoi_predicates.cpp
)

find_library(MPFR_LIBRARY
    NAMES mpfr
    HINTS /opt/homebrew/lib /usr/local/lib
)
find_path(MPFR_INCLUDE_DIR
    NAMES mpfr.h
    HINTS /opt/homebrew/include /usr/local/include
)
if(NOT MPFR_LIBRARY)
    message(FATAL_ERROR "MPFR library not found. Install with: brew install mpfr")
endif()
if(NOT MPFR_INCLUDE_DIR)
    message(FATAL_ERROR "MPFR headers not found. Install with: brew install mpfr")
endif()
target_include_directories(voronoi_predicates PRIVATE ${MPFR_INCLUDE_DIR})
target_link_libraries(voronoi_predicates PRIVATE ${MPFR_LIBRARY})

foreach(_tgt robust_predicates voronoi_predicates)
    target_include_directories(${_tgt} PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/includes/geometry
        ${CMAKE_CURRENT_SOURCE_DIR}/includes/mp11
    )
    target_include_directories(${_tgt} PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
    )
    target_compile_features(${_tgt} PRIVATE cxx_std_20)
    target_link_libraries(${_tgt} INTERFACE stdc++)
endforeach()

# Define a flag: pass -DROBUST_PREDICATES_DEBUG=ON to cmake to opt into debug mode.
option(ROBUST_PREDICATES_DEBUG "Build robust_predicates with debug info instead of -O3" OFF)

if (NOT ROBUST_PREDICATES_DEBUG)
    # Always compile at -O3 regardless of CMAKE_BUILD_TYPE.
    # -g0 explicitly disables debug info even when Debug flags are injected.
    foreach(_tgt robust_predicates voronoi_predicates)
        target_compile_options(${_tgt} PRIVATE
            -O3 -march=native -DNDEBUG -g0
        )
    endforeach()
else()
    # Debug mode: strip optimizations, enable debug symbols, keep assertions.
    foreach(_tgt robust_predicates voronoi_predicates)
        target_compile_options(${_tgt} PRIVATE
            -O0 -g
        )
    endforeach()
    message(STATUS "robust_predicates: DEBUG mode enabled (-O0 -g)")
endif()

option(ROBUST_PREDICATES_PRINT_SIZE "Emit compile-time deprecation warnings showing stage_d buffer sizes" OFF)
if(ROBUST_PREDICATES_PRINT_SIZE)
    foreach(_tgt robust_predicates voronoi_predicates)
        target_compile_definitions(${_tgt} PRIVATE ROBUST_PREDICATES_PRINT_SIZE)
    endforeach()
    message(STATUS "robust_predicates: ROBUST_PREDICATES_PRINT_SIZE enabled — expect deprecation warnings during build")
endif()

# Optionally build the executable to print expressions
option(BUILD_PRINT_EXPRESSIONS "Build the print_expressions utility" OFF)
if(BUILD_PRINT_EXPRESSIONS)
    add_executable(print_expressions print_error_bounds.cpp)
    target_include_directories(print_expressions PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/includes/geometry
        ${CMAKE_CURRENT_SOURCE_DIR}/includes/mp11
    )
    target_compile_features(print_expressions PRIVATE cxx_std_20)
    target_compile_options(print_expressions PRIVATE -O3)
endif()

