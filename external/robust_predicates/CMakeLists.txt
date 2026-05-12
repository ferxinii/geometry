cmake_minimum_required(VERSION 3.13)
project(robust_predicates LANGUAGES CXX)

add_library(robust_predicates STATIC
    robust_predicates_n1.cpp
    robust_predicates_n2.cpp
    robust_predicates_n3.cpp
)

target_include_directories(robust_predicates PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/includes/geometry
    ${CMAKE_CURRENT_SOURCE_DIR}/includes/mp11
)

target_include_directories(robust_predicates PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_compile_features(robust_predicates PRIVATE cxx_std_20)


# Define a flag: pass -DROBUST_PREDICATES_DEBUG=ON to cmake to opt into debug mode.
option(ROBUST_PREDICATES_DEBUG "Build robust_predicates with debug info instead of -O3" OFF)

if (NOT ROBUST_PREDICATES_DEBUG)
    # Always compile at -O3 regardless of CMAKE_BUILD_TYPE.
    # -g0 explicitly disables debug info even when Debug flags are injected.
    target_compile_options(robust_predicates PRIVATE
        -O3 -march=native -DNDEBUG -g0
    )
else()
    # Debug mode: strip optimizations, enable debug symbols, keep assertions.
    target_compile_options(robust_predicates PRIVATE
        -O0 -g
    )
    message(STATUS "robust_predicates: DEBUG mode enabled (-O0 -g)")
endif()


# Force any target linking against this to use the C++ linker
target_link_libraries(robust_predicates INTERFACE stdc++)

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

