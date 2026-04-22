cmake_minimum_required(VERSION 3.13)
project(robust_predicates LANGUAGES CXX)

add_library(robust_predicates STATIC
    robust_predicates.cpp
)

target_include_directories(robust_predicates PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/includes/geometry
    ${CMAKE_CURRENT_SOURCE_DIR}/includes/mp11
)

target_include_directories(robust_predicates PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_compile_definitions(robust_predicates PRIVATE NDEBUG)
target_compile_features(robust_predicates PRIVATE cxx_std_20)
target_compile_options(robust_predicates PRIVATE -O3 -march=native)

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

