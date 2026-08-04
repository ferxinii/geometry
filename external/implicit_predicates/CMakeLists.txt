cmake_minimum_required(VERSION 3.14)
project(implicit_predicates LANGUAGES CXX)

add_library(implicit_predicates STATIC
    WRAPPER.cpp
)

# Make geometry's points.h visible
target_include_directories(implicit_predicates
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/../../include
)

target_compile_features(implicit_predicates PRIVATE cxx_std_20)

target_compile_options(implicit_predicates PRIVATE
    $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-frounding-math>
    $<$<CXX_COMPILER_ID:MSVC>:/fp:strict>
)

