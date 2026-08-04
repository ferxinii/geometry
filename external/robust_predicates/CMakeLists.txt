cmake_minimum_required(VERSION 3.13)
project(robust_predicates LANGUAGES CXX)

add_library(predicates_basic STATIC
    predicates_basic_n2.cpp
    predicates_basic_n3.cpp
    predicates_basic_other.cpp
)

add_library(predicates_power EXCLUDE_FROM_ALL STATIC
    predicates_power_n1.cpp
    predicates_power_n2.cpp
    predicates_power_n3.cpp
)

add_library(predicates_lp EXCLUDE_FROM_ALL STATIC
    predicates_lp_n2.cpp
    predicates_lp_n3.cpp
)

foreach(_tgt predicates_basic predicates_power predicates_lp)
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
    foreach(_tgt predicates_basic predicates_power predicates_lp)
        target_compile_options(${_tgt} PRIVATE
            -O3 -march=native -DNDEBUG -g0
        )
    endforeach()
else()
    foreach(_tgt predicates_basic predicates_power predicates_lp)
        target_compile_options(${_tgt} PRIVATE
            -O0 -g
        )
    endforeach()
    message(STATUS "robust_predicates: DEBUG mode enabled (-O0 -g)")
endif()

option(ROBUST_PREDICATES_PRINT_SIZE "Emit compile-time deprecation warnings showing stage_d buffer sizes" OFF)
if(ROBUST_PREDICATES_PRINT_SIZE)
    foreach(_tgt predicates_basic predicates_power predicates_lp)
        target_compile_definitions(${_tgt} PRIVATE ROBUST_PREDICATES_PRINT_SIZE)
    endforeach()
    message(STATUS "robust_predicates: ROBUST_PREDICATES_PRINT_SIZE enabled — expect deprecation warnings during build")
endif()

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
