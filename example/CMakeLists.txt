cmake_minimum_required(VERSION 3.20)
project(geometry_example LANGUAGES C)

add_subdirectory(../ geometry_build)

add_executable(example_convh example_convh.c)

target_link_libraries(example_convh PUBLIC geometry)

target_compile_options(example_convh
    PRIVATE
        -Wall -Werror -Wextra -O3
)


add_custom_command(TARGET example_convh POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E create_symlink
        ${CMAKE_CURRENT_SOURCE_DIR}/lobes
        $<TARGET_FILE_DIR:example_convh>/lobes
    COMMENT "Linking lobes directory into build"
)
