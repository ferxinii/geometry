cmake_minimum_required(VERSION 3.20)
project(geometry_example LANGUAGES C)

add_subdirectory(../ geometry_build)

# example_gtests:
add_executable(example_gtests example_gtests.c)
target_link_libraries(example_gtests PUBLIC geometry)
target_compile_options(example_gtests
    PRIVATE
        -Wall -Werror -Wextra -O3
)


# example_convh:
add_executable(example_convh example_convh.c)
target_link_libraries(example_convh PUBLIC geometry)
target_compile_options(example_convh
    PRIVATE
        -Wall -Werror -Wextra -O3
)


# example_intersect:
add_executable(example_intersect example_intersect.c)
target_link_libraries(example_intersect PUBLIC geometry)
target_compile_options(example_intersect
    PRIVATE
        -Wall -Werror -Wextra -O3
)



