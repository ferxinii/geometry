
#include "points.h"
#include "gtests.h"
#include <stdio.h>

int main() {
    s_point triangle[3] = { {{{0, -1, 0}}},
                            {{{1, 1, 1}}},
                            {{{-1, 1, 1}}} };
    s_point segment1[2] = { {{{0, 0, -2}}},
                            {{{0, 0, 2}}} };
    s_point segment2[2] = { {{{0, 0, 2}}},
                            {{{0, 0, -2}}} };


    printf("Segment crosses triangle: %d (should be 1)\n", test_segment_triangle_intersect_3D(segment1, triangle, 0, 0) == INTERSECT_NONDEGENERATE);
    printf("Segment crosses triangle: %d (should be 1)\n", test_segment_triangle_intersect_3D(segment2, triangle, 0, 0) == INTERSECT_NONDEGENERATE);

    s_point triangle2[3] = { {{{-1, 0, 0}}},
                             {{{1, 0, 0}}},
                             {{{0, 0, 1}}} };
    s_point p1 = {{{0.5, 0, 0}}};
    s_point p2 = {{{2, 0, 0}}};
    printf("point in boundary of triangle: %d (Should be 1)\n", test_point_in_triangle_3D(triangle2, p1, 0, 0) == TEST_BOUNDARY);
    printf("point outside triangle: %d (Should be 1)\n", test_point_in_triangle_3D(triangle2, p2, 0, 0) == TEST_OUT);

}
