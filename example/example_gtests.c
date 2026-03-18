
#include "points.h"
#include "gtests.h"
#include <stdio.h>
#include <math.h>

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



    // INSPHERE_WEIGHTED
    printf("\n\n\n");

    // All weights 0:
    double pa[] = {1/sqrt(2), 0, -1/sqrt(2)};
    double pb[] = {-1/sqrt(2), 0, -1/sqrt(2)};
    double pc[] = {0, 1/sqrt(2), 1/sqrt(2)};
    double pd[] = {0, -1/sqrt(2), 1/sqrt(2)};
    double w[] = {0, 0, 0, 0};
    double pe_in[]  = {0, 0, 0, 0};   // origin — inside circumsphere
    double pe_out[] = {10, 0, 0, 0};  // far away — outside
    
    printf("Test wi=0\n");
    printf("inside: %g (regular insphere: %g)\n", insphere_weighted(pa, w[0], pb, w[1], pc, w[2], pd, w[3], pe_in, 0), insphere(pa,pb,pc,pd,pe_in));
    printf("outside: %g (regular insphere: %g)\n", insphere_weighted(pa, w[0], pb, w[1], pc, w[2], pd, w[3], pe_out, 0), insphere(pa,pb,pc,pd,pe_out));

    
    // pe exactly on the circumsphere, and then change weight:
    double pe_0[] = {1.0, 0.0, 0.0, 0.0};
    double pe_in1[]  = {1.0, 0.0, 0.0,  1.0};  // w = +R² = 1 : inside
    double pe_out1[] = {1.0, 0.0, 0.0, -1.0};  // w = -R² = -1 : outside
    printf("Test pe exactly on circumsphere, and then change weight\n");
    printf("exactly (w=0): %g (regular insphere: %g)\n", 
                insphere_weighted(pa, w[0], pb, w[1], pc, w[2], pd, w[3], pe_0, 0), insphere(pa,pb,pc,pd,pe_0));
    printf("exactly (w=1): %g (regular insphere: %g)\n",
                insphere_weighted(pa, w[0], pb, w[1], pc, w[2], pd, w[3], pe_in1, 0), insphere(pa,pb,pc,pd,pe_in1));
    printf("exactly (w=-1): %g (regular insphere: %g)\n",
                insphere_weighted(pa, w[0], pb, w[1], pc, w[2], pd, w[3], pe_out1, 0), insphere(pa,pb,pc,pd,pe_out1));


}
