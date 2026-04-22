
#include "points.h"
#include "gtests.h"
#include <stdio.h>
#include <math.h>

int main() 
{
    // ==========================================================================
    // ORIENT2D
    // ==========================================================================
    {
        printf("\n=== orient2d ===\n");
        // counterclockwise triangle
        double a[2] = {0, 0};
        double b[2] = {1, 0};
        double c[2] = {0, 1};
        printf("ccw:        %d (should be  1)\n", orient2d(a[0],a[1], b[0],b[1], c[0],c[1]));
        printf("cw:         %d (should be -1)\n", orient2d(a[0],a[1], c[0],c[1], b[0],b[1]));
        printf("collinear:  %d (should be  0)\n", orient2d(a[0],a[1], b[0],b[1], b[0],b[1]));
    }

    // ==========================================================================
    // ORIENT3D
    // ==========================================================================
    {
        printf("\n=== orient3d ===\n");
        // positive orientation tetrahedron
        s_point a = {{{0, 0, 0}}};
        s_point b = {{{1, 0, 0}}};
        s_point c = {{{0, 1, 0}}};
        s_point d = {{{0, 0, 1}}};  // same side as normal: negative
        s_point e = {{{0, 0,-1}}};  // opposite side: positive
        printf("same side as plane normal:   %d (should be  -1)\n",
            orient3d(a.x,a.y,a.z, b.x,b.y,b.z, c.x,c.y,c.z, d.x,d.y,d.z));
        printf("opposite side as plane normal:   %d (should be +1)\n",
            orient3d(a.x,a.y,a.z, b.x,b.y,b.z, c.x,c.y,c.z, e.x,e.y,e.z));
        printf("coplanar:   %d (should be  0)\n",
            orient3d(a.x,a.y,a.z, b.x,b.y,b.z, c.x,c.y,c.z, a.x,a.y,a.z));
    }

    // ==========================================================================
    // INCIRCLE
    // ==========================================================================
    {
        printf("\n=== incircle ===\n");
        // unit circle: a=(1,0), b=(0,1), c=(-1,0)
        double a[2] = { 1,  0};
        double b[2] = { 0,  1};
        double c[2] = {-1,  0};
        double d_on[2]  = { 0, -1};  // on circle
        double d_in[2]  = { 0,  0};  // inside
        double d_out[2] = { 2,  0};  // outside
        printf("on circle:  %d (should be  0)\n",
            incircle(a[0],a[1], b[0],b[1], c[0],c[1], d_on[0],d_on[1]));
        printf("inside:     %d (should be  1)\n",
            incircle(a[0],a[1], b[0],b[1], c[0],c[1], d_in[0],d_in[1]));
        printf("outside:    %d (should be -1)\n",
            incircle(a[0],a[1], b[0],b[1], c[0],c[1], d_out[0],d_out[1]));
        // reverse orientation flips sign
        printf("inside cw:  %d (should be -1)\n",
            incircle(a[0],a[1], c[0],c[1], b[0],b[1], d_in[0],d_in[1]));
    }

    // ==========================================================================
    // INSPHERE
    // ==========================================================================
    {
        printf("\n=== insphere ===\n");
        s_point a = {{{ 1,  0,  0}}};
        s_point b = {{{-1,  0,  0}}};
        s_point c = {{{ 0,  1,  0}}};
        s_point d = {{{ 0,  0,  1}}};
        s_point e_on  = {{{ 0, -1,  0}}};
        s_point e_in  = {{{ 0,  0,  0}}};
        s_point e_out = {{{ 2,  0,  0}}};
        printf("on sphere:  %d (should be  0)\n",
            insphere(a.x,a.y,a.z, b.x,b.y,b.z, c.x,c.y,c.z, d.x,d.y,d.z,
                     e_on.x,e_on.y,e_on.z));
        printf("inside:     %d (should be  1)\n",
            insphere(a.x,a.y,a.z, b.x,b.y,b.z, c.x,c.y,c.z, d.x,d.y,d.z,
                     e_in.x,e_in.y,e_in.z));
        printf("outside:    %d (should be -1)\n",
            insphere(a.x,a.y,a.z, b.x,b.y,b.z, c.x,c.y,c.z, d.x,d.y,d.z,
                     e_out.x,e_out.y,e_out.z));
    }

    // ==========================================================================
    // POWERTEST1D / test_orthosegment
    // ==========================================================================
    {
        printf("\n=== powertest1d / test_orthosegment ===\n");

        // Segment xa=0,xb=2 with zero weights.
        // Orthosegment centered at midpoint 1, radius 1.
        // Boundary points: xc=0 (=xa) and xc=2 (=xb).
        // xc=0.5: inside   xc=3: outside   xc=0: on boundary
        double x[2]  = {0, 2};
        double wx[2] = {0, 0};

        // raw predicate: +1 means outside, -1 means inside (for xa < xb)
        printf("--- raw powertest1d ---\n");
        printf("inside:      %d (should be -1)\n", powertest1d(x[0],wx[0], x[1],wx[1], 0.5, 0));
        printf("outside:     %d (should be +1)\n", powertest1d(x[0],wx[0], x[1],wx[1], 3,   0));
        printf("on boundary: %d (should be  0)\n", powertest1d(x[0],wx[0], x[1],wx[1], 0,   0));
        printf("wq=+10:      %d (should be -1)\n", powertest1d(x[0],wx[0], x[1],wx[1], 3,  10));
        printf("wq=-10:      %d (should be +1)\n", powertest1d(x[0],wx[0], x[1],wx[1], 0.5,-10));

        // wrapper: +1 means inside, -1 means outside
        printf("--- test_orthosegment ---\n");
        printf("inside:      %d (should be  1)\n", test_orthosegment(x, wx, 0.5, 0));
        printf("outside:     %d (should be -1)\n", test_orthosegment(x, wx, 3,   0));
        printf("on boundary: %d (should be  0)\n", test_orthosegment(x, wx, 0,   0));
        printf("wq=+10:      %d (should be  1)\n", test_orthosegment(x, wx, 3,  10));
        printf("wq=-10:      %d (should be -1)\n", test_orthosegment(x, wx, 0.5,-10));

        // reversed segment — same geometric result due to orientation handling
        double x_rev[2]  = {2, 0};
        double wx_rev[2] = {0, 0};
        printf("--- reversed segment ---\n");
        printf("inside:      %d (should be  1)\n", test_orthosegment(x_rev, wx_rev, 0.5, 0));
        printf("outside:     %d (should be -1)\n", test_orthosegment(x_rev, wx_rev, 3,   0));
        printf("on boundary: %d (should be  0)\n", test_orthosegment(x_rev, wx_rev, 0,   0));

        // Verify test_orthosegment agrees with test_orthocircle
        // for a 1D slice: segment from (-1,0) to (1,0) on x-axis
        double seg[2]    = {-1, 1};
        double wseg[2]   = {0, 0};
        printf("1d slice inside:  %d (should be  1)\n", test_orthosegment(seg, wseg,  0,   0));
        printf("1d slice outside: %d (should be -1)\n", test_orthosegment(seg, wseg,  2,   0));
        printf("1d slice on:      %d (should be  0)\n", test_orthosegment(seg, wseg, -1,   0));
    }

    // ==========================================================================
    // POWERTEST2D / test_orthocircle
    // ==========================================================================
    {
        printf("\n=== powertest2d / test_orthocircle ===\n");
        // unit circle, zero weights
        double a[2] = { 1,  0};
        double b[2] = { 0,  1};
        double c[2] = {-1,  0};
        double q_on[2]  = { 0, -1};
        double q_in[2]  = { 0,  0};
        double q_out[2] = { 2,  0};
        printf("on circle:  %d (should be  0)\n",
            test_orthocircle(a,0, b,0, c,0, q_on,0));
        printf("inside:     %d (should be  1)\n",
            test_orthocircle(a,0, b,0, c,0, q_in,0));
        printf("outside:    %d (should be -1)\n",
            test_orthocircle(a,0, b,0, c,0, q_out,0));
        // weight pulls q_on inside
        printf("on+wq=1:    %d (should be  1)\n",
            test_orthocircle(a,0, b,0, c,0, q_on,1));
        printf("on+wq=-1:   %d (should be -1)\n",
            test_orthocircle(a,0, b,0, c,0, q_on,-1));
        // agree with incircle when weights zero
        printf("agree incircle inside:  %d vs %d (should match)\n",
            test_orthocircle(a,0, b,0, c,0, q_in,0),
            test_incircle(a, b, c, q_in));
        printf("agree incircle outside: %d vs %d (should match)\n",
            test_orthocircle(a,0, b,0, c,0, q_out,0),
            test_incircle(a, b, c, q_out));
    }

    // ==========================================================================
    // POWERTEST3D / test_orthosphere
    // ==========================================================================
    {
        printf("\n=== powertest3d / test_orthosphere ===\n");
        s_point a = {{{ 1,  0,  0}}};
        s_point b = {{{-1,  0,  0}}};
        s_point c = {{{ 0,  1,  0}}};
        s_point d = {{{ 0,  0,  1}}};
        double w[] = {0, 0, 0, 0};
        s_point q_on  = {{{ 0, -1,  0}}};
        s_point q_in  = {{{ 0,  0,  0}}};
        s_point q_out = {{{ 2,  0,  0}}};
        printf("on sphere:  %d (should be  0)\n",
            test_orthosphere((s_point[]){a,b,c,d}, w, q_on,  0));
        printf("inside:     %d (should be  1)\n",
            test_orthosphere((s_point[]){a,b,c,d}, w, q_in,  0));
        printf("outside:    %d (should be -1)\n",
            test_orthosphere((s_point[]){a,b,c,d}, w, q_out, 0));
        // weight pulls q_on inside/outside
        printf("on+wq=+1:   %d (should be  1)\n",
            test_orthosphere((s_point[]){a,b,c,d}, w, q_on,  1));
        printf("on+wq=-1:   %d (should be -1)\n",
            test_orthosphere((s_point[]){a,b,c,d}, w, q_on, -1));
        // agree with insphere when weights zero
        printf("agree insphere inside:  %d vs %d (should match)\n",
            test_orthosphere((s_point[]){a,b,c,d}, w, q_in,  0),
            test_insphere((s_point[]){a,b,c,d}, q_in));
        printf("agree insphere outside: %d vs %d (should match)\n",
            test_orthosphere((s_point[]){a,b,c,d}, w, q_out, 0),
            test_insphere((s_point[]){a,b,c,d}, q_out));
        printf("agree insphere on:      %d vs %d (should match)\n",
            test_orthosphere((s_point[]){a,b,c,d}, w, q_on,  0),
            test_insphere((s_point[]){a,b,c,d}, q_on));
    }


    // // RAW INSPHERE
    // {
    // s_point pa = {{{1,  0,  0}}};
    // s_point pb = {{{-1, 0,  0}}};
    // s_point pc = {{{0,  1,  0}}};
    // s_point pd = {{{0,  0,  1}}};  // ← z != 0, breaks coplanarity
    // s_point pe_on  = {{{0, -1,  0}}};  // on sphere, should give 0
    // s_point pe_in  = {{{0,  0,  0}}};  // inside, should give nonzero
    // s_point pe_out = {{{2,  0,  0}}};  // outside, should give nonzero
    //
    // printf("orientation: %d\n",
    //     orient3d(pa.x,pa.y,pa.z, pb.x,pb.y,pb.z,
    //              pc.x,pc.y,pc.z, pd.x,pd.y,pd.z));
    // printf("on sphere:   %d (should be 0)\n",
    //     insphere(pa.x,pa.y,pa.z, pb.x,pb.y,pb.z,
    //              pc.x,pc.y,pc.z, pd.x,pd.y,pd.z,
    //              pe_on.x,pe_on.y,pe_on.z));
    // printf("inside:      %d\n",
    //     insphere(pa.x,pa.y,pa.z, pb.x,pb.y,pb.z,
    //              pc.x,pc.y,pc.z, pd.x,pd.y,pd.z,
    //              pe_in.x,pe_in.y,pe_in.z));
    // printf("outside:     %d\n",
    //     insphere(pa.x,pa.y,pa.z, pb.x,pb.y,pb.z,
    //              pc.x,pc.y,pc.z, pd.x,pd.y,pd.z,
    //              pe_out.x,pe_out.y,pe_out.z));
    // printf("powertest3d inside: %d\n",
    //     powertest3d(pa.x,pa.y,pa.z,0, pb.x,pb.y,pb.z,0,
    //                 pc.x,pc.y,pc.z,0, pd.x,pd.y,pd.z,0,
    //                 pe_in.x,pe_in.y,pe_in.z,0));
    // }
    //
    //
    // // INSPHERE_WEIGHTED
    // {
    // s_point pa = {{{1,  0,  0}}};
    // s_point pb = {{{-1, 0,  0}}};
    // s_point pc = {{{0,  1,  0}}};
    // s_point pd = {{{0,  0,  1}}};
    // double w[] = {0, 0, 0, 0};
    // s_point pe_on  = {{{0, -1,  0}}};  // on unit sphere
    // s_point pe_in  = {{{0,  0,  0}}};  // inside
    // s_point pe_out = {{{2,  0,  0}}};  // outside
    //
    // printf("Test exact integer points, wi=0\n");
    // printf("on sphere: %d (should be 0)\n",
    //     test_orthosphere((s_point[]){pa,pb,pc,pd}, w, pe_on, 0));
    // printf("inside:    %d (should be 1)\n",
    //     test_orthosphere((s_point[]){pa,pb,pc,pd}, w, pe_in, 0));
    // printf("outside:   %d (should be -1)\n",
    //     test_orthosphere((s_point[]){pa,pb,pc,pd}, w, pe_out, 0));
    //
    // // now test weights: pe_on has distance 0, adding wq moves it in or out
    // printf("on sphere wq=+1: %d (should be 1, weight pulls in)\n",
    //     test_orthosphere((s_point[]){pa,pb,pc,pd}, w, pe_on, 1));
    // printf("on sphere wq=-1: %d (should be -1, weight pushes out)\n",
    //     test_orthosphere((s_point[]){pa,pb,pc,pd}, w, pe_on, -1));
    // }
    //


    puts("\n\n\n");
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
