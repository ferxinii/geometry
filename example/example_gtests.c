#include "points.h"
#include "gtests.h"
#include <stdio.h>
#include <math.h>

// ==========================================================================
// ORIENT2D
// sign convention: +1 if p is to the LEFT of the directed line a->b
//                  (i.e. a,b,p form a CCW triangle)
// ==========================================================================
static void test_orient2d(void)
{
    printf("\n=== test_orientation_2d ===\n");

    s_point2d a = {{0, 0}};
    s_point2d b = {{1, 0}};

    // c = (0,1): left of a->b (CCW) => +1
    s_point2d left     = {{ 0,  1}};
    // c = (0,-1): right of a->b (CW) => -1
    s_point2d right    = {{ 0, -1}};
    // c = (2,0): collinear => 0
    s_point2d col      = {{ 2,  0}};

    printf("left  of a->b:   %d (want  1)\n", test_orientation_2d((s_point2d[]){a,b}, left));
    printf("right of a->b:   %d (want -1)\n", test_orientation_2d((s_point2d[]){a,b}, right));
    printf("collinear:       %d (want  0)\n", test_orientation_2d((s_point2d[]){a,b}, col));

    // Swap endpoints: flips sign
    printf("left  of b->a:   %d (want -1)\n", test_orientation_2d((s_point2d[]){b,a}, left));
    printf("right of b->a:   %d (want  1)\n", test_orientation_2d((s_point2d[]){b,a}, right));

    // Different line: p->q along +y axis
    s_point2d p = {{0, 0}};
    s_point2d q = {{0, 1}};
    s_point2d r_left  = {{-1, 0}};  // left  of p->q
    s_point2d r_right = {{ 1, 0}};  // right of p->q
    printf("left  of p->q:   %d (want  1)\n", test_orientation_2d((s_point2d[]){p,q}, r_left));
    printf("right of p->q:   %d (want -1)\n", test_orientation_2d((s_point2d[]){p,q}, r_right));
}

// ==========================================================================
// ORIENT3D
// From output: test_orientation returns -1 for a point ABOVE the plane
// defined by CCW-ordered (a,b,c), +1 for BELOW, 0 for coplanar.
// ==========================================================================
static void test_orient3d(void)
{
    printf("\n=== test_orientation ===\n");

    s_point a = {{{0, 0, 0}}};
    s_point b = {{{1, 0, 0}}};
    s_point c = {{{0, 1, 0}}};
    // Normal of (a,b,c) by right-hand rule: +z
    s_point above = {{{0, 0,  1}}};   // above xy-plane => -1
    s_point below = {{{0, 0, -1}}};   // below xy-plane => +1
    s_point cop   = {{{1, 1,  0}}};   // coplanar       =>  0

    printf("above plane: %d (want -1)\n", test_orientation((s_point[]){a,b,c}, above));
    printf("below plane: %d (want  1)\n", test_orientation((s_point[]){a,b,c}, below));
    printf("coplanar:    %d (want  0)\n", test_orientation((s_point[]){a,b,c}, cop));

    // Swap b,c flips normal => flips sign
    printf("swap b,c above: %d (want  1)\n", test_orientation((s_point[]){a,c,b}, above));
    printf("swap b,c below: %d (want -1)\n", test_orientation((s_point[]){a,c,b}, below));

    // Tilted plane: x+y+z=1, CCW normal=(1,1,1)/sqrt(3)
    s_point p0 = {{{1, 0, 0}}};
    s_point p1 = {{{0, 1, 0}}};
    s_point p2 = {{{0, 0, 1}}};
    s_point above2 = {{{1, 1, 1}}};   // positive-normal side => -1
    s_point below2 = {{{0, 0, 0}}};   // negative-normal side => +1
    printf("tilted above: %d (want -1)\n", test_orientation((s_point[]){p0,p1,p2}, above2));
    printf("tilted below: %d (want  1)\n", test_orientation((s_point[]){p0,p1,p2}, below2));
}

// ==========================================================================
// INCIRCLE
// Order-independent: +1=inside, -1=outside, 0=on circle.
// ==========================================================================
static void test_incircle_fn(void)
{
    printf("\n=== test_incircle ===\n");

    // Unit circle: a=(1,0), b=(0,1), c=(-1,0)
    s_point2d a = {{ 1,  0}};
    s_point2d b = {{ 0,  1}};
    s_point2d c = {{-1,  0}};
    s_point2d on  = {{ 0, -1}};
    s_point2d in  = {{ 0,  0}};
    s_point2d out = {{ 2,  0}};

    printf("on circle:      %d (want  0)\n", test_incircle((s_point2d[]){a,b,c}, on));
    printf("inside:         %d (want  1)\n", test_incircle((s_point2d[]){a,b,c}, in));
    printf("outside:        %d (want -1)\n", test_incircle((s_point2d[]){a,b,c}, out));

    // Order-independent
    printf("CW order in:    %d (want  1)\n", test_incircle((s_point2d[]){a,c,b}, in));
    printf("CW order out:   %d (want -1)\n", test_incircle((s_point2d[]){a,c,b}, out));
    printf("rot. order in:  %d (want  1)\n", test_incircle((s_point2d[]){b,c,a}, in));
    printf("rot. order out: %d (want -1)\n", test_incircle((s_point2d[]){b,c,a}, out));

    // Equilateral triangle with circumradius=1: (1,0), (-1/2,sqrt(3)/2), (-1/2,-sqrt(3)/2)
    double sq3 = sqrt(3.0);
    s_point2d e0 = {{ 1.0,      0.0     }};
    s_point2d e1 = {{-0.5,  sq3/2.0}};
    s_point2d e2 = {{-0.5, -sq3/2.0}};
    s_point2d ecen = {{0,   0  }};   // center, inside
    s_point2d eout = {{1.5, 0  }};   // outside
    printf("equil. inside:  %d (want  1)\n", test_incircle((s_point2d[]){e0,e1,e2}, ecen));
    printf("equil. on:      %d (want  0)\n", test_incircle((s_point2d[]){e0,e1,e2}, e0));
    printf("equil. outside: %d (want -1)\n", test_incircle((s_point2d[]){e0,e1,e2}, eout));
}

// ==========================================================================
// INSPHERE
// Order-independent. +1=inside, -1=outside, 0=on sphere.
// ==========================================================================
static void test_insphere_fn(void)
{
    printf("\n=== test_insphere ===\n");

    s_point a = {{{ 1,  0,  0}}};
    s_point b = {{{-1,  0,  0}}};
    s_point c = {{{ 0,  1,  0}}};
    s_point d = {{{ 0,  0,  1}}};
    s_point on  = {{{ 0, -1,  0}}};
    s_point in  = {{{ 0,  0,  0}}};
    s_point out = {{{ 2,  0,  0}}};

    printf("on sphere:      %d (want  0)\n", test_insphere((s_point[]){a,b,c,d}, on));
    printf("inside:         %d (want  1)\n", test_insphere((s_point[]){a,b,c,d}, in));
    printf("outside:        %d (want -1)\n", test_insphere((s_point[]){a,b,c,d}, out));

    // Order-independent
    printf("swap a,b in:    %d (want  1)\n", test_insphere((s_point[]){b,a,c,d}, in));
    printf("swap a,b out:   %d (want -1)\n", test_insphere((s_point[]){b,a,c,d}, out));
    printf("rot. order in:  %d (want  1)\n", test_insphere((s_point[]){b,c,d,a}, in));
    printf("rot. order out: %d (want -1)\n", test_insphere((s_point[]){b,c,d,a}, out));

    s_point barely_in  = {{{ 0.999,  0,  0}}};
    s_point barely_out = {{{ 1.001,  0,  0}}};
    printf("barely inside:  %d (want  1)\n", test_insphere((s_point[]){a,b,c,d}, barely_in));
    printf("barely outside: %d (want -1)\n", test_insphere((s_point[]){a,b,c,d}, barely_out));
}

// ==========================================================================
// TEST_ORTHOSEGMENT (1D)
//
// k=2: wrapper returns +1=inside (pi<0), -1=outside (pi>0), 0=on.
// k=1: returns sign(pi) directly: -1=pi<0=inside, +1=pi>0=outside, 0=on.
//      (no orientation correction for k=1)
// ==========================================================================
static void test_orthosegment_fn(void)
{
    printf("\n=== test_orthosegment ===\n");

    // --- k=1 ---
    // Orthosegment of {(xa=1, wa=4)}: S = (1, -4), so pi(xp,wp) = (xp-1)^2 + 4 - wp.
    // pi < 0 <=> wp > (xp-1)^2 + 4.
    // With wp=0: pi = (xp-1)^2 + 4 >= 4 > 0 always => always outside.
    // Need large wp to get inside.
    {
        double c[1]  = {1};
        double wc[1] = {4};
        printf("-- k=1, c={1}, wc={4}: pi(xp,wp) = (xp-1)^2 + 4 - wp --\n");
        // wp=0: pi = (xp-1)^2 + 4 > 0 always => outside (-1)
        printf("xp=1, wp=0:   %d (want -1)\n", test_orthosegment(1, c, wc,  1,  0, 0));
        printf("xp=0, wp=0:   %d (want -1)\n", test_orthosegment(1, c, wc,  0,  0, 0));
        printf("xp=5, wp=0:   %d (want -1)\n", test_orthosegment(1, c, wc,  5,  0, 0));
        // wp = (xp-1)^2 + 4 => pi=0 => on boundary
        // xp=1: pi = 0 + 4 - wp = 0 => wp=4
        printf("xp=1, wp=4:   %d (want  0)\n", test_orthosegment(1, c, wc,  1,  4, 0));
        // xp=3: pi = 4 + 4 - wp = 0 => wp=8
        printf("xp=3, wp=8:   %d (want  0)\n", test_orthosegment(1, c, wc,  3,  8, 0));
        // xp=1, wp=5: pi = 0 + 4 - 5 = -1 < 0 => inside
        printf("xp=1, wp=5:   %d (want  1)\n", test_orthosegment(1, c, wc,  1,  5, 0));
        // xp=3, wp=9: pi = 4 + 4 - 9 = -1 < 0 => inside
        printf("xp=3, wp=9:   %d (want  1)\n", test_orthosegment(1, c, wc,  3,  9, 0));
        // xp=3, wp=7: pi = 4 + 4 - 7 = 1 > 0 => outside
        printf("xp=3, wp=7:   %d (want -1)\n", test_orthosegment(1, c, wc,  3,  7, 0));

        // alpha variant: test pi < alpha, i.e. (xp-1)^2 + 4 - wp < alpha
        // xp=1, wp=0: pi=4. Test pi < alpha.
        printf("-- k=1, alpha variant (pi=4 at xp=1,wp=0) --\n");
        printf("pi=4 < alpha=5: %d (want  1)\n", test_orthosegment(1, c, wc, 1, 0, 5));
        printf("Now calling problematic...\n");
        printf("pi=4 < alpha=4: %d (want  0)\n", test_orthosegment(1, c, wc, 1, 0, 4));
        printf("pi=4 < alpha=3: %d (want -1)\n", test_orthosegment(1, c, wc, 1, 0, 3));
    }

    // --- k=2, zero weights ---
    // Orthosegment of {0,2}: center=1, r^2=1. pi(xp,wp) = (xp-1)^2 - 1 - wp.
    {
        double c[2]  = {0, 2};
        double wc[2] = {0, 0};
        printf("-- k=2, zero weights: center=1, r=1 --\n");
        printf("inside  (xp=0.5):   %d (want  1)\n", test_orthosegment(2, c, wc, 0.5, 0,   0));
        printf("outside (xp=3):     %d (want -1)\n", test_orthosegment(2, c, wc, 3,   0,   0));
        printf("on bdry (xp=0):     %d (want  0)\n", test_orthosegment(2, c, wc, 0,   0,   0));
        printf("on bdry (xp=2):     %d (want  0)\n", test_orthosegment(2, c, wc, 2,   0,   0));
        // wp=10: pi(xp=3) = 4-1-10 = -7 < 0 => inside
        printf("out+wp=10 -> in:    %d (want  1)\n", test_orthosegment(2, c, wc, 3,   10,  0));
        // wp=-10: pi(xp=0.5) = 0.25-1+10 = 9.25 > 0 => outside
        printf("in+wp=-10 -> out:   %d (want -1)\n", test_orthosegment(2, c, wc, 0.5,-10,  0));
    }

    // --- k=2, reversed endpoints: same result ---
    {
        double c[2]  = {2, 0};
        double wc[2] = {0, 0};
        printf("-- k=2 reversed --\n");
        printf("inside  (xp=0.5):   %d (want  1)\n", test_orthosegment(2, c, wc, 0.5, 0, 0));
        printf("outside (xp=3):     %d (want -1)\n", test_orthosegment(2, c, wc, 3,   0, 0));
        printf("on bdry (xp=0):     %d (want  0)\n", test_orthosegment(2, c, wc, 0,   0, 0));
    }

    // --- k=2, wider segment: c={0,4}, center=2, r=2 ---
    {
        double c[2]  = {0, 4};
        double wc[2] = {0, 0};
        printf("-- k=2, c={0,4}, center=2, r=2 --\n");
        printf("inside  (xp=2): %d (want  1)\n", test_orthosegment(2, c, wc, 2, 0, 0));
        printf("on bdry (xp=0): %d (want  0)\n", test_orthosegment(2, c, wc, 0, 0, 0));
        printf("on bdry (xp=4): %d (want  0)\n", test_orthosegment(2, c, wc, 4, 0, 0));
        printf("outside (xp=5): %d (want -1)\n", test_orthosegment(2, c, wc, 5, 0, 0));
    }

    // --- alpha variant ---
    // pi at xp=3 (c={0,2}, w=0) = (3-1)^2 - 1 = 3. Test pi < alpha.
    {
        double c[2]  = {0, 2};
        double wc[2] = {0, 0};
        printf("-- k=2, alpha variant (pi=3 at xp=3) --\n");
        printf("pi=3 < alpha=4: %d (want  1)\n", test_orthosegment(2, c, wc, 3, 0, 4));
        printf("pi=3 < alpha=3: %d (want  0)\n", test_orthosegment(2, c, wc, 3, 0, 3));
        printf("pi=3 < alpha=2: %d (want -1)\n", test_orthosegment(2, c, wc, 3, 0, 2));
    }
}

// ==========================================================================
// TEST_ORTHOCIRCLE (2D)
//
// k=1: returns sign(pi) directly: -1=inside, +1=outside, 0=on.
// k=2: returns +1=inside, -1=outside, 0=on.
// k=3: returns +1=inside, -1=outside, 0=on.
// ==========================================================================
static void test_orthocircle_fn(void)
{
    printf("\n=== test_orthocircle ===\n");

    // --- k=3, unit circle, zero weights ---
    // Center=(0,0), r^2=1. pi(p,wp) = |p|^2 - 1 - wp.
    {
        s_point2d a = {{ 1,  0}};
        s_point2d b = {{ 0,  1}};
        s_point2d c = {{-1,  0}};
        double wc[3] = {0, 0, 0};
        s_point2d on  = {{ 0, -1}};
        s_point2d in  = {{ 0,  0}};
        s_point2d out = {{ 2,  0}};

        printf("-- k=3, unit circle, zero weights --\n");
        printf("on:        %d (want  0)\n", test_orthocircle(3,(s_point2d[]){a,b,c},wc, on,  0, 0));
        printf("inside:    %d (want  1)\n", test_orthocircle(3,(s_point2d[]){a,b,c},wc, in,  0, 0));
        printf("outside:   %d (want -1)\n", test_orthocircle(3,(s_point2d[]){a,b,c},wc, out, 0, 0));
        // wp shifts pi: pi(on)=0; +wp => pi=-wp<0 => inside; -wp => pi=+>0 => outside
        printf("on+wp=+1:  %d (want  1)\n", test_orthocircle(3,(s_point2d[]){a,b,c},wc, on,  1, 0));
        printf("on+wp=-1:  %d (want -1)\n", test_orthocircle(3,(s_point2d[]){a,b,c},wc, on, -1, 0));

        // Agree with test_incircle at zero weights
        int ic_on  = test_incircle((s_point2d[]){a,b,c}, on);
        int ic_in  = test_incircle((s_point2d[]){a,b,c}, in);
        int ic_out = test_incircle((s_point2d[]){a,b,c}, out);
        printf("agree incircle on:      ortho=%d incircle=%d %s\n",
            test_orthocircle(3,(s_point2d[]){a,b,c},wc, on,  0, 0), ic_on,
            (test_orthocircle(3,(s_point2d[]){a,b,c},wc, on,  0, 0)==ic_on) ? "(match)" : "*** MISMATCH ***");
        printf("agree incircle inside:  ortho=%d incircle=%d %s\n",
            test_orthocircle(3,(s_point2d[]){a,b,c},wc, in,  0, 0), ic_in,
            (test_orthocircle(3,(s_point2d[]){a,b,c},wc, in,  0, 0)==ic_in) ? "(match)" : "*** MISMATCH ***");
        printf("agree incircle outside: ortho=%d incircle=%d %s\n",
            test_orthocircle(3,(s_point2d[]){a,b,c},wc, out, 0, 0), ic_out,
            (test_orthocircle(3,(s_point2d[]){a,b,c},wc, out, 0, 0)==ic_out) ? "(match)" : "*** MISMATCH ***");

        // Alpha: pi at out=(2,0) = 4-1 = 3
        printf("pi<4:      %d (want  1)\n", test_orthocircle(3,(s_point2d[]){a,b,c},wc, out, 0, 4));
        printf("pi<3:      %d (want  0)\n", test_orthocircle(3,(s_point2d[]){a,b,c},wc, out, 0, 3));
        printf("pi<2:      %d (want -1)\n", test_orthocircle(3,(s_point2d[]){a,b,c},wc, out, 0, 2));
    }

    // --- k=3, nonzero equal weights ---
    // Adding w=1 to all inputs: wS_new = 2. pi_new(out) = 4-2 = 2 > 0 => outside.
    // With wp=3: pi_new - wp = 2-3 = -1 < 0 => inside. With wp=2: pi_new-wp=0 => on.
    {
        s_point2d a = {{ 1,  0}};
        s_point2d b = {{ 0,  1}};
        s_point2d c = {{-1,  0}};
        double wc[3] = {1, 1, 1};
        s_point2d out = {{ 2,  0}};
        printf("-- k=3, w=1 on all inputs --\n");
        printf("out:           %d (want -1)\n", test_orthocircle(3,(s_point2d[]){a,b,c},wc, out, 0, 0));
        printf("out+wp=4:      %d (want  0)\n", test_orthocircle(3,(s_point2d[]){a,b,c},wc, out, 4, 0));
        printf("out+wp=5:      %d (want  1)\n", test_orthocircle(3,(s_point2d[]){a,b,c},wc, out, 5, 0));
    }

    // --- k=2 ---
    // Orthocircle of a=(0,0,w=0), b=(2,0,w=0): center=(1,0), r^2=1.
    // pi(p,wp) = (px-1)^2 + py^2 - 1 - wp.
    // On-sphere: (px-1)^2 + py^2 = 1. Points: (0,0),(2,0),(1,1),(1,-1).
    {
        s_point2d a = {{0, 0}};
        s_point2d b = {{2, 0}};
        double wc[2] = {0, 0};
        s_point2d in     = {{1,  0}};   // center, pi=-1
        s_point2d on_a   = {{0,  0}};  // pi=0
        s_point2d on_b   = {{2,  0}};  // pi=0
        s_point2d on_top = {{1,  1}};  // pi=0
        s_point2d on_bot = {{1, -1}};  // pi=0
        s_point2d out    = {{3,  0}};  // pi=3

        printf("-- k=2, c={(0,0),(2,0)}, center=(1,0), r=1 --\n");
        printf("inside (center):  %d (want  1)\n", test_orthocircle(2,(s_point2d[]){a,b},wc, in,     0, 0));
        printf("on (=a):          %d (want  0)\n", test_orthocircle(2,(s_point2d[]){a,b},wc, on_a,   0, 0));
        printf("on (=b):          %d (want  0)\n", test_orthocircle(2,(s_point2d[]){a,b},wc, on_b,   0, 0));
        printf("on (top):         %d (want  0)\n", test_orthocircle(2,(s_point2d[]){a,b},wc, on_top, 0, 0));
        printf("on (bot):         %d (want  0)\n", test_orthocircle(2,(s_point2d[]){a,b},wc, on_bot, 0, 0));
        printf("outside:          %d (want -1)\n", test_orthocircle(2,(s_point2d[]){a,b},wc, out,    0, 0));
        printf("on+wp=+1:         %d (want  1)\n", test_orthocircle(2,(s_point2d[]){a,b},wc, on_a,   1, 0));
        printf("on+wp=-1:         %d (want -1)\n", test_orthocircle(2,(s_point2d[]){a,b},wc, on_a,  -1, 0));
        // Reversed input: same result
        printf("inside (rev):     %d (want  1)\n", test_orthocircle(2,(s_point2d[]){b,a},wc, in,     0, 0));
        printf("outside (rev):    %d (want -1)\n", test_orthocircle(2,(s_point2d[]){b,a},wc, out,    0, 0));
        // Alpha: pi at out = 3
        printf("pi<4:             %d (want  1)\n", test_orthocircle(2,(s_point2d[]){a,b},wc, out, 0, 4));
        printf("pi<3:             %d (want  0)\n", test_orthocircle(2,(s_point2d[]){a,b},wc, out, 0, 3));
        printf("pi<2:             %d (want -1)\n", test_orthocircle(2,(s_point2d[]){a,b},wc, out, 0, 2));
    }

    // --- k=1 ---
    // S centered at a=(0,0), r^2=wa=4. pi(p,S) = |p|^2 + 4 - wp.
    {
        s_point2d a  = {{0, 0}};
        double wc[1] = {4};
        s_point2d in  = {{1, 0}};  double win = 6;   // pi = 1+4-6 = -1
        s_point2d on  = {{2, 0}};  double won = 8;   // pi = 4+4-8 = 0
        s_point2d out = {{3, 0}};  double wout = 0;  // pi = 9+4 = 13
        printf("-- k=1, center=(0,0), w=4 (sign(pi) directly: -1=in, +1=out) --\n");
        printf("inside  (pi=-1): %d (want  1)\n", test_orthocircle(1,(s_point2d[]){a},wc, in,  win, 0));
        printf("on bdry (pi=0):  %d (want  0)\n", test_orthocircle(1,(s_point2d[]){a},wc, on,  won, 0));
        printf("outside (pi=13):  %d (want -1)\n", test_orthocircle(1,(s_point2d[]){a},wc, out, wout, 0));
    }
}

// ==========================================================================
// TEST_ORTHOSPHERE (3D)
//
// k=1: returns sign(pi) directly: -1=inside, +1=outside, 0=on.
// k=2: returns +1=inside, -1=outside, 0=on.
// k=3: returns +1=inside, -1=outside, 0=on.
// k=4: returns +1=inside, -1=outside, 0=on.
// ==========================================================================
static void test_orthosphere_fn(void)
{
    printf("\n=== test_orthosphere ===\n");

    // --- k=4, unit sphere, zero weights ---
    // Center=(0,0,0), r^2=1. pi(p,wp) = |p|^2 - 1 - wp.
    {
        s_point a = {{{ 1,  0,  0}}};
        s_point b = {{{-1,  0,  0}}};
        s_point c = {{{ 0,  1,  0}}};
        s_point d = {{{ 0,  0,  1}}};
        double wc[4] = {0, 0, 0, 0};
        s_point on  = {{{ 0, -1,  0}}};
        s_point in  = {{{ 0,  0,  0}}};
        s_point out = {{{ 2,  0,  0}}};

        printf("-- k=4, unit sphere, zero weights --\n");
        printf("on:             %d (want  0)\n", test_orthosphere(4,(s_point[]){a,b,c,d},wc, on,  0, 0));
        printf("inside:         %d (want  1)\n", test_orthosphere(4,(s_point[]){a,b,c,d},wc, in,  0, 0));
        printf("outside:        %d (want -1)\n", test_orthosphere(4,(s_point[]){a,b,c,d},wc, out, 0, 0));
        printf("on+wp=+1:       %d (want  1)\n", test_orthosphere(4,(s_point[]){a,b,c,d},wc, on,  1, 0));
        printf("on+wp=-1:       %d (want -1)\n", test_orthosphere(4,(s_point[]){a,b,c,d},wc, on, -1, 0));

        // Agree with test_insphere
        int is_on  = test_insphere((s_point[]){a,b,c,d}, on);
        int is_in  = test_insphere((s_point[]){a,b,c,d}, in);
        int is_out = test_insphere((s_point[]){a,b,c,d}, out);
        printf("agree insphere on:      ortho=%d insphere=%d %s\n",
            test_orthosphere(4,(s_point[]){a,b,c,d},wc, on,  0, 0), is_on,
            (test_orthosphere(4,(s_point[]){a,b,c,d},wc, on,  0, 0)==is_on) ? "(match)" : "*** MISMATCH ***");
        printf("agree insphere inside:  ortho=%d insphere=%d %s\n",
            test_orthosphere(4,(s_point[]){a,b,c,d},wc, in,  0, 0), is_in,
            (test_orthosphere(4,(s_point[]){a,b,c,d},wc, in,  0, 0)==is_in) ? "(match)" : "*** MISMATCH ***");
        printf("agree insphere outside: ortho=%d insphere=%d %s\n",
            test_orthosphere(4,(s_point[]){a,b,c,d},wc, out, 0, 0), is_out,
            (test_orthosphere(4,(s_point[]){a,b,c,d},wc, out, 0, 0)==is_out) ? "(match)" : "*** MISMATCH ***");

        // Alpha: pi at out=(2,0,0) = 4-1 = 3
        printf("pi<4:           %d (want  1)\n", test_orthosphere(4,(s_point[]){a,b,c,d},wc, out, 0, 4));
        printf("pi<3:           %d (want  0)\n", test_orthosphere(4,(s_point[]){a,b,c,d},wc, out, 0, 3));
        printf("pi<2:           %d (want -1)\n", test_orthosphere(4,(s_point[]){a,b,c,d},wc, out, 0, 2));

        s_point barely_in  = {{{ 0.999,  0,  0}}};
        s_point barely_out = {{{ 1.001,  0,  0}}};
        printf("barely inside:  %d (want  1)\n", test_orthosphere(4,(s_point[]){a,b,c,d},wc, barely_in,  0, 0));
        printf("barely outside: %d (want -1)\n", test_orthosphere(4,(s_point[]){a,b,c,d},wc, barely_out, 0, 0));
    }

    // --- k=3 ---
    // Triangle a=(1,0,0), b=(0,1,0), c=(-1,0,0) in the xy-plane.
    // Circumcenter in xy: cx=0, cy=0 (equidistant from all three at dist=1).
    // Orthosphere center = (0,0,0), wS=r^2=1. pi(p,wp) = |p|^2 - 1 - wp.
    {
        s_point a = {{{ 1,  0,  0}}};
        s_point b = {{{ 0,  1,  0}}};
        s_point c = {{{-1,  0,  0}}};
        double wc[3] = {0, 0, 0};
        s_point on_a    = {{{ 1,  0,  0}}};   // =a, |p|^2=1 => pi=0
        s_point on_z    = {{{ 0,  0,  1}}};   // |p|^2=1 => pi=0
        s_point on_neg  = {{{ 0,  0, -1}}};   // |p|^2=1 => pi=0
        s_point inside  = {{{ 0.5, 0, 0}}};   // |p|^2=0.25 < 1 => inside
        s_point origin  = {{{ 0,  0,  0}}};   // |p|^2=0 => pi=-1 < 0 => inside
        s_point outside = {{{ 0,  0,  2}}};   // |p|^2=4 > 1 => outside

        printf("-- k=3, circumsphere center=origin, r=1 --\n");
        printf("on (=a):           %d (want  0)\n", test_orthosphere(3,(s_point[]){a,b,c},wc, on_a,   0, 0));
        printf("on (0,0,1):        %d (want  0)\n", test_orthosphere(3,(s_point[]){a,b,c},wc, on_z,   0, 0));
        printf("on (0,0,-1):       %d (want  0)\n", test_orthosphere(3,(s_point[]){a,b,c},wc, on_neg, 0, 0));
        printf("inside (0.5,0,0):  %d (want  1)\n", test_orthosphere(3,(s_point[]){a,b,c},wc, inside, 0, 0));
        printf("inside (origin):   %d (want  1)\n", test_orthosphere(3,(s_point[]){a,b,c},wc, origin, 0, 0));
        printf("outside (0,0,2):   %d (want -1)\n", test_orthosphere(3,(s_point[]){a,b,c},wc, outside,0, 0));
        printf("on+wp=+1:          %d (want  1)\n", test_orthosphere(3,(s_point[]){a,b,c},wc, on_a,   1, 0));
        printf("on+wp=-1:          %d (want -1)\n", test_orthosphere(3,(s_point[]){a,b,c},wc, on_a,  -1, 0));
    }

    // --- k=2 along X: |bx-ax| largest => ej1=ey, ej2=ez branch ---
    // a=(0,0,0), b=(2,0,0): center=(1,0,0), r^2=1.
    // pi(p,wp) = (px-1)^2 + py^2 + pz^2 - 1 - wp.
    // On-sphere: (px-1)^2 + py^2 + pz^2 = 1.
    {
        s_point a = {{{0, 0, 0}}};
        s_point b = {{{2, 0, 0}}};
        double wc[2] = {0, 0};
        s_point in    = {{{1, 0, 0}}};   // center, pi=-1
        s_point on_a  = {{{0, 0, 0}}};  // pi=0
        s_point on_b  = {{{2, 0, 0}}};  // pi=0
        s_point on_y  = {{{1, 1, 0}}};  // 0+1+0=1 => pi=0
        s_point on_z  = {{{1, 0, 1}}};  // 0+0+1=1 => pi=0
        s_point out   = {{{3, 0, 0}}};  // pi=4-1=3

        printf("-- k=2, along X (ej1=ey,ej2=ez branch) --\n");
        printf("inside (center): %d (want  1)\n", test_orthosphere(2,(s_point[]){a,b},wc, in,   0, 0));
        printf("on (=a):         %d (want  0)\n", test_orthosphere(2,(s_point[]){a,b},wc, on_a, 0, 0));
        printf("on (=b):         %d (want  0)\n", test_orthosphere(2,(s_point[]){a,b},wc, on_b, 0, 0));
        printf("on (y-side):     %d (want  0)\n", test_orthosphere(2,(s_point[]){a,b},wc, on_y, 0, 0));
        printf("on (z-side):     %d (want  0)\n", test_orthosphere(2,(s_point[]){a,b},wc, on_z, 0, 0));
        printf("outside:         %d (want -1)\n", test_orthosphere(2,(s_point[]){a,b},wc, out,  0, 0));
        printf("on+wp=+1:        %d (want  1)\n", test_orthosphere(2,(s_point[]){a,b},wc, on_a, 1, 0));
        printf("on+wp=-1:        %d (want -1)\n", test_orthosphere(2,(s_point[]){a,b},wc, on_a,-1, 0));
        // Alpha: pi at out = 3
        printf("pi<4:            %d (want  1)\n", test_orthosphere(2,(s_point[]){a,b},wc, out, 0, 4));
        printf("pi<3:            %d (want  0)\n", test_orthosphere(2,(s_point[]){a,b},wc, out, 0, 3));
        printf("pi<2:            %d (want -1)\n", test_orthosphere(2,(s_point[]){a,b},wc, out, 0, 2));
    }

    // --- k=2 along Y: |by-ay| largest => ej1=ez, ej2=ex branch ---
    // a=(0,0,0), b=(0,2,0): center=(0,1,0), r^2=1.
    {
        s_point a = {{{0, 0, 0}}};
        s_point b = {{{0, 2, 0}}};
        double wc[2] = {0, 0};
        s_point in   = {{{0, 1, 0}}};
        s_point on_a = {{{0, 0, 0}}};
        s_point on_x = {{{1, 1, 0}}};   // 1+0+0=1 => on
        s_point on_z = {{{0, 1, 1}}};   // 0+0+1=1 => on
        s_point out  = {{{0, 3, 0}}};

        printf("-- k=2, along Y (ej1=ez,ej2=ex branch) --\n");
        printf("inside:  %d (want  1)\n", test_orthosphere(2,(s_point[]){a,b},wc, in,   0, 0));
        printf("on (=a): %d (want  0)\n", test_orthosphere(2,(s_point[]){a,b},wc, on_a, 0, 0));
        printf("on (x):  %d (want  0)\n", test_orthosphere(2,(s_point[]){a,b},wc, on_x, 0, 0));
        printf("on (z):  %d (want  0)\n", test_orthosphere(2,(s_point[]){a,b},wc, on_z, 0, 0));
        printf("outside: %d (want -1)\n", test_orthosphere(2,(s_point[]){a,b},wc, out,  0, 0));
    }

    // --- k=2 along Z: |bz-az| largest => ej1=ex, ej2=ey branch ---
    // a=(0,0,0), b=(0,0,2): center=(0,0,1), r^2=1.
    {
        s_point a = {{{0, 0, 0}}};
        s_point b = {{{0, 0, 2}}};
        double wc[2] = {0, 0};
        s_point in   = {{{0, 0, 1}}};
        s_point on_a = {{{0, 0, 0}}};
        s_point on_x = {{{1, 0, 1}}};   // 1+0+0=1 => on
        s_point on_y = {{{0, 1, 1}}};   // 0+1+0=1 => on
        s_point out  = {{{0, 0, 3}}};

        printf("-- k=2, along Z (ej1=ex,ej2=ey branch) --\n");
        printf("inside:  %d (want  1)\n", test_orthosphere(2,(s_point[]){a,b},wc, in,   0, 0));
        printf("on (=a): %d (want  0)\n", test_orthosphere(2,(s_point[]){a,b},wc, on_a, 0, 0));
        printf("on (x):  %d (want  0)\n", test_orthosphere(2,(s_point[]){a,b},wc, on_x, 0, 0));
        printf("on (y):  %d (want  0)\n", test_orthosphere(2,(s_point[]){a,b},wc, on_y, 0, 0));
        printf("outside: %d (want -1)\n", test_orthosphere(2,(s_point[]){a,b},wc, out,  0, 0));
    }

    // --- k=2, diagonal b-a=(1,1,0): |bz-az|=0 smallest => ej1=ex, ej2=ey branch ---
    // a=(0,0,0), b=(1,1,0): center=(0.5,0.5,0), r^2=0.5.
    // pi(p,wp) = (px-0.5)^2 + (py-0.5)^2 + pz^2 - 0.5 - wp.
    {
        s_point a    = {{{0,   0,   0  }}};
        s_point b    = {{{1,   1,   0  }}};
        double wc[2] = {0, 0};
        s_point in   = {{{0.5, 0.5, 0  }}};               // center, pi=-0.5
        s_point on_a = {{{0,   0,   0  }}};               // =a, pi=0
        s_point on_z = {{{0, 1, 0}}};         // pz^2=0.5 => pi=0
        s_point out  = {{{2,   2,   0  }}};               // pi=1.5^2+1.5^2-0.5=4

        printf("-- k=2, diagonal (1,1,0), z-smallest branch --\n");
        printf("inside (center):  %d (want  1)\n", test_orthosphere(2,(s_point[]){a,b},wc, in,   0, 0));
        printf("on (=a):          %d (want  0)\n", test_orthosphere(2,(s_point[]){a,b},wc, on_a, 0, 0));
        printf("on (z-offset):    %d (want  0)\n", test_orthosphere(2,(s_point[]){a,b},wc, on_z, 0, 0));
        printf("outside:          %d (want -1)\n", test_orthosphere(2,(s_point[]){a,b},wc, out,  0, 0));
    }

    // --- k=1 ---
    // S centered at a=(0,0,0), r^2=wa=4. pi(p,wp) = |p|^2 + 4 - wp.
    // Returns sign(pi) directly: -1=inside, +1=outside, 0=on.
    {
        s_point a  = {{{0, 0, 0}}};
        double wc[1] = {4};
        s_point in  = {{{1, 0, 0}}};  double win = 6;   // pi=1+4-6=-1;
        s_point on  = {{{2, 0, 0}}};  double won = 8;   // pi=4+4-8=0
        s_point out = {{{3, 0, 0}}};  double wout = 0;  // pi=9+4=113

        printf("-- k=1, center=(0,0,0), w=4 (sign(pi) directly: -1=in, +1=out) --\n");
        printf("inside  (pi=-1): %d (want 1)\n", test_orthosphere(1,(s_point[]){a},wc, in,  win, 0));
        printf("on bdry (pi=0):  %d (want  0)\n", test_orthosphere(1,(s_point[]){a},wc, on,  won, 0));
        printf("outside (pi=13):  %d (want  -1)\n", test_orthosphere(1,(s_point[]){a},wc, out, wout, 0));
    }
}

// ==========================================================================
// CROSS-DIMENSION CONSISTENCY
// ==========================================================================
static void test_cross_consistency(void)
{
    printf("\n=== cross-dimension consistency ===\n");

    // k=3 orthocircle vs incircle (zero weights)
    {
        s_point2d a = {{ 1,  0}};
        s_point2d b = {{ 0,  1}};
        s_point2d c = {{-1,  0}};
        double wc[3] = {0, 0, 0};
        s_point2d pts[5] = {{{ 0,-1}}, {{ 0, 0}}, {{ 2, 0}}, {{ 0.5, 0.3}}, {{-0.5, 0.5}}};
        const char *names[5] = {"on", "in", "out", "diag1", "diag2"};
        for (int i = 0; i < 5; i++) {
            int oc = test_orthocircle(3,(s_point2d[]){a,b,c},wc, pts[i], 0, 0);
            int ic = test_incircle((s_point2d[]){a,b,c}, pts[i]);
            printf("orthocircle vs incircle [%s]: %d vs %d %s\n",
                names[i], oc, ic, (oc==ic) ? "(match)" : "*** MISMATCH ***");
        }
    }

    // k=4 orthosphere vs insphere (zero weights)
    {
        s_point a = {{{ 1,  0,  0}}};
        s_point b = {{{-1,  0,  0}}};
        s_point c = {{{ 0,  1,  0}}};
        s_point d = {{{ 0,  0,  1}}};
        double wc[4] = {0, 0, 0, 0};
        s_point pts[5] = {
            {{{ 0,-1, 0}}}, {{{0, 0, 0}}}, {{{2, 0, 0}}},
            {{{0.5, 0.3, 0.2}}}, {{{0, 0.5, 0.5}}}
        };
        const char *names[5] = {"on", "in", "out", "diag1", "diag2"};
        for (int i = 0; i < 5; i++) {
            int os = test_orthosphere(4,(s_point[]){a,b,c,d},wc, pts[i], 0, 0);
            int is = test_insphere((s_point[]){a,b,c,d}, pts[i]);
            printf("orthosphere vs insphere [%s]: %d vs %d %s\n",
                names[i], os, is, (os==is) ? "(match)" : "*** MISMATCH ***");
        }
    }
}

// ==========================================================================
// MAIN
// ==========================================================================
int main(void)
{
    test_orient2d();
    test_orient3d();
    test_incircle_fn();
    test_insphere_fn();
    test_orthosegment_fn();
    test_orthocircle_fn();
    test_orthosphere_fn();
    test_cross_consistency();
    printf("\nDone.\n");
    return 0;
}



// #include "points.h"
// #include "gtests.h"
// #include <stdio.h>
// #include <math.h>
//
// int main() 
// {
//     // ==========================================================================
//     // ORIENT2D
//     // ==========================================================================
//     {
//         printf("\n=== orient2d ===\n");
//         // counterclockwise triangle
//         double a[2] = {0, 0};
//         double b[2] = {1, 0};
//         double c[2] = {0, 1};
//         printf("ccw:        %d (should be  1)\n", orient2d(a[0],a[1], b[0],b[1], c[0],c[1]));
//         printf("cw:         %d (should be -1)\n", orient2d(a[0],a[1], c[0],c[1], b[0],b[1]));
//         printf("collinear:  %d (should be  0)\n", orient2d(a[0],a[1], b[0],b[1], b[0],b[1]));
//     }
//
//     // ==========================================================================
//     // ORIENT3D
//     // ==========================================================================
//     {
//         printf("\n=== orient3d ===\n");
//         // positive orientation tetrahedron
//         s_point a = {{{0, 0, 0}}};
//         s_point b = {{{1, 0, 0}}};
//         s_point c = {{{0, 1, 0}}};
//         s_point d = {{{0, 0, 1}}};  // same side as normal: negative
//         s_point e = {{{0, 0,-1}}};  // opposite side: positive
//         printf("same side as plane normal:   %d (should be  -1)\n",
//             orient3d(a.x,a.y,a.z, b.x,b.y,b.z, c.x,c.y,c.z, d.x,d.y,d.z));
//         printf("opposite side as plane normal:   %d (should be +1)\n",
//             orient3d(a.x,a.y,a.z, b.x,b.y,b.z, c.x,c.y,c.z, e.x,e.y,e.z));
//         printf("coplanar:   %d (should be  0)\n",
//             orient3d(a.x,a.y,a.z, b.x,b.y,b.z, c.x,c.y,c.z, a.x,a.y,a.z));
//     }
//
//     // ==========================================================================
//     // INCIRCLE
//     // ==========================================================================
//     {
//         printf("\n=== incircle ===\n");
//         // unit circle: a=(1,0), b=(0,1), c=(-1,0)
//         double a[2] = { 1,  0};
//         double b[2] = { 0,  1};
//         double c[2] = {-1,  0};
//         double d_on[2]  = { 0, -1};  // on circle
//         double d_in[2]  = { 0,  0};  // inside
//         double d_out[2] = { 2,  0};  // outside
//         printf("on circle:  %d (should be  0)\n",
//             incircle(a[0],a[1], b[0],b[1], c[0],c[1], d_on[0],d_on[1]));
//         printf("inside:     %d (should be  1)\n",
//             incircle(a[0],a[1], b[0],b[1], c[0],c[1], d_in[0],d_in[1]));
//         printf("outside:    %d (should be -1)\n",
//             incircle(a[0],a[1], b[0],b[1], c[0],c[1], d_out[0],d_out[1]));
//         // reverse orientation flips sign
//         printf("inside cw:  %d (should be -1)\n",
//             incircle(a[0],a[1], c[0],c[1], b[0],b[1], d_in[0],d_in[1]));
//     }
//
//     // ==========================================================================
//     // INSPHERE
//     // ==========================================================================
//     {
//         printf("\n=== insphere ===\n");
//         s_point a = {{{ 1,  0,  0}}};
//         s_point b = {{{-1,  0,  0}}};
//         s_point c = {{{ 0,  1,  0}}};
//         s_point d = {{{ 0,  0,  1}}};
//         s_point e_on  = {{{ 0, -1,  0}}};
//         s_point e_in  = {{{ 0,  0,  0}}};
//         s_point e_out = {{{ 2,  0,  0}}};
//         printf("on sphere:  %d (should be  0)\n",
//             insphere(a.x,a.y,a.z, b.x,b.y,b.z, c.x,c.y,c.z, d.x,d.y,d.z,
//                      e_on.x,e_on.y,e_on.z));
//         printf("inside:     %d (should be  1)\n",
//             insphere(a.x,a.y,a.z, b.x,b.y,b.z, c.x,c.y,c.z, d.x,d.y,d.z,
//                      e_in.x,e_in.y,e_in.z));
//         printf("outside:    %d (should be -1)\n",
//             insphere(a.x,a.y,a.z, b.x,b.y,b.z, c.x,c.y,c.z, d.x,d.y,d.z,
//                      e_out.x,e_out.y,e_out.z));
//     }
//
//     // ==========================================================================
//     // POWERTEST1D / test_orthosegment
//     // ==========================================================================
//     {
//         printf("\n=== powertest1d / test_orthosegment ===\n");
//
//         // Segment xa=0,xb=2 with zero weights.
//         // Orthosegment centered at midpoint 1, radius 1.
//         // Boundary points: xc=0 (=xa) and xc=2 (=xb).
//         // xc=0.5: inside   xc=3: outside   xc=0: on boundary
//         double x[2]  = {0, 2};
//         double wx[2] = {0, 0};
//
//         // raw predicate: +1 means outside, -1 means inside (for xa < xb)
//         printf("--- raw powertest1d ---\n");
//         printf("inside:      %d (should be -1)\n", powertest1d(x[0],wx[0], x[1],wx[1], 0.5, 0));
//         printf("outside:     %d (should be +1)\n", powertest1d(x[0],wx[0], x[1],wx[1], 3,   0));
//         printf("on boundary: %d (should be  0)\n", powertest1d(x[0],wx[0], x[1],wx[1], 0,   0));
//         printf("wq=+10:      %d (should be -1)\n", powertest1d(x[0],wx[0], x[1],wx[1], 3,  10));
//         printf("wq=-10:      %d (should be +1)\n", powertest1d(x[0],wx[0], x[1],wx[1], 0.5,-10));
//
//         // wrapper: +1 means inside, -1 means outside
//         printf("--- test_orthosegment ---\n");
//         printf("inside:      %d (should be  1)\n", test_orthosegment(x, wx, 0.5, 0));
//         printf("outside:     %d (should be -1)\n", test_orthosegment(x, wx, 3,   0));
//         printf("on boundary: %d (should be  0)\n", test_orthosegment(x, wx, 0,   0));
//         printf("wq=+10:      %d (should be  1)\n", test_orthosegment(x, wx, 3,  10));
//         printf("wq=-10:      %d (should be -1)\n", test_orthosegment(x, wx, 0.5,-10));
//
//         // reversed segment — same geometric result due to orientation handling
//         double x_rev[2]  = {2, 0};
//         double wx_rev[2] = {0, 0};
//         printf("--- reversed segment ---\n");
//         printf("inside:      %d (should be  1)\n", test_orthosegment(x_rev, wx_rev, 0.5, 0));
//         printf("outside:     %d (should be -1)\n", test_orthosegment(x_rev, wx_rev, 3,   0));
//         printf("on boundary: %d (should be  0)\n", test_orthosegment(x_rev, wx_rev, 0,   0));
//
//         // Verify test_orthosegment agrees with test_orthocircle
//         // for a 1D slice: segment from (-1,0) to (1,0) on x-axis
//         double seg[2]    = {-1, 1};
//         double wseg[2]   = {0, 0};
//         printf("1d slice inside:  %d (should be  1)\n", test_orthosegment(seg, wseg,  0,   0));
//         printf("1d slice outside: %d (should be -1)\n", test_orthosegment(seg, wseg,  2,   0));
//         printf("1d slice on:      %d (should be  0)\n", test_orthosegment(seg, wseg, -1,   0));
//     }
//
//     // ==========================================================================
//     // POWERTEST2D / test_orthocircle
//     // ==========================================================================
//     {
//         printf("\n=== powertest2d / test_orthocircle ===\n");
//         // unit circle, zero weights
//         double a[2] = { 1,  0};
//         double b[2] = { 0,  1};
//         double c[2] = {-1,  0};
//         double q_on[2]  = { 0, -1};
//         double q_in[2]  = { 0,  0};
//         double q_out[2] = { 2,  0};
//         printf("on circle:  %d (should be  0)\n",
//             test_orthocircle(a,0, b,0, c,0, q_on,0));
//         printf("inside:     %d (should be  1)\n",
//             test_orthocircle(a,0, b,0, c,0, q_in,0));
//         printf("outside:    %d (should be -1)\n",
//             test_orthocircle(a,0, b,0, c,0, q_out,0));
//         // weight pulls q_on inside
//         printf("on+wq=1:    %d (should be  1)\n",
//             test_orthocircle(a,0, b,0, c,0, q_on,1));
//         printf("on+wq=-1:   %d (should be -1)\n",
//             test_orthocircle(a,0, b,0, c,0, q_on,-1));
//         // agree with incircle when weights zero
//         printf("agree incircle inside:  %d vs %d (should match)\n",
//             test_orthocircle(a,0, b,0, c,0, q_in,0),
//             test_incircle(a, b, c, q_in));
//         printf("agree incircle outside: %d vs %d (should match)\n",
//             test_orthocircle(a,0, b,0, c,0, q_out,0),
//             test_incircle(a, b, c, q_out));
//     }
//
//     // ==========================================================================
//     // POWERTEST3D / test_orthosphere
//     // ==========================================================================
//     {
//         printf("\n=== powertest3d / test_orthosphere ===\n");
//         s_point a = {{{ 1,  0,  0}}};
//         s_point b = {{{-1,  0,  0}}};
//         s_point c = {{{ 0,  1,  0}}};
//         s_point d = {{{ 0,  0,  1}}};
//         double w[] = {0, 0, 0, 0};
//         s_point q_on  = {{{ 0, -1,  0}}};
//         s_point q_in  = {{{ 0,  0,  0}}};
//         s_point q_out = {{{ 2,  0,  0}}};
//         printf("on sphere:  %d (should be  0)\n",
//             test_orthosphere((s_point[]){a,b,c,d}, w, q_on,  0));
//         printf("inside:     %d (should be  1)\n",
//             test_orthosphere((s_point[]){a,b,c,d}, w, q_in,  0));
//         printf("outside:    %d (should be -1)\n",
//             test_orthosphere((s_point[]){a,b,c,d}, w, q_out, 0));
//         // weight pulls q_on inside/outside
//         printf("on+wq=+1:   %d (should be  1)\n",
//             test_orthosphere((s_point[]){a,b,c,d}, w, q_on,  1));
//         printf("on+wq=-1:   %d (should be -1)\n",
//             test_orthosphere((s_point[]){a,b,c,d}, w, q_on, -1));
//         // agree with insphere when weights zero
//         printf("agree insphere inside:  %d vs %d (should match)\n",
//             test_orthosphere((s_point[]){a,b,c,d}, w, q_in,  0),
//             test_insphere((s_point[]){a,b,c,d}, q_in));
//         printf("agree insphere outside: %d vs %d (should match)\n",
//             test_orthosphere((s_point[]){a,b,c,d}, w, q_out, 0),
//             test_insphere((s_point[]){a,b,c,d}, q_out));
//         printf("agree insphere on:      %d vs %d (should match)\n",
//             test_orthosphere((s_point[]){a,b,c,d}, w, q_on,  0),
//             test_insphere((s_point[]){a,b,c,d}, q_on));
//     }
//
//
//     // // RAW INSPHERE
//     // {
//     // s_point pa = {{{1,  0,  0}}};
//     // s_point pb = {{{-1, 0,  0}}};
//     // s_point pc = {{{0,  1,  0}}};
//     // s_point pd = {{{0,  0,  1}}};  // ← z != 0, breaks coplanarity
//     // s_point pe_on  = {{{0, -1,  0}}};  // on sphere, should give 0
//     // s_point pe_in  = {{{0,  0,  0}}};  // inside, should give nonzero
//     // s_point pe_out = {{{2,  0,  0}}};  // outside, should give nonzero
//     //
//     // printf("orientation: %d\n",
//     //     orient3d(pa.x,pa.y,pa.z, pb.x,pb.y,pb.z,
//     //              pc.x,pc.y,pc.z, pd.x,pd.y,pd.z));
//     // printf("on sphere:   %d (should be 0)\n",
//     //     insphere(pa.x,pa.y,pa.z, pb.x,pb.y,pb.z,
//     //              pc.x,pc.y,pc.z, pd.x,pd.y,pd.z,
//     //              pe_on.x,pe_on.y,pe_on.z));
//     // printf("inside:      %d\n",
//     //     insphere(pa.x,pa.y,pa.z, pb.x,pb.y,pb.z,
//     //              pc.x,pc.y,pc.z, pd.x,pd.y,pd.z,
//     //              pe_in.x,pe_in.y,pe_in.z));
//     // printf("outside:     %d\n",
//     //     insphere(pa.x,pa.y,pa.z, pb.x,pb.y,pb.z,
//     //              pc.x,pc.y,pc.z, pd.x,pd.y,pd.z,
//     //              pe_out.x,pe_out.y,pe_out.z));
//     // printf("powertest3d inside: %d\n",
//     //     powertest3d(pa.x,pa.y,pa.z,0, pb.x,pb.y,pb.z,0,
//     //                 pc.x,pc.y,pc.z,0, pd.x,pd.y,pd.z,0,
//     //                 pe_in.x,pe_in.y,pe_in.z,0));
//     // }
//     //
//     //
//     // // INSPHERE_WEIGHTED
//     // {
//     // s_point pa = {{{1,  0,  0}}};
//     // s_point pb = {{{-1, 0,  0}}};
//     // s_point pc = {{{0,  1,  0}}};
//     // s_point pd = {{{0,  0,  1}}};
//     // double w[] = {0, 0, 0, 0};
//     // s_point pe_on  = {{{0, -1,  0}}};  // on unit sphere
//     // s_point pe_in  = {{{0,  0,  0}}};  // inside
//     // s_point pe_out = {{{2,  0,  0}}};  // outside
//     //
//     // printf("Test exact integer points, wi=0\n");
//     // printf("on sphere: %d (should be 0)\n",
//     //     test_orthosphere((s_point[]){pa,pb,pc,pd}, w, pe_on, 0));
//     // printf("inside:    %d (should be 1)\n",
//     //     test_orthosphere((s_point[]){pa,pb,pc,pd}, w, pe_in, 0));
//     // printf("outside:   %d (should be -1)\n",
//     //     test_orthosphere((s_point[]){pa,pb,pc,pd}, w, pe_out, 0));
//     //
//     // // now test weights: pe_on has distance 0, adding wq moves it in or out
//     // printf("on sphere wq=+1: %d (should be 1, weight pulls in)\n",
//     //     test_orthosphere((s_point[]){pa,pb,pc,pd}, w, pe_on, 1));
//     // printf("on sphere wq=-1: %d (should be -1, weight pushes out)\n",
//     //     test_orthosphere((s_point[]){pa,pb,pc,pd}, w, pe_on, -1));
//     // }
//     //
//
//
//     puts("\n\n\n");
//     s_point triangle[3] = { {{{0, -1, 0}}},
//                             {{{1, 1, 1}}},
//                             {{{-1, 1, 1}}} };
//     s_point segment1[2] = { {{{0, 0, -2}}},
//                             {{{0, 0, 2}}} };
//     s_point segment2[2] = { {{{0, 0, 2}}},
//                             {{{0, 0, -2}}} };
//
//
//     printf("Segment crosses triangle: %d (should be 1)\n", test_segment_triangle_intersect_3D(segment1, triangle, 0, 0) == INTERSECT_NONDEGENERATE);
//     printf("Segment crosses triangle: %d (should be 1)\n", test_segment_triangle_intersect_3D(segment2, triangle, 0, 0) == INTERSECT_NONDEGENERATE);
//
//     s_point triangle2[3] = { {{{-1, 0, 0}}},
//                              {{{1, 0, 0}}},
//                              {{{0, 0, 1}}} };
//     s_point p1 = {{{0.5, 0, 0}}};
//     s_point p2 = {{{2, 0, 0}}};
//     printf("point in boundary of triangle: %d (Should be 1)\n", test_point_in_triangle_3D(triangle2, p1, 0, 0) == TEST_BOUNDARY);
//     printf("point outside triangle: %d (Should be 1)\n", test_point_in_triangle_3D(triangle2, p2, 0, 0) == TEST_OUT);
//
//     
//
// }
