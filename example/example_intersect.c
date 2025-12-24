
#include "points.h"
#include "convh.h"
#include "convh.h"
#include <stdlib.h>
#include <stdio.h>


int main() 
{
    double TOL = 1e-12, EPS_degenerate = 1e-12, min_vol_I = 1e-9;
    // double TOL = 0, EPS_degenerate = 0, min_vol_I = 0;
    /* CUBES */
    s_point cube1[8] = { {{{0, 0, 0}}},
                         {{{1, 0, 0}}},
                         {{{0, 1, 0}}},
                         {{{1, 1, 0}}},
                         {{{0, 0, 1}}},
                         {{{1, 0, 1}}},
                         {{{0, 1, 1}}},
                         {{{1, 1, 1}}} };
    s_points p_cube1 = {8, cube1};
    s_convh ch1; convhull_from_points(&p_cube1, EPS_degenerate, TOL, &ch1);
    
    s_point translate = {{{0.5, 0.5, 0.5}}};
    s_point cube2[8] = { sum_points(cube1[0], translate),
                         sum_points(cube1[1], translate),
                         sum_points(cube1[2], translate),
                         sum_points(cube1[3], translate),
                         sum_points(cube1[4], translate),
                         sum_points(cube1[5], translate),
                         sum_points(cube1[6], translate),
                         sum_points(cube1[7], translate) };
    s_points p_cube2 = {8, cube2};
    s_convh ch2; convhull_from_points(&p_cube2, EPS_degenerate, TOL, &ch2);

    s_convh I;
    intersection_convhulls(&ch1, &ch2, EPS_degenerate, TOL, &I);

    write_convhull_to_m(&ch1, "cube1.m");
    write_convhull_to_m(&ch2, "cube2.m");
    write_convhull_to_m(&I, "intersection_cubes.m");
    free_convhull(&ch1);
    free_convhull(&ch2);
    free_convhull(&I);

    /* LUNG LOBES */
    s_points pL2 = read_points_from_csv("lobes/L2.txt");
    s_points pL3 = read_points_from_csv("lobes/L3.txt");
    s_convh L2; convhull_from_points(&pL2, EPS_degenerate, TOL, &L2);
    s_convh L3; convhull_from_points(&pL3, EPS_degenerate, TOL, &L3);
    
    free_convhull(&I);
    intersection_convhulls(&L2, &L3, EPS_degenerate, TOL, &I);

    write_convhull_to_m(&L2, "L2.m");
    write_convhull_to_m(&L3, "L3.m");
    write_convhull_to_m(&I, "I_L2_L3.m");

    remove_intersection_convhulls(&L2, &L3, EPS_degenerate, TOL, min_vol_I);
    write_convhull_to_m(&L2, "L2_post.m");
    write_convhull_to_m(&L3, "L3_post.m");

    free_convhull(&L2);
    free_convhull(&L3);
    free_convhull(&I);
}


