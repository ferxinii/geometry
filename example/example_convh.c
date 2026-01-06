
#include "points.h"
#include "convh.h"
#include <stdio.h>

double EPS_degenerate = 1e-12;

int main(void) {

    /* Standard cube */
    // s_point cube[8] = {
    //         {{{0, 0, 0}}},
    //         {{{1, 0, 0}}},
    //         {{{0, 1, 0}}},
    //         {{{1, 1, 0}}},
    //         {{{0, 0, 1}}},
    //         {{{1, 0, 1}}},
    //         {{{0, 1, 1}}},
    //         {{{1, 1, 1}}}
    // };
    // s_points p_cube = {8, cube};
    //
    // s_convh ch_cube;
    // if (convhull_from_points(&p_cube, EPS_degenerate, 0, &ch_cube) != 1) {
    //     printf("Error making convhull of degenerate cube.\n");
    //     return 0;
    // }
    // printf("CUBE: Volume = %g, area = %g\n", volume_convhull(&ch_cube), surface_area_convhull(&ch_cube));
    //
    // write_convhull_to_m(&ch_cube, "cube.m");
    //
    // 
    // /* Duplicate vertices */
    // s_point cube_duplicates[12] = {
    //         {{{0, 0, 0}}},
    //         {{{1, 0, 0}}},
    //         {{{0, 1, 0}}},
    //         {{{1, 1, 0}}},
    //         {{{0, 0, 1}}},
    //         {{{1, 0, 1}}},
    //         {{{0, 1, 1}}},
    //         {{{1, 1, 1}}},
    //         /* Duplicates of 0, 2, 4, 6 */
    //         {{{0, 0, 0}}},
    //         {{{0, 1, 0}}},
    //         {{{0, 0, 1}}},
    //         {{{0, 1, 1}}},
    // };
    // s_points p_cube_duplicates = {12, cube_duplicates};
    //
    // s_convh ch_cube_duplicates;
    // if (convhull_from_points(&p_cube_duplicates, EPS_degenerate, 0, &ch_cube_duplicates) != 1) {
    //     printf("Error making convhull of degenerate cube.\n");
    //     return 0;
    // }
    // printf("CUBE DUPLICATES: Volume = %g, area = %g\n", volume_convhull(&ch_cube_duplicates), surface_area_convhull(&ch_cube_duplicates));
    //
    // write_convhull_to_m(&ch_cube_duplicates, "cube_duplicates.m");
    // 
    //
    // /* Degeneracies */
    // s_point cube_degenerate[23] = { 
    //         /* --- Original cube vertices (8) --- */
    //         {{{0, 0, 0}}},
    //         {{{1, 0, 0}}},
    //         {{{0, 1, 0}}},
    //         {{{1, 1, 0}}},
    //         {{{0, 0, 1}}},
    //         {{{1, 0, 1}}},
    //         {{{0, 1, 1}}},
    //         {{{1, 1, 1}}},
    //         /* --- Duplicates of some vertices (4) --- */
    //         {{{0, 0, 0}}},   /* duplicate of vertex 0 */
    //         {{{1, 1, 1}}},   /* duplicate of vertex 7 */
    //         {{{1, 0, 0}}},   /* duplicate of vertex 1 */
    //         {{{0, 1, 0}}},   /* duplicate of vertex 2 */
    //         /* --- Colinear along X axis (4) --- */
    //         {{{0.25, 0, 0}}},
    //         {{{0.4,  0, 0}}},
    //         {{{0.55,  0, 0}}},
    //         {{{0.75,  0, 0}}},
    //         /* --- Colinear along the main diagonal (3) --- */
    //         {{{0.25, 0.25, 0.25}}},
    //         {{{0.5, 0.5, 0.5}}},
    //         {{{0.6, 0.6, 0.6}}},
    //         /* --- Coplanar set: z = 0 plane (4) --- */
    //         {{{0.2, 0.8, 0}}},
    //         {{{0.8, 0.2, 0}}},
    //         {{{0.3, 0.6, 0}}},
    //         {{{0.6, 0.3, 0}}}, 
    // };
    // s_points p_cube_degenerate = {23, cube_degenerate};
    //
    // s_convh ch_cube_degenerate;
    // if (convhull_from_points(&p_cube_degenerate, EPS_degenerate, 0, &ch_cube_degenerate) != 1) {
    //     printf("Error making convhull of degenerate cube.\n");
    //     return 0;
    // }
    // printf("CUBE DEGENERATE: Volume = %g, area = %g\n", volume_convhull(&ch_cube_degenerate), surface_area_convhull(&ch_cube_degenerate));
    //
    // write_convhull_to_m(&ch_cube_degenerate, "cube_degenerate.m");
    //
    
    /* Noise */
    s_point cube_noise[12] = {
            {{{0, 0, 0}}},
            {{{1, 0, 0}}},
            {{{0, 1, 0}}},
            {{{1, 1, 0}}},
            {{{0, 0, 1}}},
            {{{1, 0, 1}}},
            {{{0, 1, 1}}},
            {{{1, 1, 1}}},
            /* Duplicates of 0, 2, 4, 6 */
            {{{-1e-8, 0, 0}}},
            {{{0, 1, -1e-10}}},
            {{{0, -1e-9, 1}}},
            {{{0, 1 +1e-10, 1+1e-9}}},
    };
    s_points p_cube_noise = {12, cube_noise};

    s_convh ch_cube_noise;
    if (convhull_from_points(&p_cube_noise, EPS_degenerate, 0, &ch_cube_noise) != 1) {
        printf("Error making convhull of noisy cube.\n");
        return 0;
    }
    printf("CUBE NOISY: Volume = %g, area = %g\n", volume_convhull(&ch_cube_noise), surface_area_convhull(&ch_cube_noise));

    write_convhull_to_m(&ch_cube_noise, "cube_noise.m");


    /* Segment */
    s_convh S1; 

    if (convhull_from_csv("lobes/S1.csv", EPS_degenerate, 0, &S1) != 1) {
        printf("Error making convhull of lung segment.\n");
        return 0;
    }
    printf("LUNG SEGMENT: Volume = %g, area = %g\n", volume_convhull(&S1), surface_area_convhull(&S1));

    write_convhull_to_m(&S1, "lung_segment.m");


    /* Non-manifold */
    s_convh nm; 

    if (convhull_from_csv("../problematic.csv", EPS_degenerate, 0, &nm) != 1) {
        printf("Error making convhull of nonmanifold.\n");
        return 0;
    }
    printf("nonmanifold: Volume = %g, area = %g\n", volume_convhull(&nm), surface_area_convhull(&nm));

    write_convhull_to_m(&nm, "nm.m");
}
