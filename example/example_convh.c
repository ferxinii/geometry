
#include "points.h"
#include "convh.h"
#include <stdio.h>

double EPS_degenerate = 1e-12;

int main(void) {

    /* Standard cube */
    s_point cube[8] = {
            {{{0, 0, 0}}},
            {{{1, 0, 0}}},
            {{{0, 1, 0}}},
            {{{1, 1, 0}}},
            {{{0, 0, 1}}},
            {{{1, 0, 1}}},
            {{{0, 1, 1}}},
            {{{1, 1, 1}}}
    };
    s_points p_cube = {8, cube};

    s_convh ch_cube;
    if (convhull_from_points(&p_cube, EPS_degenerate, 0, &ch_cube) != 1) {
        printf("Error making convhull of degenerate cube.\n");
        return 0;
    }
    printf("CUBE: Volume = %g, area = %g\n", volume_convhull(&ch_cube), surface_area_convhull(&ch_cube));

    write_convhull_to_m(&ch_cube, "cube.m");


    /* Duplicate vertices */
    puts("");
    s_point cube_duplicates[12] = {
            {{{0, 0, 0}}},
            {{{1, 0, 0}}},
            {{{0, 1, 0}}},
            {{{1, 1, 0}}},
            {{{0, 0, 1}}},
            {{{1, 0, 1}}},
            {{{0, 1, 1}}},
            {{{1, 1, 1}}},
            /* Duplicates of 0, 2, 4, 6 */
            {{{0, 0, 0}}},
            {{{0, 1, 0}}},
            {{{0, 0, 1}}},
            {{{0, 1, 1}}},
    };
    s_points p_cube_duplicates = {12, cube_duplicates};

    s_convh ch_cube_duplicates;
    if (convhull_from_points(&p_cube_duplicates, EPS_degenerate, 0, &ch_cube_duplicates) != 1) {
        printf("Error making convhull of degenerate cube.\n");
        return 0;
    }
    printf("CUBE DUPLICATES: Volume = %g, area = %g\n", volume_convhull(&ch_cube_duplicates), surface_area_convhull(&ch_cube_duplicates));

    write_convhull_to_m(&ch_cube_duplicates, "cube_duplicates.m");


    /* Degeneracies */
    puts("");
    s_point cube_degenerate[23] = { 
            /* --- Original cube vertices (8) --- */
            {{{0, 0, 0}}},
            {{{1, 0, 0}}},
            {{{0, 1, 0}}},
            {{{1, 1, 0}}},
            {{{0, 0, 1}}},
            {{{1, 0, 1}}},
            {{{0, 1, 1}}},
            {{{1, 1, 1}}},
            /* --- Duplicates of some vertices (4) --- */
            {{{0, 0, 0}}},   /* duplicate of vertex 0 */
            {{{1, 1, 1}}},   /* duplicate of vertex 7 */
            {{{1, 0, 0}}},   /* duplicate of vertex 1 */
            {{{0, 1, 0}}},   /* duplicate of vertex 2 */
            /* --- Colinear along X axis (4) --- */
            {{{0.25, 0, 0}}},
            {{{0.4,  0, 0}}},
            {{{0.55,  0, 0}}},
            {{{0.75,  0, 0}}},
            /* --- Colinear along the main diagonal (3) --- */
            {{{0.25, 0.25, 0.25}}},
            {{{0.5, 0.5, 0.5}}},
            {{{0.6, 0.6, 0.6}}},
            /* --- Coplanar set: z = 0 plane (4) --- */
            {{{0.2, 0.8, 0}}},
            {{{0.8, 0.2, 0}}},
            {{{0.3, 0.6, 0}}},
            {{{0.6, 0.3, 0}}}, 
    };
    s_points p_cube_degenerate = {23, cube_degenerate};

    s_convh ch_cube_degenerate;
    if (convhull_from_points(&p_cube_degenerate, EPS_degenerate, 0, &ch_cube_degenerate) != 1) {
        printf("Error making convhull of degenerate cube.\n");
        return 0;
    }
    printf("CUBE DEGENERATE: Volume = %g, area = %g\n", volume_convhull(&ch_cube_degenerate), surface_area_convhull(&ch_cube_degenerate));

    write_convhull_to_m(&ch_cube_degenerate, "cube_degenerate.m");


    /* Noise */
    puts("");
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
    puts("");
    s_convh S1; 
    if (convhull_from_csv("lobes/S1.csv", EPS_degenerate, 0, &S1) != 1) {
        printf("Error making convhull of lung segment.\n");
        return 0;
    }
    printf("LUNG SEGMENT: Volume = %g, area = %g\n", volume_convhull(&S1), surface_area_convhull(&S1));
    write_convhull_to_m(&S1, "lung_segment.m");


    // /* Lung */
    // puts("");
    // s_convh lung; 
    // if (convhull_from_csv("lobes/lung.csv", EPS_degenerate, 0, &lung) != 1) {
    //     printf("Error making convhull of lung segment.\n");
    //     return 0;
    // }
    // printf("LUNG SEGMENT: Volume = %g, area = %g\n", volume_convhull(&lung), surface_area_convhull(&lung));
    // write_convhull_to_m(&lung, "lung_segment.m");
    //

    /* Non-manifold */
    puts("");
    s_convh nm; 
    if (convhull_from_csv("../problematic.csv", EPS_degenerate, 0, &nm) != 1) {
        printf("Error making convhull of nonmanifold.\n");
        return 0;
    }
    printf("nonmanifold: Volume = %g, area = %g\n", volume_convhull(&nm), surface_area_convhull(&nm));
    write_convhull_to_m(&nm, "nm.m");


    /* Problematic 2: adding collinear faces */
    puts("");
    s_convh nm2; 
    if (convhull_from_csv("../problematic_2.csv", 1e-14, 1e-14, &nm2) != 1) {
        printf("Error making convhull of nonmanifold.\n");
        return 0;
    }
    printf("nonmanifold2: Volume = %g, area = %g\n", volume_convhull(&nm2), surface_area_convhull(&nm2));
    write_convhull_to_m(&nm2, "nm2.m");

    /* Problematic 3: adding collinear faces */
    puts("");
    s_convh nm3; 
    if (convhull_from_csv("../problematic_3.csv", 1e-14, 1e-14, &nm3) != 1) {
        printf("Error making convhull of nonmanifold.\n");
        return 0;
    }
    printf("nonmanifold3: Volume = %g, area = %g\n", volume_convhull(&nm3), surface_area_convhull(&nm3));
    write_convhull_to_m(&nm3, "nm3.m");

    /* Problematic 4: adding collinear faces */
    puts("");
    s_convh nm4; 
    if (convhull_from_csv("../problematic_4.csv", 1e-14, 1e-14, &nm4) != 1) {
        printf("Error making convhull of nonmanifold.\n");
        return 0;
    }
    printf("nonmanifold4: Volume = %g, area = %g\n", volume_convhull(&nm4), surface_area_convhull(&nm4));
    write_convhull_to_m(&nm4, "nm4.m");

    /* Problematic 5: adding collinear faces */
    puts("");
    s_convh nm5; 
    if (convhull_from_csv("../problematic_5.csv", 1e-14, 1e-14, &nm5) != 1) {
        printf("Error making convhull of nonmanifold.\n");
        return 0;
    }
    printf("nonmanifold5: Volume = %g, area = %g\n", volume_convhull(&nm5), surface_area_convhull(&nm5));
    write_convhull_to_m(&nm5, "nm5.m");

    /* Problematic 6: adding collinear faces */
    puts("");
    s_convh nm6; 
    if (convhull_from_csv("../problematic_6.csv", 1e-14, 1e-14, &nm6) != 1) {
        printf("Error making convhull of nonmanifold.\n");
        return 0;
    }
    printf("nonmanifold6: Volume = %g, area = %g\n", volume_convhull(&nm6), surface_area_convhull(&nm6));
    write_convhull_to_m(&nm6, "nm6.m");

    /* Problematic 7*/
    puts("");
    s_convh nm7; 
    if (convhull_from_csv("../problematic_7.csv", 1e-14, 1e-14, &nm7) != 1) {
        printf("Error making convhull of nonmanifold.\n");
        return 0;
    }
    printf("nonmanifold7: Volume = %g, area = %g\n", volume_convhull(&nm7), surface_area_convhull(&nm7));
    write_convhull_to_m(&nm7, "nm7.m");

    /* Overlapping faces */
    puts("");
    s_convh of; 
    if (convhull_from_csv("../overlapping_faces.csv", 1e-14, 1e-14, &of) != 1) {
        printf("Error making convhull of overlapping_faces.\n");
        return 0;
    }
    printf("overlapping_faces: Volume = %g, area = %g\n", volume_convhull(&of), surface_area_convhull(&of));
    write_convhull_to_m(&of, "of.m");

    /* Problematic 8 */
    puts("");
    s_convh nm8; 
    if (convhull_from_csv("../problematic_8.csv", 1e-14, 1e-14, &nm8) != 1) {
        printf("Error making convhull of nm8.\n");
        return 0;
    }
    printf("overlapping_faces: Volume = %g, area = %g\n", volume_convhull(&nm8), surface_area_convhull(&nm8));
    write_convhull_to_m(&nm8, "nm8.m");

    /* Problematic 9 */
    puts("");
    s_convh nm9; 
    if (convhull_from_csv("../problematic_9.csv", 1e-14, 1e-14, &nm9) != 1) {
        printf("Error making convhull of nm9.\n");
        return 0;
    }
    printf("overlapping_faces: Volume = %g, area = %g\n", volume_convhull(&nm9), surface_area_convhull(&nm9));
    write_convhull_to_m(&nm9, "nm9.m");
    
    /* Overlapping faces 2*/
    puts("");
    s_convh of2; 
    if (convhull_from_csv("../of_2.csv", 1e-14, 1e-14, &of2) != 1) {
        printf("Error making convhull of overlapping_faces_2.\n");
        return 0;
    }
    printf("overlapping_faces: Volume = %g, area = %g\n", volume_convhull(&of2), surface_area_convhull(&of2));
    write_convhull_to_m(&of2, "of2.m");



}
