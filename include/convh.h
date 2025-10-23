#ifndef CONVH_H
#define CONVH_H

#include "geometry.h"

typedef struct convhull {
    const_s_points points;
    int Nf;
    int *faces;
    s_point *fnormals;  // Unnormalized
} s_convhull;


s_convhull convhull_from_points(const_s_points points); 
void free_convhull(s_convhull *convh);

int is_inside_convhull(s_convhull convh, s_point query);  // 1: inside, 0: outise, -1: in boundary
int is_in_boundary_convhull(s_convhull convh, int point_id);
int mark_inside_convhull(s_convhull convh, const_s_points query, int *out_mark);

s_point random_point_inside_convhull(s_convhull convh, s_point min, s_point max);

double volume_convhull(s_convhull convh);
double volume_convhull_from_points(const_s_points points);

#endif
