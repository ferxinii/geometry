#ifndef CONVH_H
#define CONVH_H

#include "geometry.h"


typedef struct convhull {
    s_points points;
    int Nf;
    int *faces;
    s_point *fnormals;  // Unnormalized
} s_convhull;


s_convhull convhull_from_points(const s_points *points); // is all 0 / NULL if error 
void free_convhull(s_convhull *convh);
s_convhull copy_convhull(const s_convhull *in);

int is_inside_convhull(const s_convhull *convh, s_point query);  // 1: inside, 0: outise, -1: in boundary
int is_in_boundary_convhull(const s_convhull *convh, int point_id);
int mark_boundary_convhull(const s_convhull *convh, int out[convh->points.N]);
s_points boundary_convhull(const s_convhull *convh);
int points_inside_convhull(const s_convhull *convh, const s_points query, int out_mark[query.N]);

s_point random_point_inside_convhull(const s_convhull *convh, s_point min, s_point max);

double volume_convhull(const s_convhull *convh);
double volume_convhull_from_points(const s_points *points);


void convh_get_face(const s_convhull *convh, int id, s_point out[3]);
int mark_faces_incident_to_vertex(const s_convhull *C, int vid, int out[C->Nf]);

s_convhull intersection_convhulls(const s_convhull *A, const s_convhull *B);



#endif
