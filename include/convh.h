#ifndef CONVH_H
#define CONVH_H

#include "points.h"
#include "gtests.h"
#include <stdlib.h>

typedef struct convhull {
    s_points points;
    int Nf;
    int *faces;
    s_point *fnormals;  // Unnormalized
} s_convh;


#define convhull_NAN (s_convh){0}
int convhull_is_valid(const s_convh *convh);

int convhull_from_points(const s_points *points, double EPS_degenerate_face, double TOL, s_convh *out); 
int convhull_from_csv(const char *filename, double EPS_degenerate_face, double TOL, s_convh *out);
/* 1: OK, 0: Error related with nature of points, -1: Error of memory or other serious error. */
void free_convhull(s_convh *convh);
s_convh copy_convhull(const s_convh *in);

int serialize_convhull(const s_convh *convh, uint8_t *buff_write, size_t *size, uint8_t **data);
int deserialize_convhull(const uint8_t *data, s_convh *out, size_t *bytes_read);

typedef struct list s_list;  /* Definition in lists.h */
int list_edges_convhull(const s_convh *C, s_list *out_edges);

e_geom_test test_point_in_convhull(const s_convh *convh, s_point query, double EPS_degenerate, double TOL_boundary);  
s_points_test test_points_in_convhull(const s_convh *convh, const s_points *query, double EPS_degenerate, double TOL_boundary, e_geom_test buff[query->N]);
int convhull_id_in_boundary(const s_convh *convh, int point_id);
int mark_boundary_convhull(const s_convh *convh, int out[convh->points.N]);
s_points points_boundary_convhull(const s_convh *convh, int buff_isboundary[convh->points.N]);

s_point random_point_inside_convhull(const s_convh *convh, double EPS_degenerate_face, s_point min, s_point max);  /* min max can be point_NAN: computed inside */
s_point closest_point_on_convhull_boundary(const s_convh *convh, s_point query, double EPS_degenerate);



double volume_convhull(const s_convh *convh);
double surface_area_convhull(const s_convh *convh);
s_point convhull_volume_centroid(const s_convh *convh, double EPS_degenerate);
int convex_hull_winding_valid(const s_convh *convh);
void convh_get_face(const s_convh *convh, int id, s_point out[3]);
int mark_faces_incident_to_vertex(const s_convh *C, int vid, int out[C->Nf]);

void write_convhull_to_m(const s_convh *vertices, const char *m_filename);

/* ch_intersect.c */
e_intersect_type test_segment_convhull_surface_intersect(const s_convh *C, const s_point segment[2], double EPS_degenerate, double TOL);
s_segment_intersect segment_convhull_surface_intersect(const s_convh *C, const s_point segment[2], double EPS_degenerate, double TOL);
int clip_convhull_halfspace(const s_convh *C, s_point plane[3], double EPS_degenerate, double TOL, s_convh *out); /* 1: Clipped. 0: Non-clipped. -1: Error */
int intersection_convhulls(const s_convh *A, const s_convh *B, double EPS_degenerate, double TOL, s_convh *out_I);
int remove_intersection_convhulls(s_convh *A, s_convh *B, double EPS_degenerate, double TOL, double min_vol_I); /* 1: OK, 0: Error points, -1: Error serious. */
int clip_convhull_convhull(s_convh *C, const s_convh *clipper, double EPS_degenerate, double TOL, double min_vol_I);






#endif
