#ifndef GEOMETRY_GTESTS_H
#define GEOMETRY_GTESTS_H

#include "points.h"

typedef enum geom_test {  
    IN,
    OUT,
    BOUNDARY,
    ERROR
}   e_geom_test;

typedef struct points_test {
    int Nin;
    int Nbdy;
    int Nout;
    int Nerr;
    e_geom_test *indicator;
} s_points_test;


/* My wrappers, +1, 0, -1 */
int orientation_robust(const s_point p[3], s_point q);
int insphere_robust(const s_point p[4], s_point q);
/* Shewchuck's robust geometric predicates */
extern void exactinit(void);  
extern double orient2d(const double *pa, const double *pb, const double *pc);
extern double orient3d(const double *pa, const double *pb, const double *pc, const double *pd);
extern double incircle(const double *pa, const double *pb, const double *pc, const double *pd);
extern double insphere(const double *pa, const double *pb, const double *pc, const double *pd, const double *pe);


e_geom_test test_point_in_interval_1D(double x, double a, double b, double EPS_degenerate, double TOL);
e_geom_test test_point_in_triangle_2D(const double a[2], const double b[2], const double c[2], const double p[2], double EPS_degenerate, double TOL);
e_geom_test test_point_in_triangle_3D(const s_point triangle[3], s_point p, double EPS_degenerate, double TOL);
e_geom_test test_point_in_tetrahedron(const s_point tetra[4], s_point query, double EPS_degenerate, double TOL);
s_points_test test_points_in_halfspace(const s_point plane_ordered[3], const s_points *points, double EPS_degenerate, double TOL, e_geom_test out_buff[points->N]);
e_geom_test test_segment_segment_intersect_2D_robust(const double A1[2], const double A2[2], const double B1[2], const double B2[2]);
e_geom_test test_segment_triangle_intersect_2D_robust(const double S1[2], const double S2[2], const double A[2], const double B[2], const double C[2]);
e_geom_test test_segment_triangle_intersect_3D_robust(const s_point segment[2], const s_point triangle[3]);

e_geom_test test_segment_triangle_intersect_robust(const s_point segment[2], const s_point triangle[3]);


int segment_segment_intersection_2D(const double A1[2], const double A2[2], const double B1[2], const double B2[2], double EPS_degenerate, double TOL, double out[4]);
int segment_plane_intersection(const s_point seg[2], const s_point plane[3], double EPS_degenerate, double TOL, s_point out[2]);
int segment_triangle_intersection_2D(const double s1[2], const double s2[2], const double a[2], const double b[2], const double c[2], double EPS_degenerate, double TOL, double out[4]);
int segment_triangle_intersection_3D(const s_point segment[2], const s_point triangle[3], double EPS_degenerate, double TOL, s_point out[2]);


#endif
