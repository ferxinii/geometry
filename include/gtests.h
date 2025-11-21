#ifndef GEOMETRY_GTESTS_H
#define GEOMETRY_GTESTS_H
#include "points.h"

typedef enum geom_test {  
    TEST_IN,
    TEST_OUT,
    TEST_BOUNDARY,
    TEST_DEGENERATE,
    TEST_ERROR
}   e_geom_test;

typedef struct points_test {
    int Nin;
    int Nbdy;
    int Nout;
    int Nerr;
    e_geom_test *indicator;
} s_points_test;

typedef enum intersect_type {
    INTERSECT_EMPTY,
    INTERSECT_NONDEGENERATE,
    INTERSECT_DEGENERATE,
    INTERSECT_ERROR
} e_intersect_type;

typedef struct segment_intersect {
    int N;
    s_point coords[2];
    e_intersect_type type;
} s_segment_intersect;


/* My wrappers, +1, 0, -1 */
int orientation_robust(const s_point p[3], s_point q);
int insphere_robust(const s_point p[4], s_point q);
/* Shewchuck's robust geometric predicates */
extern void exactinit(void);  
extern double orient2d(const double *pa, const double *pb, const double *pc);
extern double orient3d(const double *pa, const double *pb, const double *pc, const double *pd);
extern double incircle(const double *pa, const double *pb, const double *pc, const double *pd);
extern double insphere(const double *pa, const double *pb, const double *pc, const double *pd, const double *pe);

/* EPS_degenerate is a scale for the minimum value of an object to be non-degenerate. Avoids division by 0, ignore too small / degenerate triangles, ...
 * TOL_boundary is the distance from a point to the object's boundary to be considered as belonging to it. If ==0, tests are ROBUST. */
e_geom_test test_point_in_interval_1D(double x, double a, double b, double EPS_degenerate, double TOL_boundary);
e_geom_test test_point_in_triangle_2D(const double a[2], const double b[2], const double c[2], const double p[2], double EPS_degenerate, double TOL_boundary);
e_geom_test test_point_in_triangle_3D(const s_point triangle[3], s_point p, double EPS_degenerate, double TOL_boundary);
e_geom_test test_point_in_tetrahedron(const s_point tetra[4], s_point query, double EPS_degenerate, double TOL_boundary);
s_points_test test_points_in_halfspace(const s_point plane_ordered[3], const s_points *points, double EPS_degenerate, double TOL_boundary, e_geom_test out_buff[points->N]);

/* Segment intersections */
e_intersect_type test_segment_segment_intersect_2D(const double A1[2], const double A2[2], const double B1[2], const double B2[2], double EPS_degenerate, double TOL_boundary);
e_intersect_type test_segment_plane_intersect(const s_point seg[2], const s_point plane[3], double EPS_degenerate, double TOL_boundary);
e_intersect_type test_segment_triangle_intersect_2D(const double S1[2], const double S2[2], const double A[2], const double B[2], const double C[2], double EPS_degenerate, double TOL_boundary);
e_intersect_type test_segment_triangle_intersect_3D(const s_point segment[2], const s_point triangle[3], double EPS_degenerate, double TOL_boundary);

s_segment_intersect segment_segment_intersect_2D(const double A1[2], const double A2[2], const double B1[2], const double B2[2], double EPS_degenerate, double TOL_boundary);
s_segment_intersect segment_plane_intersect(const s_point seg[2], const s_point plane[3], double EPS_degenerate, double TOL_boundary);
s_segment_intersect segment_triangle_intersect_2D(const double S1[2], const double S2[2], const double A[2], const double B[2], const double C[2], double EPS_degenerate, double TOL_boundary); s_segment_intersect segment_triangle_intersect_3D(const s_point segment[2], const s_point triangle[3], double EPS_degenerate, double TOL_boundary);


#endif
