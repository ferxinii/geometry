#ifndef GEOMETRY_GTESTS_H
#define GEOMETRY_GTESTS_H
#include "points.h"

typedef enum geom_test {  
    TEST_IN,
    TEST_OUT,
    TEST_BOUNDARY,
    TEST_DEGENERATE,
    TEST_ERROR
} e_geom_test;

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


/* My wrappers (considers orientation of first points) */
int test_orientation_2d(const s_point2d line[2], s_point2d p);
int test_orientation(const s_point plane[3], s_point p);
int test_incircle(const s_point2d circle[3], s_point2d p);
int test_insphere(const s_point sph[4], s_point p);
int test_orthosegment(int k, const double c[k], const double wc[k], 
                      double xp, double wp);
int test_orthocircle(int k, const s_point2d c[k], const double wc[k], 
                     s_point2d p, double wp);
int test_orthosphere(int k, const s_point c[k], const double wc[k],
                     s_point p, double wp);
int test_orthosegment_w(int k, const double c[k], const double wc[k], 
                        double alpha);
int test_orthocircle_w(int k, const s_point2d c[k], const double wc[k], 
                       double alpha);
int test_orthosphere_w(int k, const s_point c[k], const double wc[k],
                       double alpha);





/* EPS_degenerate is a scale for the minimum value of an object to be non-degenerate. Avoids division by 0, ignore too small / degenerate triangles, ...
 * TOL_boundary is the distance from a point to the object's boundary to be considered as belonging to it. If ==0, tests are ROBUST. */
e_geom_test test_point_in_interval_1D(double x, double a, double b, double EPS_degenerate, double TOL_boundary);
e_geom_test test_point_in_triangle_2D(const s_point2d triangle[3], s_point2d p, double EPS_degenerate, double TOL_boundary);
e_geom_test test_point_in_triangle_3D(const s_point triangle[3], s_point p, double EPS_degenerate, double TOL_boundary);
e_geom_test test_point_in_tetrahedron(const s_point tetra[4], s_point query, double EPS_degenerate, double TOL_boundary);
s_points_test test_points_in_halfspace(const s_point plane_ordered[3], const s_points *points, double EPS_degenerate, double TOL_boundary, e_geom_test out_buff[points->N]);
/* Segment intersections */
e_intersect_type test_segment_segment_intersect_2D(const s_point2d s1[2], const s_point2d s2[2], double EPS_degenerate, double TOL_boundary);
e_intersect_type test_segment_plane_intersect(const s_point seg[2], const s_point plane[3], double EPS_degenerate, double TOL_boundary);
e_intersect_type test_segment_triangle_intersect_2D(const s_point2d s[2], const s_point2d tri[3], double EPS_degenerate, double TOL_boundary);
e_intersect_type test_segment_triangle_intersect_3D(const s_point segment[2], const s_point triangle[3], double EPS_degenerate, double TOL_boundary);
s_segment_intersect segment_segment_intersect_2D(const s_point2d s1[2], const s_point2d s2[2], double EPS_degenerate, double TOL_boundary);
s_segment_intersect segment_plane_intersect(const s_point seg[2], const s_point plane[3], double EPS_degenerate, double TOL_boundary);
s_segment_intersect segment_triangle_intersect_2D(const s_point2d s[2], const s_point2d tri[3], double EPS_degenerate, double TOL_boundary);
s_segment_intersect segment_triangle_intersect_3D(const s_point segment[2], const s_point triangle[3], double EPS_degenerate, double TOL_boundary);

#endif
