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
int test_orientation_2d(const double a[2], const double b[2], const double p[2]);
int test_orientation(const s_point plane[3], s_point p);
int test_incircle(const double c1[2], const double c2[2], const double c3[2],
                  const double p[2]);
int test_insphere(const s_point sph[4], s_point q);
int test_orthosegment(const double x[2], double wx[2], 
                      double xp, double wp);
int test_orthocircle(const double c1[2], double wc1,
                     const double c2[2], double wc2,
                     const double c3[2], double wc3,
                     const double p[2], double wp);
int test_orthosphere(const s_point sph[4], const double wsph[4],
                     s_point p, double wp);

/* Robust Geometric Predicates external/robust_predicates */
extern int orient2d(double ax, double ay,
                    double bx, double by,
                    double cx, double cy);
extern int orient3d(double ax, double ay, double az,
                    double bx, double by, double bz,
                    double cx, double cy, double cz,
                    double dx, double dy, double dz);
extern int incircle(double ax, double ay,
                    double bx, double by,
                    double cx, double cy,
                    double dx, double dy);
extern int insphere(double ax, double ay, double az,
                    double bx, double by, double bz,
                    double cx, double cy, double cz,
                    double dx, double dy, double dz,
                    double ex, double ey, double ez);
extern int powertest1d(double xa, double wa,
                       double xb, double wb,
                       double xc, double wc);
extern int powertest2d(double ax, double ay, double wa,
                       double bx, double by, double wb,
                       double cx, double cy, double wc,
                       double dx, double dy, double wd);
extern int powertest3d(double ax, double ay, double az, double wa,
                       double bx, double by, double bz, double wb,
                       double cx, double cy, double cz, double wc,
                       double dx, double dy, double dz, double wd,
                       double ex, double ey, double ez, double we);

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
s_segment_intersect segment_triangle_intersect_2D(const double S1[2], const double S2[2], const double A[2], const double B[2], const double C[2], double EPS_degenerate, double TOL_boundary); 
s_segment_intersect segment_triangle_intersect_3D(const s_point segment[2], const s_point triangle[3], double EPS_degenerate, double TOL_boundary);


#endif
