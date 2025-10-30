#ifndef GEOMETRY_H
#define GEOMETRY_H

typedef struct point {
    union {
        double coords[3];
        struct {
            double x, y, z;
        };
    };
} s_point;

typedef struct points {
    int N;
    s_point *p;
} s_points;


// Basic s_points
s_points copy_points(const s_points *points);
void free_points(s_points *points);
s_points remove_duplicate_points(const s_points *points, double tol_dist);
s_points read_points_from_csv(const char *file);
int write_points_to_csv(const char *file, const char *f_access_mode, const s_points *points);  // f_access_mode: "w" or "a"
void print_points(const s_points *points);


// Basic geometrical operations
s_point sum_points(s_point p1, s_point p2);
s_point subtract_points(s_point u, s_point v);
s_point spherical_to_cartesian(double r, double theta, double phi);
double dot_prod(s_point u, s_point v);
s_point cross_prod(s_point u, s_point v);
double norm_squared(s_point v);
double norm(s_point v);
double distance_squared(s_point a, s_point b);
double distance(s_point a, s_point b);
double max_distance(const s_points *points, s_point query);
s_point normalize_3d(s_point v);

s_point point_average(const s_points *points);
int coord_with_largest_component_3d(s_point n);
s_point random_point_uniform_3d(s_point min, s_point max);


// Tests and closest on:
int segment_crosses_triangle_3d(const s_point triangle[3], s_point a, s_point b);
int between_1d(double x, double a, double b, double eps);
int segments_intersect_2d(const s_point AB[2], const s_point pd[2]);
s_point closest_point_on_triangle(const s_point triangle[3], s_point p);
s_point closest_point_on_segment(const s_point segment[2], s_point p);
int point_in_triangle_2d(const s_point triangle[3], s_point p);
int point_in_triangle_3d(const s_point triangle[3], s_point p);
int point_in_tetrahedron(const s_point tetra[4], s_point query);


// Predicates:
int orientation(const s_point p[3], s_point q);
int in_sphere(const s_point p[4], s_point q);
extern void exactinit(void);
extern double orient2d(const double *pa, const double *pb, const double *pc);
extern double orient3d(const double *pa, const double *pb, const double *pc, const double *pd);
extern double incircle(const double *pa, const double *pb, const double *pc, const double *pd);
extern double insphere(const double *pa, const double *pb, const double *pc, const double *pd, const double *pe);


double volume_tetrahedron_approx(s_point p1, s_point p2, s_point p3, s_point p4);

#endif
