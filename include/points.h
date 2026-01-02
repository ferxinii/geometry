#ifndef GEOMETRY_POINTS_H
#define GEOMETRY_POINTS_H

#include <stdbool.h>
#include <stddef.h>

typedef struct point {
    union {
        double coords[3];
        struct {
            double x, y, z;
        };
    };
} s_point;

_Static_assert(sizeof(s_point) == 3 * sizeof(double), "size mismatch");

typedef struct points {
    int N;
    s_point *p;
} s_points;

#define point_NAN (s_point){{{ NAN, NAN, NAN }}}
int point_is_valid(s_point p);  /* Only way to check if point is NAN */
#define points_NAN (s_points){0, NULL}
int points_is_valid(const s_points *p);

void free_points(s_points *points);
s_points copy_points(const s_points *points);
void homotethy_points(s_points *points, double s, s_point pivot);
int mark_duplicate_points(const s_points *points, double TOL, bool mark[points->N]);
s_points copy_points_remove_duplicates(const s_points *points, double tol_d);
s_points read_points_from_csv(const char *file);
int write_points_to_csv(const char *file, const char *f_access_mode, const s_points *points);  /* f_access_mode: "w" or "a" */
void print_points(const s_points *points);

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
s_point scale_point(s_point a, double s);
s_point normalize_vec(s_point v, double EPS);  /* Returns point_NAN if norm < EPS */
s_point interpolate_points(s_point a, s_point b, double t);
void bounding_box_points(const s_points *points, s_point *min_out, s_point *max_out);
s_point span_points(const s_points *points);
double absolute_tolerance_from_scale(const s_points *points, double rel_tol);
s_point point_average(const s_points *points);
int coord_with_largest_component_3D(s_point n);
int coord_with_smallest_component_3D(s_point v);
s_point random_point_uniform_3D(s_point min, s_point max);
double area_triangle(const s_point face[3]);
double signed_volume_tetra(const s_point tetra[4]);
int basis_vectors_plane(const s_point plane[3], double EPS_degenerate, s_point *out_n, s_point *out_t1, s_point *out_t2);
int plane_from_point_normal(const s_point p0, const s_point n, double EPS_degenerate, s_point out[3]);
int plane_equation_from_points(const s_point plane[3], double EPS_degenerate, s_point *abc_out, double *d_out);
int circumcentre_tetrahedron(const s_point p[4], double EPS_degenerate, s_point *out);

/* The following return point_NAN if object is EPS-degenerate */
s_point project_point_to_plane(s_point p, const s_point plane[3], double EPS_degenerate);  
double signed_distance_point_to_plane(s_point p, const s_point plane[3], double EPS_degenerate);
s_point closest_point_on_segment(const s_point segment[2], double EPS_degenerate, s_point p); 
s_point closest_point_on_triangle(const s_point triangle[3], double EPS_degenerate, s_point p);  


#endif

