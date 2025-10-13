#ifndef POINT_H
#define POINT_H

typedef struct point {
    union {
        double coords[3];
        struct {
            double x, y, z;
        };
    };
} s_point;


int orientation(const s_point *p3, s_point q);
int in_sphere(const s_point *p4, s_point q);

extern double orient2d(const double * pa, const double * pb, const double * pc);
extern double orient3d(const double * pa, const double * pb, const double * pc, const double * pd);
extern double incircle(const double * pa, const double * pb, const double * pc, const double * pd);
extern double insphere(const double * pa, const double * pb, const double * pc, const double * pd, const double * pe);

s_point sum_points(s_point p1, s_point p2);  // From BM
s_point subtract_points(s_point u, s_point v);
s_point spherical_to_cartesian(double r, double theta, double phi);

double dot_prod(s_point u, s_point v);
s_point cross_prod(s_point u, s_point v);

double norm_squared(s_point v);
double norm(s_point v);
double distance_squared(s_point a, s_point b);
double distance(s_point a, s_point b);
double max_distance(const s_point *p, int N, s_point query);

s_point normalize_3d(s_point v);

s_point find_center_mass(const s_point *in, int N_points);
int coord_with_largest_component_3d(s_point n);
s_point random_point_uniform_3d(s_point min, s_point max);

int segment_crosses_triangle_3d(const s_point *triangle, s_point a, s_point b);
int between_1d(double x, double a, double b, double eps);
int segments_intersect_2d(const s_point *AB, const s_point *pd);
s_point closest_point_on_triangle(const s_point *triangle, s_point p);
s_point closest_point_on_segment(const s_point *segment, s_point p);
int point_in_triangle_2d(const s_point *triangle, s_point p);
int point_in_triangle_3d(const s_point *triangle, s_point p);


#endif
