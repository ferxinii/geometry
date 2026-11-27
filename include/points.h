#ifndef _GEOMETRY_POINTS_H
#define _GEOMETRY_POINTS_H

#include <stdbool.h>
#include <stddef.h>
#include <math.h>

typedef struct point {
    union {
        double coords[3];
        struct {
            double x, y, z;
        };
    };
} s_point;

typedef struct point2d {
    union {
        double coords[2];
        struct {
            double x, y;
        };
    };
} s_point2d;

_Static_assert(sizeof(s_point) == 3 * sizeof(double), "size mismatch");

typedef struct points {
    int N;
    s_point *p;
} s_points;


#define point_NAN (s_point){{{ NAN, NAN, NAN }}}
bool point_is_valid(s_point p);  /* Only way to check if point is NAN */
/* Degeneracy is decided RELATIVE to the object's own size (scale invariant):
 * GEOM_EPS_REL is dimensionless and never user-set.  The checks compare squared
 * quantities (no sqrt). */
#define GEOM_EPS_REL 1e-12

/* |b-a| is zero at float resolution: |b-a|^2 <= eps^2 * max(|a|^2, |b|^2). */
static inline bool vec_is_degenerate(s_point a, s_point b)
{
    double dx = b.x - a.x, dy = b.y - a.y, dz = b.z - a.z;
    double na = a.x*a.x + a.y*a.y + a.z*a.z, nb = b.x*b.x + b.y*b.y + b.z*b.z;
    return dx*dx + dy*dy + dz*dz <= GEOM_EPS_REL * GEOM_EPS_REL * (na > nb ? na : nb);
}

/* Zero area relative to the edges: |ab x ac|^2 <= eps^2 * max(|ab|^2, |ac|^2)^2. */
static inline bool tri_is_degenerate(s_point a, s_point b, s_point c)
{
    double ux = b.x - a.x, uy = b.y - a.y, uz = b.z - a.z;
    double vx = c.x - a.x, vy = c.y - a.y, vz = c.z - a.z;
    double nx = uy*vz - uz*vy, ny = uz*vx - ux*vz, nz = ux*vy - uy*vx;
    double lu = ux*ux + uy*uy + uz*uz, lv = vx*vx + vy*vy + vz*vz;
    double M = lu > lv ? lu : lv;
    return nx*nx + ny*ny + nz*nz <= GEOM_EPS_REL * GEOM_EPS_REL * M * M;
}

/* Zero volume relative to the edges from p[0]: det^2 <= eps^2 * M^3, M = max |edge|^2. */
static inline bool tet_is_degenerate(const s_point p[4])
{
    double a[3], b[3], c[3];
    for (int k = 0; k < 3; k++) {
        a[k] = p[1].coords[k] - p[0].coords[k];
        b[k] = p[2].coords[k] - p[0].coords[k];
        c[k] = p[3].coords[k] - p[0].coords[k];
    }
    double det = a[0]*(b[1]*c[2] - b[2]*c[1]) - a[1]*(b[0]*c[2] - b[2]*c[0])
               + a[2]*(b[0]*c[1] - b[1]*c[0]);
    double la = a[0]*a[0] + a[1]*a[1] + a[2]*a[2];
    double lb = b[0]*b[0] + b[1]*b[1] + b[2]*b[2];
    double lc = c[0]*c[0] + c[1]*c[1] + c[2]*c[2];
    double M = la > lb ? la : lb;  M = M > lc ? M : lc;
    return det*det <= GEOM_EPS_REL * GEOM_EPS_REL * M * M * M;
}

#define points_NAN (s_points){0, NULL}
bool points_is_valid(const s_points *p);

void free_points(s_points *points);
s_points copy_points(const s_points *points);
void homotethy_points(s_points *points, double s, s_point pivot);
int mark_duplicate_points(const s_points *points, double tol_dup, bool mark[points->N]);
s_points copy_points_remove_duplicates(const s_points *points, double tol_dup);
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
s_point normalize_vec(s_point v);  /* point_NAN iff |v| is 0 or not finite */
s_point interpolate_points(s_point a, s_point b, double t);
void bounding_box_points(const s_points *points, s_point *min_out, s_point *max_out);
/* Fill perm[0..N-1] (caller-allocated) with a spatial Morton/Z-order ordering:
 * points->p[perm[0]], p[perm[1]], ... visits the cloud along a Z-order curve,
 * so consecutive entries are spatially close.  Intended as an insertion order
 * for incremental algorithms with locality (walk-hinted Delaunay: near-O(1)
 * point location per insert).  Deterministic: keys quantize to 21 bits/axis on
 * a uniform (isotropic) scale over the bbox; exact key ties break by original
 * index, so equal inputs always yield the same permutation.  Returns 1, or 0
 * on allocation failure (perm untouched). */
int points_morton_permutation(const s_points *points, int *perm);
s_point span_points(const s_points *points);
double absolute_tolerance_from_scale(const s_points *points, double rel_tol);
s_point point_average(const s_points *points);
int coord_with_largest_component_3D(s_point n);
int coord_with_smallest_component_3D(s_point v);
s_point random_point_uniform_3D(double (*randd01)(void *rctx), void *rctx,
                                s_point min, s_point max);
double area_triangle(const s_point face[3]);
double signed_volume_tetra(const s_point tetra[4]);
int basis_vectors_plane(const s_point plane[3], s_point *out_n, s_point *out_t1, s_point *out_t2);
int plane_from_point_normal(const s_point p0, const s_point n, s_point out[3]);
int plane_equation_from_points(const s_point plane[3], s_point *abc_out, double *d_out);
int circumcentre_tetrahedron(const s_point p[4], s_point *out);
/* Power circumcentre (equal power distance |x-p_i|^2 - w_i to all four). */
int circumcentre_tetrahedron_w(const s_point p[4], const double w[4] /* NULL = unweighted */,
                               s_point *out);

/* The following return point_NAN / NAN if the object is degenerate (see above) */
s_point project_point_to_plane(s_point p, const s_point plane[3]);  
void project_point_to_plane_2D(s_point p, s_point p0, s_point t1, s_point t2, double out[2]);
double signed_distance_point_to_plane(s_point p, const s_point plane[3]);
double signed_distance_point_to_plane_v2(s_point p, s_point plane_n, double plane_d);
s_point closest_point_on_segment(const s_point segment[2], s_point p);
double distance_sq_point_line(const s_point line[2], s_point p);
/* out_pa / out_pb: closest points on each segment; either may be NULL. */
double distance_segment_segment(const s_point seg_a[2], const s_point seg_b[2],
                                s_point *out_pa, s_point *out_pb);
s_point closest_point_on_triangle(const s_point triangle[3], s_point p);  


#endif

