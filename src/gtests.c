/* tol_bdry: distance to a boundary that still counts as on it.
 * tol_bdry == 0: EXACT path, decided by exact predicates only; no epsilon test of
 *   any kind (callers such as the DT rely on this for near-flat faces).
 * tol_bdry > 0: degenerate inputs (relative to their own size, see points.h)
 *   give TEST_DEGENERATE / INTERSECT_ERROR.
 * Divisions that compute an output point after the predicates established the
 * configuration only guard against an exact zero / non-finite denominator. */

#include "points.h"
#include "gtests.h"
#include "predicates_basic.h"
#ifdef GEOMETRY_COMPILE_TESTS_POWER
#include "predicates_power.h"
#endif
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <stdatomic.h>
#include <float.h>

/* PREDICATES (wrappers of external/robust_predicates) */
int test_orientation_2d(const s_point2d line[2], s_point2d p)
{
    return orient2d(line[0].x, line[0].y, line[1].x, line[1].y, p.x, p.y);
}

int test_orientation(const s_point plane[3], s_point q)
{
    return orient3d(plane[0].x, plane[0].y, plane[0].z,
                    plane[1].x, plane[1].y, plane[1].z,
                    plane[2].x, plane[2].y, plane[2].z,
                    q.x, q.y, q.z);
}

int test_incircle(const s_point2d circle[3], s_point2d p)
{   
    return - powertest_n2_k3_unweighted(circle[0].x, circle[0].y,
                                        circle[1].x, circle[1].y,
                                        circle[2].x, circle[2].y,
                                        p.x, p.y);
}

int test_insphere(const s_point sph[4], s_point p)
{   
    return - powertest_n3_k4_unweighted(sph[0].x, sph[0].y, sph[0].z,
                                        sph[1].x, sph[1].y, sph[1].z,
                                        sph[2].x, sph[2].y, sph[2].z,
                                        sph[3].x, sph[3].y, sph[3].z,
                                        p.x, p.y, p.z);
}

bool test_points_concyclic(s_point a, s_point b, s_point c, s_point d)
{
    if (test_orientation((s_point[]){ a, b, c }, d) != 0) return false;   /* not coplanar */
    return incircle3d(a.x, a.y, a.z, b.x, b.y, b.z,
                      c.x, c.y, c.z, d.x, d.y, d.z) == 0;
}

bool test_points_colinear(s_point a, s_point b, s_point c)
{
    if (orient3d(a.x,a.y,a.z, b.x,b.y,b.z, c.x,c.y,c.z, a.x+1,a.y,a.z) == 0 &&
        orient3d(a.x,a.y,a.z, b.x,b.y,b.z, c.x,c.y,c.z, a.x,a.y+1,a.z) == 0 &&
        orient3d(a.x,a.y,a.z, b.x,b.y,b.z, c.x,c.y,c.z, a.x,a.y,a.z+1) == 0)
        return true;
    else return false;
}


#ifdef GEOMETRY_COMPILE_TESTS_POWER
int test_orthosegment(int k, const double c[k], const double wc[k], 
                      double xp, double wp)
{
    assert(1 <= k && k <= 2);
    
    switch (k) {
        case 1:
            return - powertest_n1_k1(c[0], wc[0], xp, wp); 
        case 2:
            return - powertest_n1_k2(c[0], wc[0], c[1], wc[1], xp, wp);
        default:
            printf("Unsupported case in test_orthosegment (%d not in {1,2})\n", k);
            exit(1);
    }
}

int test_orthocircle(int k, const s_point2d c[k], const double wc[k], 
                     s_point2d p, double wp)
{   
    assert(1 <= k && k <= 3);
    
    switch (k) {
        case 1:
            return - powertest_n2_k1(c[0].x, c[0].y, wc[0],
                                     p.x, p.y, wp);
        case 2:
            return - powertest_n2_k2(c[0].x, c[0].y, wc[0],
                                     c[1].x, c[1].y, wc[1],
                                     p.x, p.y, wp);
        case 3:
            return - powertest_n2_k3(c[0].x, c[0].y, wc[0],
                                     c[1].x, c[1].y, wc[1],
                                     c[2].x, c[2].y, wc[2],
                                     p.x, p.y, wp);
        default:
            printf("Unsupported case in test_orthocircle (%d not in {1,2,3})\n", k);
            exit(1);
    }
 }

int test_orthosphere(int k, const s_point c[k], const double wc[k],
                     s_point p, double wp)
{   
    assert(1 <= k && k <= 4);

    switch (k) {
        case 1:
            return - powertest_n3_k1(c[0].x, c[0].y, c[0].z, wc[0],
                                     p.x, p.y, p.z, wp);
        case 2:
            return - powertest_n3_k2(c[0].x, c[0].y, c[0].z, wc[0],
                                     c[1].x, c[1].y, c[1].z, wc[1],
                                     p.x, p.y, p.z, wp);
        case 3:
            return - powertest_n3_k3(c[0].x, c[0].y, c[0].z, wc[0],
                                     c[1].x, c[1].y, c[1].z, wc[1],
                                     c[2].x, c[2].y, c[2].z, wc[2],
                                     p.x, p.y, p.z, wp);
        case 4:
            return - powertest_n3_k4(c[0].x, c[0].y, c[0].z, wc[0],
                                     c[1].x, c[1].y, c[1].z, wc[1],
                                     c[2].x, c[2].y, c[2].z, wc[2],
                                     c[3].x, c[3].y, c[3].z, wc[3],
                                     p.x, p.y, p.z, wp);
        default:
            printf("Unsupported case in test_orthosphere (%d not in {1,2,3,4})\n", k);
            exit(1);
    }
}

int test_orthosegment_w(int k, const double c[k], const double wc[k], 
                        double alpha)
{
    assert(1 <= k && k <= 2);
    
    switch (k) {
        case 1:
            return orthow_n1_k1(c[0],wc[0], alpha); 
        case 2:
            return orthow_n1_k2(c[0],wc[0], c[1],wc[1], alpha);
        default:
            printf("Unsupported case in test_orthosegment_w (%d not in {1,2})\n", k);
            exit(1);
    }
}

int test_orthocircle_w(int k, const s_point2d c[k], const double wc[k], 
                       double alpha)
{   
    assert(1 <= k && k <= 3);
    
    switch (k) {
        case 1:
            return orthow_n2_k1(c[0].x, c[0].y, wc[0],
                                  alpha);
        case 2:
            return orthow_n2_k2(c[0].x, c[0].y, wc[0],
                                c[1].x, c[1].y, wc[1],
                                alpha);
        case 3:
            return orthow_n2_k3(c[0].x, c[0].y, wc[0],
                                c[1].x, c[1].y, wc[1],
                                c[2].x, c[2].y, wc[2],
                                alpha);
        default:
            printf("Unsupported case in test_orthocircle_w (%d not in {1,2,3})\n", k);
            exit(1);
    }
 }

int test_orthosphere_w(int k, const s_point c[k], const double wc[k],
                       double alpha)
{   
    assert(1 <= k && k <= 4);

    switch (k) {
        case 1:
            return orthow_n3_k1(c[0].x, c[0].y, c[0].z, wc[0],
                                alpha);
        case 2:
            return orthow_n3_k2(c[0].x, c[0].y, c[0].z, wc[0],
                                c[1].x, c[1].y, c[1].z, wc[1],
                                alpha);
        case 3:
            return orthow_n3_k3(c[0].x, c[0].y, c[0].z, wc[0],
                                c[1].x, c[1].y, c[1].z, wc[1],
                                c[2].x, c[2].y, c[2].z, wc[2],
                                alpha);
        case 4:
            return orthow_n3_k4(c[0].x, c[0].y, c[0].z, wc[0],
                                c[1].x, c[1].y, c[1].z, wc[1],
                                c[2].x, c[2].y, c[2].z, wc[2],
                                c[3].x, c[3].y, c[3].z, wc[3],
                                alpha);
        default:
            printf("Unsupported case in test_orthosphere_w (%d not in {1,2,3,4})\n", k);
            exit(1);
    }
}

bool test_points_coorthocircular(s_point a, double wa, s_point b, double wb,
                                 s_point c, double wc, s_point d, double wd)
{   /* Coplanarity gate is unweighted (purely geometric); power_incircle3d == 0
     * iff d is on the orthogonal circle of (a,wa),(b,wb),(c,wc). */
    if (test_orientation((s_point[]){ a, b, c }, d) != 0) return false;   /* not coplanar */
    return power_incircle3d(a.x, a.y, a.z, wa, b.x, b.y, b.z, wb,
                            c.x, c.y, c.z, wc, d.x, d.y, d.z, wd) == 0;
}
#endif  // GEOMETRY_COMPILE_TESTS_POWER


/* HELPERS */
static int coord_to_drop_from_plane(const s_point plane[3])
{
    s_point AB = subtract_points(plane[1], plane[0]);
    s_point AC = subtract_points(plane[2], plane[0]);
    s_point n  = cross_prod(AB, AC);
    return coord_with_largest_component_3D(n);
}

static inline bool seg2_is_degenerate(s_point2d a, s_point2d b)
{   /* 2D vec_is_degenerate */
    double dx = b.x - a.x, dy = b.y - a.y;
    double na = a.x*a.x + a.y*a.y, nb = b.x*b.x + b.y*b.y;
    return dx*dx + dy*dy <= GEOM_EPS_REL * GEOM_EPS_REL * fmax(na, nb);
}

static inline bool tri2_is_degenerate(const s_point2d t[3])
{   /* 2D tri_is_degenerate */
    double ux = t[1].x - t[0].x, uy = t[1].y - t[0].y;
    double vx = t[2].x - t[0].x, vy = t[2].y - t[0].y;
    double c = ux*vy - uy*vx, M = fmax(ux*ux + uy*uy, vx*vx + vy*vy);
    return c*c <= GEOM_EPS_REL * GEOM_EPS_REL * M * M;
}

static s_point2d drop_to_2D(const s_point p, int coord_to_drop)
{
    int i1 = (coord_to_drop + 1) % 3;
    int i2 = (coord_to_drop + 2) % 3;
    s_point2d out;
    out.x = p.coords[i1];
    out.y = p.coords[i2];
    return out;
}

static s_point lift_point_from_dropped_2D(const s_point plane[3], int drop, s_point2d in)
{
    int i1 = (drop + 1) % 3;
    int i2 = (drop + 2) % 3;
    s_point out;
    out.coords[i1] = in.x;
    out.coords[i2] = in.y;
    s_point n = cross_prod(subtract_points(plane[1], plane[0]), subtract_points(plane[2], plane[0]));
    if (fabs(n.coords[drop]) > 0) {   /* drop = largest |n| component */
        double t0i = plane[0].coords[i1];
        double t0j = plane[0].coords[i2];
        double numer = - ( n.coords[i1]*(out.coords[i1]-t0i) + n.coords[i2]*(out.coords[i2]-t0j) );
        out.coords[drop] = plane[0].coords[drop] + numer / n.coords[drop];
        return out;
    }
    return point_NAN;
}

static double point_segment_dist2_2D(s_point2d p, const s_point2d ab[2])
{
    s_point2d a = ab[0], b = ab[1];

    double vx = b.x - a.x, vy = b.y - a.y;
    double wx = p.x - a.x, wy = p.y - a.y;

    double denom = vx*vx + vy*vy;
    if (seg2_is_degenerate(a, b)) {
        return (p.x-a.x)*(p.x-a.x) + (p.y-a.y)*(p.y-a.y);
    }
    double t = (wx*vx + wy*vy) / denom;

    if (t <= 0.0) return (p.x-a.x)*(p.x-a.x) + (p.y-a.y)*(p.y-a.y);
    else if (t >= 1.0) return (p.x-b.x)*(p.x-b.x) + (p.y-b.y)*(p.y-b.y);
    else {
        double projx = a.x + t*vx;
        double projy = a.y + t*vy;
        double dx = p.x - projx, dy = p.y - projy;
        return dx*dx + dy*dy;
    }
}

static int points_close_2D(s_point2d a, s_point2d b, double tol_bdry)
{
    if (tol_bdry == 0) {
        if (a.x == b.x && a.y == b.y) return 1;
        else return 0;
    } else {
        if ((a.x-b.x)*(a.x-b.x) + (a.y-b.y)*(a.y-b.y) <= tol_bdry*tol_bdry) return 1;
        else return 0;
    }
}

static int append_unique_2D_lim2(s_point2d out[2], int *count, s_point2d p, double tol_bdry)
{
    if ((*count)<2 && tol_bdry > 0) {
        for (int i=0; i<(*count); i++) {
            if (points_close_2D(out[i], p, tol_bdry)) return 0;
        }
    }
    if (out) {
        out[*count].x = p.x;
        out[*count].y = p.y;
    }
    (*count)++;
    return 1;
}



/* GEOMETRICAL TESTS */
e_geom_test test_point_in_interval_1D(double x, double a, double b, double tol_bdry)
{
    if (tol_bdry > 0 && fabs(b - a) <= GEOM_EPS_REL * fmax(fabs(a), fabs(b))) return TEST_DEGENERATE;
    if (fabs(x - a) <= tol_bdry|| fabs(x - b) <= tol_bdry) return TEST_BOUNDARY;
    if (a > b) { double tmp = a; a = b; b = tmp; }
    if (x + tol_bdry>= a && x - tol_bdry<= b) return TEST_IN;
    return TEST_OUT;
}


static e_geom_test test_point_in_triangle_2D_robust_from_orientations(int o1, int o2, int o3)
{
    /* If any non-zero orientation disagrees, p is outside */
    int ref = (o1 != 0) ? o1 : ((o2 != 0) ? o2 : o3);
    if (ref == 0) return TEST_DEGENERATE;
    if ((o1 != 0 && o1 != ref) || (o2 != 0 && o2 != ref) || (o3 != 0 && o3 != ref)) return TEST_OUT;
    if (o1 == 0 || o2 == 0 || o3 == 0) return TEST_BOUNDARY;
    return TEST_IN;
}


e_geom_test test_point_in_triangle_2D(const s_point2d triangle[3], s_point2d p, double tol_bdry)
{
    if (tol_bdry > 0 && tri2_is_degenerate(triangle)) return TEST_DEGENERATE;

    /* Robust orientation tests of p vs each edge */
    int o1 = test_orientation_2d((s_point2d[]){triangle[0], triangle[1]}, p);
    int o2 = test_orientation_2d((s_point2d[]){triangle[1], triangle[2]}, p);
    int o3 = test_orientation_2d((s_point2d[]){triangle[2], triangle[0]}, p);

    if (tol_bdry == 0) return test_point_in_triangle_2D_robust_from_orientations(o1, o2, o3);

    /* Numerical tolerance branch */
    double tol2 = tol_bdry*tol_bdry;
    if (point_segment_dist2_2D(p, (s_point2d[]){triangle[0], triangle[1]}) <= tol2) return TEST_BOUNDARY;
    if (point_segment_dist2_2D(p, (s_point2d[]){triangle[1], triangle[2]}) <= tol2) return TEST_BOUNDARY;
    if (point_segment_dist2_2D(p, (s_point2d[]){triangle[2], triangle[0]}) <= tol2) return TEST_BOUNDARY;
    
    /* If not on boundary, check robust orientation to determine TEST_IN/TEST_OUT */
    return test_point_in_triangle_2D_robust_from_orientations(o1, o2, o3);
}


static e_geom_test test_point_in_triangle_3D_robust(const s_point triangle[3], s_point p)
{
    if (test_orientation(triangle, p) != 0) return TEST_OUT;

    /* Point coplanar. Project to 2D and check */   
    int drop = coord_to_drop_from_plane(triangle);
    s_point2d A2 = drop_to_2D(triangle[0], drop);
    s_point2d B2 = drop_to_2D(triangle[1], drop);
    s_point2d C2 = drop_to_2D(triangle[2], drop);
    s_point2d p2 = drop_to_2D(p, drop);

    return test_point_in_triangle_2D((s_point2d[]){A2, B2, C2}, p2, 0);
}


e_geom_test test_point_in_triangle_3D(const s_point triangle[3], s_point p, double tol_bdry) 
{
    if (tol_bdry > 0 && tri_is_degenerate(triangle[0], triangle[1], triangle[2])) return TEST_DEGENERATE;

    if (tol_bdry == 0) return test_point_in_triangle_3D_robust(triangle, p);

    s_point closest = project_point_to_plane(p, triangle);
    if (!point_is_valid(closest)) return TEST_ERROR;
    if (distance_squared(closest, p) > tol_bdry*tol_bdry) return TEST_OUT;
    
    /* Point is tol_bdry close to plane. Project to 2D and check */
    int drop = coord_to_drop_from_plane(triangle);
    s_point2d A2 = drop_to_2D(triangle[0], drop);
    s_point2d B2 = drop_to_2D(triangle[1], drop);
    s_point2d C2 = drop_to_2D(triangle[2], drop);
    s_point2d p2 = drop_to_2D(closest, drop);
    return test_point_in_triangle_2D((s_point2d[]){A2, B2, C2}, p2, tol_bdry);
}


static e_geom_test test_point_in_tetrahedron_robust(const s_point tetra[4], s_point query)
{   /* No need to be properly oriented */
    s_point tmp[3];

    /* Reference signs: e, query signs: s */
    tmp[0] = tetra[1];   tmp[1] = tetra[2];   tmp[2] = tetra[3];
    int e0 = test_orientation(tmp, tetra[0]);
    int s0 = test_orientation(tmp, query);

    tmp[0] = tetra[0];   tmp[1] = tetra[3];   tmp[2] = tetra[2];
    int e1 = test_orientation(tmp, tetra[1]);
    int s1 = test_orientation(tmp, query);

    tmp[0] = tetra[0];   tmp[1] = tetra[1];   tmp[2] = tetra[3];
    int e2 = test_orientation(tmp, tetra[2]);
    int s2 = test_orientation(tmp, query);

    tmp[0] = tetra[0];   tmp[1] = tetra[2];   tmp[2] = tetra[1];
    int e3 = test_orientation(tmp, tetra[3]);
    int s3 = test_orientation(tmp, query);

    if (e0 == 0 || e1 == 0 || e2 == 0 || e3 == 0) return TEST_DEGENERATE;

    if ((s0 != 0 && s0 * e0 < 0) ||
        (s1 != 0 && s1 * e1 < 0) ||
        (s2 != 0 && s2 * e2 < 0) ||
        (s3 != 0 && s3 * e3 < 0))
        return TEST_OUT;
    

    if (s0 == 0 || s1 == 0 || s2 == 0 || s3 == 0) return TEST_BOUNDARY;
    return TEST_IN;
}


e_geom_test test_point_in_tetrahedron(const s_point tetra[4], s_point query, double tol_bdry)
{
    if (tol_bdry > 0 && tet_is_degenerate(tetra)) return TEST_DEGENERATE;

    if (tol_bdry == 0) return test_point_in_tetrahedron_robust(tetra, query);

    /* First check if query is tol-close to each face */
    s_point tmp[3];
    tmp[0] = tetra[1];   tmp[1] = tetra[2];   tmp[2] = tetra[3];
    e_geom_test t1 = test_point_in_triangle_3D(tmp, query, tol_bdry);
    
    tmp[0] = tetra[0];   tmp[1] = tetra[3];   tmp[2] = tetra[2];
    e_geom_test t2 = test_point_in_triangle_3D(tmp, query, tol_bdry);

    tmp[0] = tetra[0];   tmp[1] = tetra[1];   tmp[2] = tetra[3];
    e_geom_test t3 = test_point_in_triangle_3D(tmp, query, tol_bdry);

    tmp[0] = tetra[0];   tmp[1] = tetra[2];   tmp[2] = tetra[1];
    e_geom_test t4 = test_point_in_triangle_3D(tmp, query, tol_bdry);

    if (t1 == TEST_ERROR || t2 == TEST_ERROR || t3 == TEST_ERROR || t4 == TEST_ERROR) return TEST_ERROR;
    if (t1 == TEST_DEGENERATE || t2 == TEST_DEGENERATE || t3 == TEST_DEGENERATE || t4 == TEST_DEGENERATE) return TEST_DEGENERATE;
    else if (t1 == TEST_IN || t1 == TEST_BOUNDARY ||
             t2 == TEST_IN || t2 == TEST_BOUNDARY ||
             t3 == TEST_IN || t3 == TEST_BOUNDARY ||
             t4 == TEST_IN || t4 == TEST_BOUNDARY) 
        return TEST_BOUNDARY;
    
    /* Do robust test to determine if TEST_IN / TEST_OUT */
    return test_point_in_tetrahedron_robust(tetra, query);
}


s_points_test test_points_in_halfspace(const s_point plane_ordered[3], const s_points *points, double tol_bdry,
                                       e_geom_test out_buff[points->N])
{   /* Right hand rule: normal outwards */
    if (out_buff == NULL) out_buff = malloc(points->N * sizeof(e_geom_test));
    int Nin = 0, Nbdy = 0, Nout = 0;

    if (tol_bdry == 0) {  /* Robust branch */
        for (int ii=0; ii<points->N; ii++) {
            int o = test_orientation(plane_ordered, points->p[ii]);
            if (o == 1) { out_buff[ii] = TEST_IN; Nin++; } 
            else if (o == 0) { out_buff[ii] = TEST_BOUNDARY; Nbdy++; }
            else { out_buff[ii] = TEST_OUT; Nout++; }
        } 
    } else {  /* Non-robust branch */
        s_point n = cross_prod(subtract_points(plane_ordered[1], plane_ordered[0]),
                               subtract_points(plane_ordered[2], plane_ordered[0]));
        n = tri_is_degenerate(plane_ordered[0], plane_ordered[1], plane_ordered[2])
            ? point_NAN : normalize_vec(n);
        if (!point_is_valid(n)) {
            for (int ii=0; ii<points->N; ii++) out_buff[ii] = TEST_ERROR;
            return (s_points_test){ .Nin = 0, .Nbdy = 0, .Nout = 0, .Nerr = points->N, .indicator = out_buff };
        }
        double d = dot_prod(n, plane_ordered[0]);

        for (int ii=0; ii<points->N; ii++) {
            double s = dot_prod(points->p[ii], n) - d;
            if (s < -tol_bdry) { out_buff[ii] = TEST_IN; Nin++; } 
            else if (fabs(s) <= tol_bdry) { out_buff[ii] = TEST_BOUNDARY; Nbdy++; } 
            else { out_buff[ii] = TEST_OUT; Nout++; }
        }
    }
    
    return (s_points_test){ .Nin = Nin, .Nbdy = Nbdy, .Nout = Nout, .Nerr = 0, .indicator = out_buff };
}




/* INTERSECTIONS */
/* Segment segment intersection 2D */
static void compute_segment_segment_intersection_2D_nondegenerate(s_point2d A1, s_point2d A2,
                                                                   s_point2d B1, s_point2d B2,
                                                                   int *Nout, s_point2d *out)
{
    double A_a = A2.y - A1.y;
    double A_b = A1.x - A2.x;
    double A_c = A_a*A1.x + A_b*A1.y;
    double B_a = B2.y - B1.y;
    double B_b = B1.x - B2.x;
    double B_c = B_a*B1.x + B_b*B1.y;

    double det = A_a*B_b - B_a*A_b;
    if (det == 0 || !isfinite(det)) { *Nout = 0; }
    else {
        out[0].x = (B_b*A_c - A_b*B_c) / det;
        out[0].y = (A_a*B_c - B_a*A_c) / det;
        *Nout = 1;
    }
}

static void compute_segment_segment_intersection_2D_colinear(s_point2d A1, s_point2d A2,
                                                              s_point2d B1, s_point2d B2,
                                                              double tol_bdry,
                                                              int *Nout, s_point2d *out)
{
    double dx = A2.x - A1.x, dy = A2.y - A1.y;
    double denom = dx*dx + dy*dy;
    if (!(denom > 0) || !isfinite(denom)) { *Nout = 0; return; }
    double tolt = tol_bdry / sqrt(denom);   /* tol_bdry in units of the parameter t */

    double tB1 = ((B1.x-A1.x)*dx + (B1.y-A1.y)*dy) / denom;
    double tB2 = ((B2.x-A1.x)*dx + (B2.y-A1.y)*dy) / denom;

    double tmin = fmin(tB1, tB2), tmax = fmax(tB1, tB2);
    double lo = fmax(0.0, tmin), hi = fmin(1.0, tmax);

    if (hi + tolt < lo - tolt) { *Nout = 0; return; }
    if (fabs(hi - lo) <= tolt) {
        double t = 0.5 * (hi + lo);
        out[0].x = A1.x + t * dx;
        out[0].y = A1.y + t * dy;
        *Nout = 1;
        return;
    }
    if (lo < 0.0) lo = 0.0;
    if (hi > 1.0) hi = 1.0;
    out[0].x = A1.x + lo * dx;
    out[0].y = A1.y + lo * dy;
    out[1].x = A1.x + hi * dx;
    out[1].y = A1.y + hi * dy;
    *Nout = 2;
}

static e_intersect_type core_segment_segment_intersect_2D(const s_point2d s1[2], const s_point2d s2[2],
                                                          double tol_bdry,
                                                          int *Nout, s_point2d *out)
{
    if (Nout) *Nout = 0;

    s_point2d A1 = s1[0], A2 = s1[1], B1 = s2[0], B2 = s2[1];

    double tol2 = tol_bdry*tol_bdry;
    int A_is_point = (tol_bdry == 0) ? points_close_2D(A1, A2, 0) : seg2_is_degenerate(A1, A2);
    int B_is_point = (tol_bdry == 0) ? points_close_2D(B1, B2, 0) : seg2_is_degenerate(B1, B2);
    if (A_is_point && B_is_point) {
        if (points_close_2D(A1, B1, tol_bdry)) {
            if (Nout && out) {
                out[0].x = (A1.x + A2.x + B1.x + B2.x) / 4.0;
                out[0].y = (A1.y + A2.y + B1.y + B2.y) / 4.0;
                *Nout = 1;
            }
            return INTERSECT_DEGENERATE;
        }
        if (Nout) *Nout = 0;
        return INTERSECT_EMPTY;
    }
    if (A_is_point) {
        if ((tol_bdry==0) ? (test_orientation_2d((s_point2d[]){B1, B2}, A1) == 0 &&
                        test_point_in_interval_1D(A1.x, B1.x, B2.x, tol_bdry) != TEST_OUT &&
                        test_point_in_interval_1D(A1.y, B1.y, B2.y, tol_bdry) != TEST_OUT)
                     : (point_segment_dist2_2D(A1, s2) <= tol2)) {
            if (Nout && out) { out[0] = A1; *Nout = 1; }
            return INTERSECT_DEGENERATE;
        }
        if (Nout) *Nout = 0;
        return INTERSECT_EMPTY;
    }
    if (B_is_point) {
        if ((tol_bdry==0) ? (test_orientation_2d((s_point2d[]){A1, A2}, B1) == 0 &&
                        test_point_in_interval_1D(B1.x, A1.x, A2.x, tol_bdry) != TEST_OUT &&
                        test_point_in_interval_1D(B1.y, A1.y, A2.y, tol_bdry) != TEST_OUT)
                     : (point_segment_dist2_2D(B1, s1) <= tol2)) {
            if (Nout && out) { out[0] = B1; *Nout = 1; }
            return INTERSECT_DEGENERATE;
        }
        if (Nout) *Nout = 0;
        return INTERSECT_EMPTY;
    }

    double d1 = (tol_bdry==0) ? test_orientation_2d((s_point2d[]){A1, A2}, B1) : point_segment_dist2_2D(B1, s1);
    double d2 = (tol_bdry==0) ? test_orientation_2d((s_point2d[]){A1, A2}, B2) : point_segment_dist2_2D(B2, s1);
    double d3 = (tol_bdry==0) ? test_orientation_2d((s_point2d[]){B1, B2}, A1) : point_segment_dist2_2D(A1, s2);
    double d4 = (tol_bdry==0) ? test_orientation_2d((s_point2d[]){B1, B2}, A2) : point_segment_dist2_2D(A2, s2);

    if ((tol_bdry==0) ? (d1 == 0 && d2 == 0 && d3 == 0 && d4 == 0) :
                   (d1 <= tol2 && d2 <= tol2 && d3 <= tol2 && d4 <= tol2)) {
        double min_x = fmax(fmin(A1.x, A2.x), fmin(B1.x, B2.x));
        double max_x = fmin(fmax(A1.x, A2.x), fmax(B1.x, B2.x));
        double min_y = fmax(fmin(A1.y, A2.y), fmin(B1.y, B2.y));
        double max_y = fmin(fmax(A1.y, A2.y), fmax(B1.y, B2.y));
        if (min_x <= max_x + tol_bdry && min_y <= max_y + tol_bdry) {
            if (Nout && out) compute_segment_segment_intersection_2D_colinear(A1, A2, B1, B2, tol_bdry, Nout, out);
            return INTERSECT_DEGENERATE;
        } else {
            if (Nout) *Nout = 0;
            return INTERSECT_EMPTY;
        }
    }

    if ((tol_bdry==0) ? (d1==0 &&
                    test_point_in_interval_1D(B1.x, A1.x, A2.x, tol_bdry) != TEST_OUT &&
                    test_point_in_interval_1D(B1.y, A1.y, A2.y, tol_bdry) != TEST_OUT)
                 : (d1<=tol2)) {
        if (Nout && out) { out[0] = B1; *Nout = 1; }
        return INTERSECT_DEGENERATE;
    }
    if ((tol_bdry==0) ? (d2==0 &&
                    test_point_in_interval_1D(B2.x, A1.x, A2.x, tol_bdry) != TEST_OUT &&
                    test_point_in_interval_1D(B2.y, A1.y, A2.y, tol_bdry) != TEST_OUT)
                 : (d2<=tol2)) {
        if (Nout && out) { out[0] = B2; *Nout = 1; }
        return INTERSECT_DEGENERATE;
    }
    if ((tol_bdry==0) ? (d3==0 &&
                    test_point_in_interval_1D(A1.x, B1.x, B2.x, tol_bdry) != TEST_OUT &&
                    test_point_in_interval_1D(A1.y, B1.y, B2.y, tol_bdry) != TEST_OUT)
                 : (d3<=tol2)) {
        if (Nout && out) { out[0] = A1; *Nout = 1; }
        return INTERSECT_DEGENERATE;
    }
    if ((tol_bdry==0) ? (d4==0 &&
                    test_point_in_interval_1D(A2.x, B1.x, B2.x, tol_bdry) != TEST_OUT &&
                    test_point_in_interval_1D(A2.y, B1.y, B2.y, tol_bdry) != TEST_OUT)
                 : (d4<=tol2)) {
        if (Nout && out) { out[0] = A2; *Nout = 1; }
        return INTERSECT_DEGENERATE;
    }

    if (tol_bdry != 0) {
        d1 = test_orientation_2d((s_point2d[]){A1, A2}, B1);
        d2 = test_orientation_2d((s_point2d[]){A1, A2}, B2);
        d3 = test_orientation_2d((s_point2d[]){B1, B2}, A1);
        d4 = test_orientation_2d((s_point2d[]){B1, B2}, A2);
    }
    if ( ((d1>0 && d2<0) || (d1<0 && d2>0)) &&
         ((d3>0 && d4<0) || (d3<0 && d4>0)) ) {
        if (Nout && out) compute_segment_segment_intersection_2D_nondegenerate(A1, A2, B1, B2, Nout, out);
        return INTERSECT_NONDEGENERATE;
    }

    if (Nout) *Nout = 0;
    return INTERSECT_EMPTY;
}

e_intersect_type test_segment_segment_intersect_2D(const s_point2d s1[2], const s_point2d s2[2], double tol_bdry)
{
    return core_segment_segment_intersect_2D(s1, s2, tol_bdry, NULL, NULL);
}

s_segment_intersect segment_segment_intersect_2D(const s_point2d s1[2], const s_point2d s2[2], double tol_bdry)
{
    s_segment_intersect out;
    s_point2d coords[2];
    out.type = core_segment_segment_intersect_2D(s1, s2, tol_bdry, &out.N, coords);
    for (int ii=0; ii<out.N; ii++) {
        out.coords[ii].x = coords[ii].x;
        out.coords[ii].y = coords[ii].y;
        out.coords[ii].z = 0;
    }
    return out;
}





/* segment plane intersection */
static e_intersect_type core_segment_plane_intersect_robust(const s_point seg[2], const s_point plane[3],
                                                            int *Nout, s_point out[2])
{
    int o0 = test_orientation(plane, seg[0]);
    int o1 = test_orientation(plane, seg[1]);

    if (o0 == 0 && o1 == 0) { 
        if (Nout && out) { out[0] = seg[0]; out[1] = seg[1]; *Nout = 2; } 
        return INTERSECT_DEGENERATE;
    }
    if (o0 == 0) { 
        if (Nout && out) { out[0] = seg[0]; *Nout = 1; }
        return INTERSECT_DEGENERATE;
    }
    if (o1 == 0) { 
        if (Nout && out) { out[0] = seg[1]; *Nout = 1; } 
        return INTERSECT_DEGENERATE;
    }
    if (o0 == o1) {  /* same side of plane */
        if (Nout) *Nout = 0; 
        return INTERSECT_EMPTY; 
    }  
    
    if (Nout && out) {
        s_point n = cross_prod(subtract_points(plane[1], plane[0]),
                               subtract_points(plane[2], plane[0]));
        double d = dot_prod(n, plane[0]);
        double sA = dot_prod(seg[0], n) - d;
        double sB = dot_prod(seg[1], n) - d;
        if (sA == sB || !isfinite(sA - sB)) { *Nout = 0; return INTERSECT_NONDEGENERATE; }
        double t = sA / (sA - sB);
        out[0] = interpolate_points(seg[0], seg[1], t); 
        *Nout = 1;
    }
    return INTERSECT_NONDEGENERATE;
}

static e_intersect_type core_segment_plane_intersect(const s_point seg[2], const s_point plane[3],
                                                     double tol_bdry,
                                                     int *Nout, s_point out[2])
{
    if (Nout) *Nout = 0;
    if (tol_bdry == 0) return core_segment_plane_intersect_robust(seg, plane, Nout, out);

    if (tri_is_degenerate(plane[0], plane[1], plane[2])) return INTERSECT_ERROR;
    s_point n = normalize_vec(cross_prod(subtract_points(plane[1], plane[0]),
                                         subtract_points(plane[2], plane[0])));
    double d = dot_prod(n, plane[0]);   /* sA, sB below: signed distances */
    double sA = dot_prod(seg[0], n) - d;
    double sB = dot_prod(seg[1], n) - d;

    if (sA * sB > tol_bdry*tol_bdry) {  /* both endpoints on same side (and not near plane) */
        if (Nout && out) *Nout = 0;
        return INTERSECT_EMPTY;
    }

    if (fabs(sA) < tol_bdry && fabs(sB) < tol_bdry) {  /* whole segment lies in plane */
        if (Nout && out) {
            out[0] = project_point_to_plane(seg[0], plane);
            out[1] = project_point_to_plane(seg[1], plane);
            *Nout = point_is_valid(out[0]) + point_is_valid(out[1]); 
        }
        return INTERSECT_DEGENERATE;
    }
    if (fabs(sA) < tol_bdry || fabs(sB) < tol_bdry) {  /* Only one end close to plane */
        if (Nout && out) {
            out[0] = project_point_to_plane((fabs(sA)<tol_bdry ? seg[0] : seg[1]), plane);
            *Nout = point_is_valid(out[0]);
        }
        return INTERSECT_DEGENERATE;
    }

    if (Nout && out) {  /* find parametric intersection */
        double denom = sA - sB;
        if (denom == 0 || !isfinite(denom)) { *Nout = 0; return INTERSECT_NONDEGENERATE; }
        double t = sA / denom;
        double tolt = tol_bdry / distance(seg[0], seg[1]);   /* tol_bdry in units of t */
        if (t < -tolt || t > 1 + tolt) { *Nout = 0; return INTERSECT_NONDEGENERATE; }
        if (t < 0) t = 0;
        if (t > 1) t = 1;
        out[0] = interpolate_points(seg[0], seg[1], t);
        *Nout = 1;
    }
    return INTERSECT_NONDEGENERATE;
}

e_intersect_type test_segment_plane_intersect(const s_point seg[2], const s_point plane[3], double tol_bdry)
{
    return core_segment_plane_intersect(seg, plane, tol_bdry, NULL, NULL);
}

s_segment_intersect segment_plane_intersect(const s_point seg[2], const s_point plane[3], double tol_bdry)
{
    s_segment_intersect out;
    out.type = core_segment_plane_intersect(seg, plane, tol_bdry, &out.N, out.coords);
    return out;
}



/* Segment triangle intersection */
static e_intersect_type core_segment_triangle_intersect_2D(const s_point2d s[2],
                                                           const s_point2d tri[3],
                                                           double tol_bdry,
                                                           int *Nout, s_point2d *out)
{
    if (Nout) *Nout = 0;
    s_point2d S1 = s[0], S2 = s[1];

    e_geom_test i1 = test_point_in_triangle_2D(tri, S1, tol_bdry);
    e_geom_test i2 = test_point_in_triangle_2D(tri, S2, tol_bdry);
    if (i1 == TEST_ERROR || i2 == TEST_ERROR || i1 == TEST_DEGENERATE || i2 == TEST_DEGENERATE) return INTERSECT_ERROR;

    if (i1 == TEST_BOUNDARY && i2 == TEST_BOUNDARY) {
        if (Nout && out) { out[0] = S1; out[1] = S2; *Nout = 2; }
        return INTERSECT_DEGENERATE;
    }
    if ((i1 == TEST_IN || i1 == TEST_BOUNDARY) && (i2 == TEST_IN || i2 == TEST_BOUNDARY)) {
        if (Nout && out) { out[0] = S1; out[1] = S2; *Nout = 2; }
        return INTERSECT_NONDEGENERATE;
    }

    s_segment_intersect t1 = segment_segment_intersect_2D(s, (s_point2d[]){tri[0], tri[1]}, tol_bdry);
    s_segment_intersect t2 = segment_segment_intersect_2D(s, (s_point2d[]){tri[1], tri[2]}, tol_bdry);
    s_segment_intersect t3 = segment_segment_intersect_2D(s, (s_point2d[]){tri[2], tri[0]}, tol_bdry);

    if ((i1 == TEST_IN || i2 == TEST_IN) && (t1.N > 0 || t2.N > 0 || t3.N > 0)) {
        if (Nout && out) {
            out[0] = (i1 == TEST_IN) ? S1 : S2;
            s_point *coords = t1.N > 0 ? t1.coords : (t2.N > 0 ? t2.coords : t3.coords);
            out[1].x = coords[0].x; out[1].y = coords[0].y;
            *Nout = 2;
        }
        return INTERSECT_NONDEGENERATE;
    }

    if ((t1.type == INTERSECT_DEGENERATE && t1.N == 2) ||
        (t2.type == INTERSECT_DEGENERATE && t2.N == 2) ||
        (t3.type == INTERSECT_DEGENERATE && t3.N == 2)) {
        if (Nout && out) {
            s_point *coords = t1.N == 2 ? t1.coords : (t2.N == 2 ? t2.coords : t3.coords);
            out[0].x = coords[0].x; out[0].y = coords[0].y;
            out[1].x = coords[1].x; out[1].y = coords[1].y;
            *Nout = 2;
        }
        return INTERSECT_DEGENERATE;
    }

    if (Nout && out) {
        int h = 0;
        for (int ii=0; ii<t1.N; ii++) append_unique_2D_lim2(out, &h, (s_point2d){.x=t1.coords[ii].x, .y=t1.coords[ii].y}, tol_bdry);
        for (int ii=0; ii<t2.N; ii++) append_unique_2D_lim2(out, &h, (s_point2d){.x=t2.coords[ii].x, .y=t2.coords[ii].y}, tol_bdry);
        for (int ii=0; ii<t3.N; ii++) append_unique_2D_lim2(out, &h, (s_point2d){.x=t3.coords[ii].x, .y=t3.coords[ii].y}, tol_bdry);
        *Nout = h;
        assert(*Nout <= t1.N + t2.N + t3.N);
    }
    if (t1.N + t2.N + t3.N == 1) return INTERSECT_DEGENERATE;
    if (t1.N + t2.N + t3.N == 2) return INTERSECT_NONDEGENERATE;
    return INTERSECT_EMPTY;
}

e_intersect_type test_segment_triangle_intersect_2D(const s_point2d s[2], const s_point2d tri[3], double tol_bdry)
{
    return core_segment_triangle_intersect_2D(s, tri, tol_bdry, NULL, NULL);
}

s_segment_intersect segment_triangle_intersect_2D(const s_point2d s[2], const s_point2d tri[3], double tol_bdry)
{
    s_segment_intersect out;
    s_point2d coords[2];
    out.type = core_segment_triangle_intersect_2D(s, tri, tol_bdry, &out.N, coords);
    for (int ii=0; ii<out.N; ii++) {
        out.coords[ii].x = coords[ii].x;
        out.coords[ii].y = coords[ii].y;
        out.coords[ii].z = 0;
    }
    return out;
}


/* segment triangle 3D */
static void compute_segment_triangle_intersection_3D_nondegenerate(const s_point segment[2],
                                                                   const s_point triangle[3],
                                                                   int *Nout, s_point *out)
{   /* Moller's algorithm */
    s_point D = subtract_points(segment[1], segment[0]);
    s_point E1 = subtract_points(triangle[1], triangle[0]);
    s_point E2 = subtract_points(triangle[2], triangle[0]);
    s_point T = subtract_points(segment[0], triangle[0]);

    double denom = dot_prod(cross_prod(D, E2), E1);
    if (denom == 0 || !isfinite(denom)) { *Nout = 0; return; }
    double u = dot_prod(cross_prod(D, E2), T) / denom;
    double v = dot_prod(cross_prod(T, E1), D) / denom;

    *out = sum_points(sum_points(scale_point(triangle[0], 1-u-v), 
                                 scale_point(triangle[1], u)), 
                                 scale_point(triangle[2], v));
    *Nout = 1;
}

static e_intersect_type core_segment_triangle_intersect_3D_robust(const s_point segment[2],
                                                                  const s_point triangle[3],
                                                                  int *Nout, s_point *out)
{
    int o1 = test_orientation(triangle, segment[0]);
    int o2 = test_orientation(triangle, segment[1]);

    if (o1 == 0 && o2 == 0) {
        int drop = coord_to_drop_from_plane(triangle);
        s_point2d tri2[3] = { drop_to_2D(triangle[0], drop),
                              drop_to_2D(triangle[1], drop),
                              drop_to_2D(triangle[2], drop) };
        s_point2d seg2[2] = { drop_to_2D(segment[0], drop),
                              drop_to_2D(segment[1], drop) };

        s_point2d coords[2];
        int N2 = 0;
        e_intersect_type type = core_segment_triangle_intersect_2D(seg2, tri2, 0, &N2, coords);
        if (type == INTERSECT_ERROR) { if (Nout) *Nout = 0; return INTERSECT_ERROR; }
        if (type == INTERSECT_NONDEGENERATE || type == INTERSECT_DEGENERATE) {
            if (Nout && out) {
                *Nout = 0;
                for (int ii=0; ii<N2; ii++) {
                    out[ii] = lift_point_from_dropped_2D(triangle, drop, coords[ii]);
                    if (point_is_valid(out[ii])) (*Nout)++;
                }
            }
            return INTERSECT_DEGENERATE;
        }
        if (Nout) *Nout = 0;
        return INTERSECT_EMPTY;
    }
    else if (o1 == 0 || o2 == 0) {
        e_geom_test test = test_point_in_triangle_3D(triangle, (o1==0 ? segment[0] : segment[1]), 0);
        if (test == TEST_ERROR) { if (Nout) *Nout = 0; return INTERSECT_ERROR; }
        if (test == TEST_IN || test == TEST_BOUNDARY) {
            if (Nout && out) { out[0] = (o1==0 ? segment[0] : segment[1]); *Nout = 1; }
            return INTERSECT_DEGENERATE;
        }
        if (Nout) *Nout = 0;
        return INTERSECT_EMPTY;
    }

    if (o1 == o2) { if (Nout) *Nout = 0; return INTERSECT_EMPTY; }

    s_point a = triangle[0], b = triangle[1], c = triangle[2];
    s_point q1 = segment[0], q2 = segment[1];
    int s1 = test_orientation((s_point[3]){ a, b, q1 }, q2);
    int s2 = test_orientation((s_point[3]){ b, c, q1 }, q2);
    int s3 = test_orientation((s_point[3]){ c, a, q1 }, q2);

    if ( (s1 == 0 && s2 == s3) ||
         (s2 == 0 && s1 == s3) ||
         (s3 == 0 && s1 == s2) ||
         (s1 == 0 && s2 == 0) ||
         (s1 == 0 && s3 == 0) ||
         (s2 == 0 && s3 == 0) ) {
        if (Nout && out) compute_segment_triangle_intersection_3D_nondegenerate(segment, triangle, Nout, out);
        return INTERSECT_DEGENERATE;
    }
    if (s1 == s2 && s2 == s3) {
        if (Nout && out) compute_segment_triangle_intersection_3D_nondegenerate(segment, triangle, Nout, out);
        return INTERSECT_NONDEGENERATE;
    }
    if (Nout) *Nout = 0;
    return INTERSECT_EMPTY;
}

static e_intersect_type core_segment_triangle_intersect_3D(const s_point segment[2],
                                                           const s_point triangle[3],
                                                           double tol_bdry,
                                                           int *Nout, s_point out[2])
{
    if (Nout) *Nout = 0;
    if (tol_bdry == 0) return core_segment_triangle_intersect_3D_robust(segment, triangle, Nout, out);
    if (tri_is_degenerate(triangle[0], triangle[1], triangle[2])) return INTERSECT_ERROR;

    s_segment_intersect iplane = segment_plane_intersect(segment, triangle, tol_bdry);

    if (iplane.N == 1) {
        e_geom_test t = test_point_in_triangle_3D(triangle, iplane.coords[0], tol_bdry);
        if (t == TEST_IN || t == TEST_BOUNDARY) {
            if (Nout && out) { out[0] = iplane.coords[0]; *Nout = 1; }
            return INTERSECT_NONDEGENERATE;
        }
    }

    if (iplane.N == 2) {
        int drop = coord_to_drop_from_plane(triangle);
        s_point2d tri2[3] = { drop_to_2D(triangle[0], drop),
                              drop_to_2D(triangle[1], drop),
                              drop_to_2D(triangle[2], drop) };
        s_point2d seg2[2] = { drop_to_2D(iplane.coords[0], drop),
                              drop_to_2D(iplane.coords[1], drop) };

        s_point2d coords[2];
        int N2 = 0;
        e_intersect_type type = core_segment_triangle_intersect_2D(seg2, tri2, tol_bdry, &N2, coords);
        if (type == INTERSECT_ERROR) { if (Nout) *Nout = 0; return INTERSECT_ERROR; }
        if (type == INTERSECT_NONDEGENERATE || type == INTERSECT_DEGENERATE) {
            if (Nout && out) {
                *Nout = 0;
                for (int ii=0; ii<N2; ii++) {
                    out[ii] = lift_point_from_dropped_2D(triangle, drop, coords[ii]);
                    if (point_is_valid(out[ii])) (*Nout)++;
                    assert(*Nout <= 2);
                }
            }
            return type;
        }
    }

    if (Nout) *Nout = 0;
    return INTERSECT_EMPTY;
}

e_intersect_type test_segment_triangle_intersect_3D(const s_point segment[2], const s_point triangle[3], double tol_bdry)
{
    return core_segment_triangle_intersect_3D(segment, triangle, tol_bdry, NULL, NULL);
}

s_segment_intersect segment_triangle_intersect_3D(const s_point segment[2], const s_point triangle[3], double tol_bdry)
{
    s_segment_intersect out;
    out.type = core_segment_triangle_intersect_3D(segment, triangle, tol_bdry, &out.N, out.coords);
    return out;
}


/* Triangle tetrahedron intersection */
static e_intersect_type test_triangle_tetrahedron_intersect_robust(
    const s_point tri[3], const s_point tet[4])
{
    static const int tet_edges[6][2] = {{0,1},{0,2},{0,3},{1,2},{1,3},{2,3}};
    static const int tet_faces[4][3] = {{1,2,3},{0,3,2},{0,1,3},{0,2,1}};

    e_geom_test      vtx_tet[3];
    e_intersect_type tet_edg[6];
    e_intersect_type tri_fac[3][4];

    /* Strict checks first: short-circuit to NONDEGENERATE on first hit */
    for (int i = 0; i < 3; i++) {
        vtx_tet[i] = test_point_in_tetrahedron_robust(tet, tri[i]);
        if (vtx_tet[i] == TEST_IN) return INTERSECT_NONDEGENERATE;
    }

    for (int e = 0; e < 6; e++) {
        s_point seg[2] = {tet[tet_edges[e][0]], tet[tet_edges[e][1]]};
        tet_edg[e] = core_segment_triangle_intersect_3D_robust(seg, tri, NULL, NULL);
        if (tet_edg[e] == INTERSECT_NONDEGENERATE) return INTERSECT_NONDEGENERATE;
    }

    for (int i = 0; i < 3; i++) {
        s_point seg[2] = {tri[i], tri[(i+1)%3]};
        for (int f = 0; f < 4; f++) {
            s_point face[3] = {tet[tet_faces[f][0]], tet[tet_faces[f][1]], tet[tet_faces[f][2]]};
            tri_fac[i][f] = core_segment_triangle_intersect_3D_robust(seg, face, NULL, NULL);
            if (tri_fac[i][f] == INTERSECT_NONDEGENERATE) return INTERSECT_NONDEGENERATE;
        }
    }

    /* Boundary checks: short-circuit to DEGENERATE on first hit */
    for (int i = 0; i < 3; i++)
        if (vtx_tet[i] == TEST_BOUNDARY) return INTERSECT_DEGENERATE;

    for (int e = 0; e < 6; e++)
        if (tet_edg[e] == INTERSECT_DEGENERATE) return INTERSECT_DEGENERATE;

    for (int i = 0; i < 3; i++)
        for (int f = 0; f < 4; f++)
            if (tri_fac[i][f] == INTERSECT_DEGENERATE) return INTERSECT_DEGENERATE;

    return INTERSECT_EMPTY;
}

e_intersect_type test_triangle_tetrahedron_intersect(const s_point tri[3], const s_point tet[4])
{
    if (tri_is_degenerate(tri[0], tri[1], tri[2])) return INTERSECT_ERROR;
    if (tet_is_degenerate(tet))                    return INTERSECT_ERROR;
    return test_triangle_tetrahedron_intersect_robust(tri, tet);
}


