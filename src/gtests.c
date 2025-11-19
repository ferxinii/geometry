/* TOL: minimum distance between objects to be different. If TOL==0, all tests are ROBUST
 * EPS_degenerate: Avoid dividing by < EPS_degenerate, minimum area / lenght of object to be non-degenerate */

#include "points.h"
#include "gtests.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <stdatomic.h>
#include <float.h>


static atomic_flag predicates_init_flag = ATOMIC_FLAG_INIT;
static inline void ensure_predicates_initialized(void)
{
    if (!atomic_flag_test_and_set(&predicates_init_flag))
        exactinit();
}

int orientation_robust(const s_point p[3], s_point q)
{
    ensure_predicates_initialized();

    double aux = orient3d(p[0].coords, p[1].coords, p[2].coords, q.coords);
    if (aux > 0) return 1;
    else if (aux < 0) return -1;
    else return 0;
}

int insphere_robust(const s_point p[4], s_point q)
{   
    ensure_predicates_initialized();

    int o = orientation_robust(p, p[3]);
    int factor;
    if (o == 0) return 0;
    else if (o == 1) factor = 1;
    else factor = -1;

    double aux = insphere(p[0].coords, p[1].coords, p[2].coords, p[3].coords, q.coords);
    
    if (aux > 0) return factor;
    else if (aux < 0) return -factor;
    else return 0;
}


/* HELPERS */
static int sign(double x)
{
    if (x < 0) return -1;
    else if (x > 0) return 1;
    else return 0;
}

static int coord_to_drop_from_plane(const s_point plane[3])
{
    s_point AB = subtract_points(plane[1], plane[0]);
    s_point AC = subtract_points(plane[2], plane[0]);
    s_point n  = cross_prod(AB, AC);
    
    return coord_with_largest_component_3D(n);
}


static void drop_to_2D(const s_point p, int coord_to_drop, double out[2])
{
    int i1 = (coord_to_drop + 1) % 3;
    int i2 = (coord_to_drop + 2) % 3;
    out[0] = p.coords[i1];
    out[1] = p.coords[i2];
}

static s_point lift_point_from_dropped_2D(const s_point plane[3], int drop, const double in[2], double EPS_degenerate) 
{
    int i1 = (drop + 1) % 3;
    int i2 = (drop + 2) % 3;
    /* set kept coordinates from p2 */
    s_point out;
    out.coords[i1] = in[0];
    out.coords[i2] = in[1];
    /* try plane equation solving for missing coord */
    s_point n = cross_prod(subtract_points(plane[1], plane[0]), subtract_points(plane[2], plane[0]));
    if (fabs(n.coords[drop]) > EPS_degenerate) {
        double t0i = plane[0].coords[i1];
        double t0j = plane[0].coords[i2];
        double numer = - ( n.coords[i1]*(out.coords[i1]-t0i) + n.coords[i2]*(out.coords[i2]-t0j) );
        out.coords[drop] = plane[0].coords[drop] + numer / n.coords[drop];
        return out;
    }
    return point_NAN;
 }

static double point_segment_dist2_2D(const double p[2], const double a[2], const double b[2], double EPS_degenerate)
{
    double vx=b[0]-a[0], vy=b[1]-a[1];
    double wx=p[0]-a[0], wy=p[1]-a[1];
    

    double denom = vx*vx + vy*vy;
    if (fabs(denom) < EPS_degenerate) {
        return (p[0]-a[0])*(p[0]-a[0]) + (p[1]-a[1])*(p[1]-a[1]);
    } 
    double t = (wx*vx + wy*vy) / denom;

    if (t <= 0.0) return (p[0]-a[0])*(p[0]-a[0]) + (p[1]-a[1])*(p[1]-a[1]);
    else if (t >= 1.0) return (p[0]-b[0])*(p[0]-b[0]) + (p[1]-b[1])*(p[1]-b[1]);
    else {
        double projx = a[0] + t*vx;
        double projy = a[1] + t*vy;
        double dx = p[0] - projx, dy = p[1] - projy;
        return dx*dx + dy*dy;
    }
}

static double area_triangle_2D(const double a[2], const double b[2], const double c[2]) 
{
    return 0.5 * fabs(((a[0]-c[0])*(b[1]-a[1]) + (a[0]-b[0])*(c[1]-a[1])));
}

static int points_close_2D(const double a[2], const double b[2], double TOL) 
{
    if (TOL == 0) {  /* Robust branch */
        if (a[0] == b[0] && a[1] == b[1]) return 1;
        else return 0;
    } else {
        if ((a[0]-b[0])*(a[0]-b[0]) + (a[1]-b[1])*(a[1]-b[1]) <= TOL*TOL) return 1;
        else return 0;
    }
}

static int append_unique_2D_lim2(double out[4], int *count, const double p[2], double TOL)
{
    if ((*count)<2 && TOL > 0) {
        for (int i=0; i<(*count); i++) {
            double q[2] = { out[i*2+0], out[i*2+1] };
            if (points_close_2D(q, p, TOL)) return 0;
        }
    }
    if (*count >= 2) return 0;  /* already full (safety) */
    if (out) {
        out[(*count)*2+0] = p[0];
        out[(*count)*2+1] = p[1];
    }
    (*count)++;
    return 1;
}


/* GEOMETRICAL TESTS */
e_geom_test test_point_in_interval_1D(double x, double a, double b, double EPS_degenerate, double TOL)
{
    if (fabs(b - a) < EPS_degenerate) return TEST_DEGENERATE;
    if (fabs(x - a) <= TOL|| fabs(x - b) <= TOL) return TEST_BOUNDARY;
    if (a > b) { double tmp = a; a = b; b = tmp; }
    if (x + TOL>= a && x - TOL<= b) return TEST_IN;
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


e_geom_test test_point_in_triangle_2D(const double a[2], const double b[2], const double c[2], const double p[2], double EPS_degenerate, double TOL)
{
    if (area_triangle_2D(a, b, c) < EPS_degenerate) return TEST_DEGENERATE;

    /* Robust orientation tests of p vs each edge */
    int o1 = sign(orient2d(a, b, p));
    int o2 = sign(orient2d(b, c, p));
    int o3 = sign(orient2d(c, a, p));

    if (TOL == 0) return test_point_in_triangle_2D_robust_from_orientations(o1, o2, o3);

    /* Numerical tolerance branch */
    double TOL2 = TOL*TOL;
    if (point_segment_dist2_2D(p, a, b, EPS_degenerate) <= TOL2) return TEST_BOUNDARY;
    if (point_segment_dist2_2D(p, b, c, EPS_degenerate) <= TOL2) return TEST_BOUNDARY;
    if (point_segment_dist2_2D(p, c, a, EPS_degenerate) <= TOL2) return TEST_BOUNDARY;
    
    /* If not on boundary, check robust orientation to determine TEST_IN/TEST_OUT */
    return test_point_in_triangle_2D_robust_from_orientations(o1, o2, o3);
}


static e_geom_test test_point_in_triangle_3D_robust(const s_point triangle[3], s_point p)
{
    if (orientation_robust(triangle, p) != 0) return TEST_OUT;

    /* Point coplanar. Project to 2D and check */   
    int drop = coord_to_drop_from_plane(triangle);
    double A2[2]; drop_to_2D(triangle[0], drop, A2);
    double B2[2]; drop_to_2D(triangle[1], drop, B2);
    double C2[2]; drop_to_2D(triangle[2], drop, C2);
    double p2[2]; drop_to_2D(p, drop, p2);

    return test_point_in_triangle_2D(A2, B2, C2, p2, 0, 0);
}


e_geom_test test_point_in_triangle_3D(const s_point triangle[3], s_point p, double EPS_degenerate, double TOL) 
{
    if (area_triangle(triangle) < EPS_degenerate) return TEST_DEGENERATE;

    if (TOL == 0) return test_point_in_triangle_3D_robust(triangle, p);

    s_point closest = project_point_to_plane(p, triangle, EPS_degenerate);
    if (!point_is_valid(closest)) return TEST_ERROR;
    if (distance_squared(closest, p) > TOL*TOL) return TEST_OUT;
    
    /* Point is TOL close to plane. Project to 2D and check */
    int drop = coord_to_drop_from_plane(triangle);
    double A2[2]; drop_to_2D(triangle[0], drop, A2);
    double B2[2]; drop_to_2D(triangle[1], drop, B2);
    double C2[2]; drop_to_2D(triangle[2], drop, C2);
    double p2[2]; drop_to_2D(closest, drop, p2);
    return test_point_in_triangle_2D(A2, B2, C2, p2, EPS_degenerate, TOL);
}


static e_geom_test test_point_in_tetrahedron_robust(const s_point tetra[4], s_point query)
{   /* No need to be properly oriented */
    s_point tmp[3];

    /* Reference signs: e, query signs: s */
    tmp[0] = tetra[1];   tmp[1] = tetra[2];   tmp[2] = tetra[3];
    int e0 = orientation_robust(tmp, tetra[0]);
    int s0 = orientation_robust(tmp, query);

    tmp[0] = tetra[0];   tmp[1] = tetra[3];   tmp[2] = tetra[2];
    int e1 = orientation_robust(tmp, tetra[1]);
    int s1 = orientation_robust(tmp, query);

    tmp[0] = tetra[0];   tmp[1] = tetra[1];   tmp[2] = tetra[3];
    int e2 = orientation_robust(tmp, tetra[2]);
    int s2 = orientation_robust(tmp, query);

    tmp[0] = tetra[0];   tmp[1] = tetra[2];   tmp[2] = tetra[1];
    int e3 = orientation_robust(tmp, tetra[3]);
    int s3 = orientation_robust(tmp, query);

    if (e0 == 0 || e1 == 0 || e2 == 0 || e3 == 0) return TEST_DEGENERATE;

    if ((s0 != 0 && s0 * e0 < 0) ||
        (s1 != 0 && s1 * e1 < 0) ||
        (s2 != 0 && s2 * e2 < 0) ||
        (s3 != 0 && s3 * e3 < 0))
        return TEST_OUT;
    

    if (s0 == 0 || s1 == 0 || s2 == 0 || s3 == 0) return TEST_BOUNDARY;
    return TEST_IN;
}


e_geom_test test_point_in_tetrahedron(const s_point tetra[4], s_point query, double EPS_degenerate, double TOL)
{   
    if (signed_volume_tetra(tetra) < EPS_degenerate) return TEST_DEGENERATE;

    if (TOL == 0) return test_point_in_tetrahedron_robust(tetra, query);

    /* First check if query is EPS-in each face */
    s_point tmp[3];
    tmp[0] = tetra[1];   tmp[1] = tetra[2];   tmp[2] = tetra[3];
    e_geom_test t1 = test_point_in_triangle_3D(tmp, query, EPS_degenerate, TOL);
    
    tmp[0] = tetra[0];   tmp[1] = tetra[3];   tmp[2] = tetra[2];
    e_geom_test t2 = test_point_in_triangle_3D(tmp, query, EPS_degenerate, TOL);

    tmp[0] = tetra[0];   tmp[1] = tetra[1];   tmp[2] = tetra[3];
    e_geom_test t3 = test_point_in_triangle_3D(tmp, query, EPS_degenerate, TOL);

    tmp[0] = tetra[0];   tmp[1] = tetra[2];   tmp[2] = tetra[1];
    e_geom_test t4 = test_point_in_triangle_3D(tmp, query, EPS_degenerate, TOL);

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


s_points_test test_points_in_halfspace(const s_point plane_ordered[3], const s_points *points, double EPS_degenerate, double TOL, e_geom_test out_buff[points->N])
{   /* Right hand rule: normal outwards */
    if (out_buff == NULL) out_buff = malloc(points->N * sizeof(e_geom_test));
    int Nin = 0, Nbdy = 0, Nout = 0;

    if (TOL == 0) {  /* Robust branch */
        for (int ii=0; ii<points->N; ii++) {
            int o = orientation_robust(plane_ordered, points->p[ii]);
            if (o == 1) { out_buff[ii] = TEST_IN; Nin++; } 
            else if (o == 0) { out_buff[ii] = TEST_BOUNDARY; Nbdy++; }
            else { out_buff[ii] = TEST_OUT; Nout++; }
        } 
    } else {  /* Non-robust branch */
        s_point n = cross_prod(subtract_points(plane_ordered[1], plane_ordered[0]),
                               subtract_points(plane_ordered[2], plane_ordered[0]));
        n = normalize_vec(n, EPS_degenerate);
        if (!point_is_valid(n)) {
            for (int ii=0; ii<points->N; ii++) out_buff[ii] = TEST_ERROR;
            return (s_points_test){ .Nin = 0, .Nbdy = 0, .Nout = 0, .Nerr = points->N, .indicator = out_buff };
        }
        double d = dot_prod(n, plane_ordered[0]);

        for (int ii=0; ii<points->N; ii++) {
            double s = dot_prod(points->p[ii], n) - d;
            if (s < -TOL) { out_buff[ii] = TEST_IN; Nin++; } 
            else if (fabs(s) <= TOL) { out_buff[ii] = TEST_BOUNDARY; Nbdy++; } 
            else { out_buff[ii] = TEST_OUT; Nout++; }
        }
    }
    
    return (s_points_test){ .Nin = Nin, .Nbdy = Nbdy, .Nout = Nout, .Nerr = 0, .indicator = out_buff };
}




/* INTERSECTIONS */
/* Segment segment intersection 2D */
static void compute_segment_segment_intersection_2D_nondegenerate(const double A1[2], const double A2[2], const double B1[2], const double B2[2], double EPS_degenerate, int *Nout, double *out)
{
    double A_a = A2[1] - A1[1];
    double A_b = A1[0] - A2[0];
    double A_c = A_a*A1[0] + A_b*A1[1];
    double B_a = B2[1] - B1[1];
    double B_b = B1[0] - B2[0];
    double B_c = B_a*B1[0] + B_b*B1[1];

    double det = A_a*B_b - B_a*A_b;
    if (fabs(det) <= EPS_degenerate) *Nout = 0;  /* nearly parallel */
    else {
        out[0] = (B_b*A_c - A_b*B_c)/det;
        out[1] = (A_a*B_c - B_a*A_c)/det;
        *Nout = 1;
    }
}

static void compute_segment_segment_intersection_2D_colinear(const double A1[2], const double A2[2], const double B1[2], const double B2[2], double EPS_degenerate, double TOL, int *Nout, double *out)
{
    double dx = A2[0] - A1[0], dy = A2[1] - A1[1];
    double denom = dx*dx + dy*dy;
    if (denom <= EPS_degenerate) { *Nout = 0; return; }

    /* Project B endpoints onto A's parameter t ∈ R */
    double tB1 = ((B1[0]-A1[0])*dx + (B1[1]-A1[1])*dy) / denom;
    double tB2 = ((B2[0]-A1[0])*dx + (B2[1]-A1[1])*dy) / denom;

    /* Interval of intersection in parametric t */
    double tmin = fmin(tB1, tB2), tmax = fmax(tB1, tB2);
    double lo = fmax(0.0, tmin), hi = fmin(1.0, tmax);

    if (hi + TOL < lo - TOL) { *Nout = 0; return; }  /* Case 1: No overlap */
    if (fabs(hi - lo) <= TOL) {  /* Case 2: Touching at exactly one point */
        double t = 0.5 * (hi + lo);
        out[0] = A1[0] + t * dx;
        out[1] = A1[1] + t * dy;
        *Nout = 1;
        return;
    }
    /* Case 3: Proper interval (2 points) */
    if (lo < 0.0) lo = 0.0;
    if (hi > 1.0) hi = 1.0;
    out[0] = A1[0] + lo * dx;
    out[1] = A1[1] + lo * dy;
    out[2] = A1[0] + hi * dx;
    out[3] = A1[1] + hi * dy;
    *Nout = 2;
}


static e_intersect_type core_segment_segment_intersect_2D(const double A1[2], const double A2[2], const double B1[2], const double B2[2], double EPS_degenerate, double TOL, int *Nout, double *out)
{
    /* 1) Check degeneracy (robust and non-robust branches) */
    double TOL2 = TOL*TOL;
    int A_is_point = points_close_2D(A1, A2, EPS_degenerate);
    int B_is_point = points_close_2D(B1, B2, EPS_degenerate);
    if (A_is_point && B_is_point) {
        if (points_close_2D(A1, B1, TOL)) {
            if (Nout && out) {  /* Return average of points */
                out[0] = (A1[0]+A2[0]+B1[0]+B2[0])/4.0; 
                out[1] = (A1[1]+A2[1]+B1[1]+B2[1])/4.0; 
                *Nout = 1;
            }
            return INTERSECT_DEGENERATE;
        } 
        if (Nout) *Nout = 0;
        return INTERSECT_EMPTY;
    }
    if (A_is_point) {
        if ((TOL==0) ? (orient2d(B1,B2,A1) == 0 &&
                        test_point_in_interval_1D(A1[0], B1[0], B2[0], EPS_degenerate, TOL) != TEST_OUT &&
                        test_point_in_interval_1D(A1[1], B1[1], B2[1], EPS_degenerate, TOL) != TEST_OUT) 
                     : (point_segment_dist2_2D(A1,B1,B2,EPS_degenerate) <= TOL2)) {
            if (Nout && out) { out[0] = A1[0]; out[1] = A1[1]; *Nout = 1; }
            return INTERSECT_DEGENERATE;
        }
        if (Nout) *Nout = 0;
        return INTERSECT_EMPTY;
    }
    if (B_is_point) {
        if ((TOL==0) ? (orient2d(A1,A2,B1)==0 &&
                        test_point_in_interval_1D(B1[0], A1[0], A2[0], EPS_degenerate, TOL) != TEST_OUT &&
                        test_point_in_interval_1D(B1[1], A1[1], A2[1], EPS_degenerate, TOL) != TEST_OUT) 
                     : (point_segment_dist2_2D(B1,A1,A2,EPS_degenerate) <= TOL2)) {
            if (Nout && out) { out[0] = B1[0]; out[1] = B1[1]; *Nout = 1; }
            return INTERSECT_DEGENERATE;
        }
        if (Nout) *Nout = 0;
        return INTERSECT_EMPTY;
    }


    /* 2) Segments are non-degenerate */
    double d1 = (TOL==0) ? orient2d(A1, A2, B1) : point_segment_dist2_2D(B1, A1, A2, EPS_degenerate); 
    double d2 = (TOL==0) ? orient2d(A1, A2, B2) : point_segment_dist2_2D(B2, A1, A2, EPS_degenerate); 
    double d3 = (TOL==0) ? orient2d(B1, B2, A1) : point_segment_dist2_2D(A1, B1, B2, EPS_degenerate); 
    double d4 = (TOL==0) ? orient2d(B1, B2, A2) : point_segment_dist2_2D(A2, B1, B2, EPS_degenerate); 
    
    /* 2a) Check if colinear (robust) or if both overlap (non-robust) */
    if ( ((TOL==0) ? (d1 == 0 && d2 == 0 && d3 == 0 && d4 == 0) :
                     (d1 <= TOL2 && d2 <= TOL2 && d3 <= TOL2 && d4 <= TOL2))) {
        double min_x = fmax(fmin(A1[0], A2[0]), fmin(B1[0], B2[0]));
        double max_x = fmin(fmax(A1[0], A2[0]), fmax(B1[0], B2[0]));
        double min_y = fmax(fmin(A1[1], A2[1]), fmin(B1[1], B2[1]));
        double max_y = fmin(fmax(A1[1], A2[1]), fmax(B1[1], B2[1]));
        if (min_x <= max_x + TOL && min_y <= max_y + TOL) {
            if (Nout && out) compute_segment_segment_intersection_2D_colinear(A1, A2, B1, B2, EPS_degenerate, TOL, Nout, out);
            return INTERSECT_DEGENERATE;  
        } else {
            if (Nout) *Nout = 0;
            return INTERSECT_EMPTY;
        }
    }

    /* 2b) Single end colinear with segment (check if it lies within the other segment) */
    if ( (TOL==0) ? (d1==0 && 
                     test_point_in_interval_1D(B1[0], A1[0], A2[0], EPS_degenerate, TOL) != TEST_OUT &&
                     test_point_in_interval_1D(B1[1], A1[1], A2[1], EPS_degenerate, TOL) != TEST_OUT) 
                  : (d1<=TOL2)) {
        if (Nout && out) { out[0] = B1[0]; out[1] = B1[1]; *Nout = 1; }
        return INTERSECT_DEGENERATE;
    }
    if ( (TOL==0) ? (d2==0 && 
                     test_point_in_interval_1D(B2[0], A1[0], A2[0], EPS_degenerate, TOL) != TEST_OUT &&
                     test_point_in_interval_1D(B2[1], A1[1], A2[1], EPS_degenerate, TOL) != TEST_OUT) 
                  : (d2<=TOL2)) {
        if (Nout && out) { out[0] = B2[0]; out[1] = B2[1]; *Nout = 1; }
        return INTERSECT_DEGENERATE;
    }
    if ( (TOL==0) ? (d3==0 &&
                     test_point_in_interval_1D(A1[0], B1[0], B2[0], EPS_degenerate, TOL) != TEST_OUT &&
                     test_point_in_interval_1D(A1[1], B1[1], B2[1], EPS_degenerate, TOL) != TEST_OUT) 
                  : (d3<=TOL2)) {
        if (Nout && out) { out[0] = A1[0]; out[1] = A1[1]; *Nout = 1; }
        return INTERSECT_DEGENERATE;
    }
    if ( (TOL==0) ? (d4==0 &&
                     test_point_in_interval_1D(A2[0], B1[0], B2[0], EPS_degenerate, TOL) != TEST_OUT &&
                     test_point_in_interval_1D(A2[1], B1[1], B2[1], EPS_degenerate, TOL) != TEST_OUT)
                  : (d4<=TOL2)) {
        if (Nout && out) { out[0] = A2[0]; out[1] = A2[1]; *Nout = 1; }
        return INTERSECT_DEGENERATE;
    }

    /* 2c) Fallback to robust detection of interior using two stradling tests */
    if (TOL != 0) {  /* Else, di already is orient2d */
        d1 = orient2d(A1, A2, B1);
        d2 = orient2d(A1, A2, B2);
        d3 = orient2d(B1, B2, A1);
        d4 = orient2d(B1, B2, A2);
    }
    if ( ((d1>0 && d2<0) || (d1<0 && d2>0)) &&
         ((d3>0 && d4<0) || (d3<0 && d4>0)) ) {
        if (Nout && out) compute_segment_segment_intersection_2D_nondegenerate(A1, A2, B1, B2, EPS_degenerate, Nout, out);
        return INTERSECT_NONDEGENERATE;
    }

    if (Nout) *Nout = 0;
    return INTERSECT_EMPTY;
}

e_intersect_type test_segment_segment_intersect_2D(const double A1[2], const double A2[2], const double B1[2], const double B2[2], double EPS_degenerate, double TOL)
{
    return core_segment_segment_intersect_2D(A1, A2, B1, B2, EPS_degenerate, TOL, NULL, NULL);
}

s_segment_intersect segment_segment_intersect_2D(const double A1[2], const double A2[2], const double B1[2], const double B2[2], double EPS_degenerate, double TOL)
{
    s_segment_intersect out;  /* out.coords is s_point*, not double* */
    double coords[4];
    out.type = core_segment_segment_intersect_2D(A1, A2, B1, B2, EPS_degenerate, TOL, &out.N, coords);
    for (int ii=0; ii<out.N; ii++) {
        out.coords[ii].x = coords[ii*2+0];
        out.coords[ii].y = coords[ii*2+1];
        out.coords[ii].z = 0;
    }
    return out;
}



/* segment plane intersection */
static e_intersect_type core_segment_plane_intersect_robust(const s_point seg[2], const s_point plane[3], double EPS_degenerate, int *Nout, s_point out[2])
{
    int o0 = orientation_robust(plane, seg[0]);
    int o1 = orientation_robust(plane, seg[1]);

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
        if (fabs(sA - sB) < EPS_degenerate) { *Nout = 0; return INTERSECT_NONDEGENERATE; }
        double t = sA / (sA - sB);
        out[0] = interpolate_points(seg[0], seg[1], t); 
        *Nout = 1;
    }
    return INTERSECT_NONDEGENERATE;
}

static e_intersect_type core_segment_plane_intersect(const s_point seg[2], const s_point plane[3], double EPS_degenerate, double TOL, int *Nout, s_point out[2])
{
    if (TOL == 0) return core_segment_plane_intersect_robust(seg, plane, EPS_degenerate, Nout, out);

    s_point n = cross_prod(subtract_points(plane[1], plane[0]),
                           subtract_points(plane[2], plane[0]));
    double d = dot_prod(n, plane[0]);
    double sA = dot_prod(seg[0], n) - d;
    double sB = dot_prod(seg[1], n) - d;

    if (sA * sB > TOL*TOL) {  /* both endpoints on same side (and not near plane) */
        if (Nout && out) *Nout = 0;
        return INTERSECT_EMPTY;
    }

    if (fabs(sA) < TOL && fabs(sB) < TOL) {  /* whole segment lies in plane */
        if (Nout && out) {
            out[0] = project_point_to_plane(seg[0], plane, EPS_degenerate);
            out[1] = project_point_to_plane(seg[1], plane, EPS_degenerate);
            *Nout = point_is_valid(out[0]) + point_is_valid(out[1]); 
        }
        return INTERSECT_DEGENERATE;
    }
    if (fabs(sA) < TOL || fabs(sB) < TOL) {  /* Only one end close to plane */
        if (Nout && out) {
            out[0] = project_point_to_plane((fabs(sA)<TOL ? seg[0] : seg[1]), plane, EPS_degenerate);
            *Nout = point_is_valid(out[0]);
        }
        return INTERSECT_DEGENERATE;
    }

    if (Nout && out) {  /* find parametric intersection */
        double denom = sA - sB;
        if (fabs(denom) < EPS_degenerate) { *Nout = 0; return INTERSECT_NONDEGENERATE; }
        double t = sA / denom;
        if (t < -TOL || t > 1 + TOL) { *Nout = 0; return INTERSECT_NONDEGENERATE; }
        if (t < 0) t = 0;
        if (t > 1) t = 1;
        out[0] = interpolate_points(seg[0], seg[1], t);
        *Nout = 1;
    }
    return INTERSECT_NONDEGENERATE;
}

e_intersect_type test_segment_plane_intersect(const s_point seg[2], const s_point plane[3], double EPS_degenerate, double TOL)
{
    return core_segment_plane_intersect(seg, plane, EPS_degenerate, TOL, NULL, NULL);
}

s_segment_intersect segment_plane_intersect(const s_point seg[2], const s_point plane[3], double EPS_degenerate, double TOL)
{
    s_segment_intersect out;
    out.type = core_segment_plane_intersect(seg, plane, EPS_degenerate, TOL, &out.N, out.coords);
    return out;
}



/* Segment triangle intersection */
e_intersect_type core_segment_triangle_intersect_2D(const double S1[2], const double S2[2], const double A[2], const double B[2], const double C[2], double EPS_degenerate, double TOL, int *Nout, double *out)
{
    e_geom_test i1 = test_point_in_triangle_2D(A, B, C, S1, EPS_degenerate, TOL);
    e_geom_test i2 = test_point_in_triangle_2D(A, B, C, S2, EPS_degenerate, TOL);
    if (i1 == TEST_ERROR || i2 == TEST_ERROR || i1 == TEST_DEGENERATE || i2 == TEST_DEGENERATE) return INTERSECT_ERROR;
    if (i1 == TEST_BOUNDARY && i2 == TEST_BOUNDARY) {  /* Both on boundary of triangle */
        if (Nout && out) {
            out[0] = S1[0]; out[1] = S1[1];
            out[2] = S2[0]; out[3] = S2[1];
            *Nout = 2;
        }
        return INTERSECT_DEGENERATE;
    }
    if ( (i1 == TEST_IN || i1 == TEST_BOUNDARY) && (i2 == TEST_IN || i2 == TEST_BOUNDARY)) {
        if (Nout && out) {
            out[0] = S1[0]; out[1] = S1[1];
            out[2] = S2[0]; out[3] = S2[1];
            *Nout = 2;
        }
        return INTERSECT_NONDEGENERATE;
    }

    /* At least one point is outside. Check intersections with triangle edges */
    s_segment_intersect t1 = segment_segment_intersect_2D(S1, S2, A, B, EPS_degenerate, TOL);
    s_segment_intersect t2 = segment_segment_intersect_2D(S1, S2, B, C, EPS_degenerate, TOL);
    s_segment_intersect t3 = segment_segment_intersect_2D(S1, S2, C, A, EPS_degenerate, TOL);

    /* One point inside + one intersection */
    if ((i1 == TEST_IN || i2 == TEST_IN) && (t1.N > 0 || t2.N > 0 || t3.N > 0)) {
        if (Nout && out) {
            out[0] = (i1 == TEST_IN) ? S1[0] : S2[0]; out[1] = (i1 == TEST_IN) ? S1[1] : S2[1];
            s_point *coords = t1.N > 0 ? t1.coords : (t2.N > 0 ? t2.coords : t3.coords);
            out[2] = coords[0].x; out[3] = coords[0].y;
        }
        return INTERSECT_NONDEGENERATE;
    }
    
    /* Two intersections on same triangle side... colinear with side */
    if ((t1.type == INTERSECT_DEGENERATE && t1.N == 2) || 
        (t2.type == INTERSECT_DEGENERATE && t2.N == 2) ||
        (t3.type == INTERSECT_DEGENERATE && t3.N == 2)) {
        if (Nout && out) {
            s_point *coords = t1.N == 2 ? t1.coords : (t2.N == 2 ? t2.coords : t3.coords);
            out[0] = coords[0].x; out[1] = coords[0].y;
            out[2] = coords[1].x; out[3] = coords[1].y;
            *Nout = 2;
        }
        return INTERSECT_DEGENERATE;
    }

    /* Two segments outside */ 
    if (Nout && out) {
        int h = 0;
        if (t1.N > 0 || t2.N > 0 || t3.N > 0) {
            for (int ii=0; ii<t1.N; ii++) append_unique_2D_lim2(out, &h, t1.coords[ii].coords, TOL);
            for (int ii=0; ii<t2.N; ii++) append_unique_2D_lim2(out, &h, t2.coords[ii].coords, TOL);
            for (int ii=0; ii<t3.N; ii++) append_unique_2D_lim2(out, &h, t3.coords[ii].coords, TOL);
        }
        *Nout = h;
        assert(*Nout <= t1.N + t2.N + t3.N);
    }
    if (t1.N + t2.N + t3.N == 1) return INTERSECT_DEGENERATE;
    if (t1.N + t2.N + t3.N == 2) return INTERSECT_NONDEGENERATE;
    return INTERSECT_EMPTY;
}

e_intersect_type test_segment_triangle_intersect_2D(const double S1[2], const double S2[2], const double A[2], const double B[2], const double C[2], double EPS_degenerate, double TOL)
{
    return core_segment_triangle_intersect_2D(S1, S2, A, B, C, EPS_degenerate, TOL, NULL, NULL);
}

s_segment_intersect segment_triangle_intersect_2D(const double S1[2], const double S2[2], const double A[2], const double B[2], const double C[2], double EPS_degenerate, double TOL)
{
    s_segment_intersect out;  /* out.coords is s_point*, not double* */
    double coords[4];
    out.type = core_segment_triangle_intersect_2D(S1, S2, A, B, C, EPS_degenerate, TOL, &out.N, coords);
    for (int ii=0; ii<out.N; ii++) {
        out.coords[ii].x = coords[ii*2+0];
        out.coords[ii].y = coords[ii*2+1];
        out.coords[ii].z = 0;
    }
    return out;
}


static void compute_segment_triangle_intersection_3D_nondegenerate(const s_point segment[2], const s_point triangle[3], double EPS_degenerate, int *Nout, s_point *out)
{   /* Moller's algorithm */
    s_point D = subtract_points(segment[1], segment[0]);
    s_point E1 = subtract_points(triangle[1], triangle[0]);
    s_point E2 = subtract_points(triangle[2], triangle[0]);
    s_point T = subtract_points(segment[0], triangle[0]);

    double denom = dot_prod(cross_prod(D, E2), E1);
    if (fabs(denom) < EPS_degenerate) { *Nout = 0; return; }
    double u = dot_prod(cross_prod(D, E2), T) / denom;
    double v = dot_prod(cross_prod(T, E1), D) / denom;

    *out = sum_points(sum_points(scale_point(triangle[0], 1-u-v), 
                                 scale_point(triangle[1], u)), 
                                 scale_point(triangle[2], v));
    *Nout = 1;
}

static e_intersect_type core_segment_triangle_intersect_3D_robust(const s_point segment[2], const s_point triangle[3], double EPS_degenerate, int *Nout, s_point *out)
{   
    /* Check if any end of segment coplanar with triangle */
    int o1 = orientation_robust(triangle, segment[0]);
    int o2 = orientation_robust(triangle, segment[1]);

    if (o1 == 0 && o2 == 0) {  /* Whole segment complanar with triangle */
        int drop = coord_to_drop_from_plane(triangle);
        double A2[2]; drop_to_2D(triangle[0], drop, A2);
        double B2[2]; drop_to_2D(triangle[1], drop, B2);
        double C2[2]; drop_to_2D(triangle[2], drop, C2);
        double P2[2]; drop_to_2D(segment[0], drop, P2);
        double D2[2]; drop_to_2D(segment[1], drop, D2);
        
        s_segment_intersect i2D = segment_triangle_intersect_2D(P2, D2, A2, B2, C2, EPS_degenerate, 0);
        if (i2D.type == INTERSECT_ERROR) return INTERSECT_ERROR;
        if (i2D.type == INTERSECT_NONDEGENERATE || i2D.type == INTERSECT_DEGENERATE) {
            if (Nout && out) {
                *Nout = 0;
                for (int ii=0; ii<i2D.N; ii++) {
                    out[ii] = lift_point_from_dropped_2D(triangle, drop, i2D.coords[ii].coords, EPS_degenerate);
                    if (point_is_valid(out[ii])) (*Nout)++;
                }
            }
            return INTERSECT_DEGENERATE;
        }
        if (Nout) *Nout = 0;
        return INTERSECT_EMPTY; 
    }
    else if (o1 == 0 || o2 == 0) {  /* A single end is coplanar */
        e_geom_test test = test_point_in_triangle_3D(triangle, (o1==0 ? segment[0] : segment[1]), EPS_degenerate, 0);
        if (test == TEST_ERROR) return INTERSECT_ERROR;
        if (test == TEST_IN || test == TEST_BOUNDARY) {
            if (Nout && out) { out[0] = (o1==0 ? segment[0] : segment[1]); *Nout = 1; }
            return INTERSECT_DEGENERATE;
        }
        if (Nout) *Nout = 0;
        return INTERSECT_EMPTY;
    }

    /* Check that segment crosses triangle's plane */
    if (o1 == o2) {
        if (Nout) *Nout = 0;
        return INTERSECT_EMPTY;
    }

    /* Segura's algorithm */
    s_point a = triangle[0], b = triangle[1], c = triangle[2];
    s_point q1 = segment[0], q2 = segment[1];
    int s1 = orientation_robust((s_point[3]){ a, b, q1 }, q2);
    int s2 = orientation_robust((s_point[3]){ b, c, q1 }, q2);
    int s3 = orientation_robust((s_point[3]){ c, a, q1 }, q2);

    if ( (s1 == 0 && s2 == s3) ||
         (s2 == 0 && s1 == s3) ||
         (s3 == 0 && s1 == s2) ||
         (s1 == 0 && s2 == 0) ||
         (s1 == 0 && s3 == 0) ||
         (s2 == 0 && s3 == 0) ) {  /* Intersection at edge */
        if (Nout && out) compute_segment_triangle_intersection_3D_nondegenerate(segment, triangle, EPS_degenerate, Nout, out);
        return INTERSECT_DEGENERATE;
    }
    if (s1 == s2 && s2 == s3 ) {
        if (Nout && out) compute_segment_triangle_intersection_3D_nondegenerate(segment, triangle, EPS_degenerate, Nout, out);
        return INTERSECT_NONDEGENERATE;  
    }
    if (Nout) *Nout = 0;
    return INTERSECT_EMPTY;  
}

static e_intersect_type core_segment_triangle_intersect_3D(const s_point segment[2], const s_point triangle[3], double EPS_degenerate, double TOL, int *Nout, s_point out[2])
{
    if (area_triangle(triangle) < EPS_degenerate) return 0;

    if (TOL == 0) return core_segment_triangle_intersect_3D_robust(segment, triangle, EPS_degenerate, Nout, out);

    /* Non-robust branch */
    s_segment_intersect iplane = segment_plane_intersect(segment, triangle, EPS_degenerate, TOL);
    
    if (iplane.N == 1) {
        e_geom_test t = test_point_in_triangle_3D(triangle, iplane.coords[0], EPS_degenerate, TOL);
        if (t == TEST_IN || t == TEST_BOUNDARY) {
            if (Nout && out) {
                out[0] = iplane.coords[0];
                *Nout = 1;
            }
            return INTERSECT_NONDEGENERATE;
        }
    }

    if (iplane.N == 2) {
        int drop = coord_to_drop_from_plane(triangle);
        double A2[2]; drop_to_2D(triangle[0], drop, A2);
        double B2[2]; drop_to_2D(triangle[1], drop, B2);
        double C2[2]; drop_to_2D(triangle[2], drop, C2);
        double P2[2]; drop_to_2D(iplane.coords[0], drop, P2);
        double D2[2]; drop_to_2D(iplane.coords[1], drop, D2);
        s_segment_intersect i2D = segment_triangle_intersect_2D(P2, D2, A2, B2, C2, EPS_degenerate, TOL);
        if (i2D.type == INTERSECT_ERROR) return INTERSECT_ERROR;
        if (i2D.type == INTERSECT_NONDEGENERATE || i2D.type == INTERSECT_DEGENERATE) {
            if (Nout && out) {
                *Nout = 0;
                for (int ii=0; ii<i2D.N; ii++) {
                    out[ii] = lift_point_from_dropped_2D(triangle, drop, i2D.coords[ii].coords, EPS_degenerate);
                    if (point_is_valid(out[ii])) (*Nout)++;
                }
            }
            return INTERSECT_DEGENERATE;
        }
    }

    if (Nout) *Nout = 0;
    return INTERSECT_EMPTY;
}

e_intersect_type test_segment_triangle_intersect_3D(const s_point segment[2], const s_point triangle[3], double EPS_degenerate, double TOL)
{
    return core_segment_triangle_intersect_3D(segment, triangle, EPS_degenerate, TOL, NULL, NULL);
}

s_segment_intersect segment_triangle_intersect_3D(const s_point segment[2], const s_point triangle[3], double EPS_degenerate, double TOL)
{
    s_segment_intersect out; 
    out.type = core_segment_triangle_intersect_3D(segment, triangle, EPS_degenerate, TOL, &out.N, out.coords);
    return out;
}


