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


static atomic_flag predicates_init_flag = ATOMIC_FLAG_TEST_INIT;
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
static int coord_to_drop_from_plane(const s_point plane[3])
{
    s_point AB = subtract_points(plane[1], plane[0]);
    s_point AC = subtract_points(plane[2], plane[0]);
    s_point n  = cross_prod(AB, AC);

    double ax = fabs(n.coords[0]);
    double ay = fabs(n.coords[1]);
    double az = fabs(n.coords[2]);

    if (ax >= ay && ax >= az) return 0;  // drop X
    else if (ay >= ax && ay >= az) return 1;  // drop Y
    else return 2;  // drop Z
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

static double point_segment_dist2_2D(const double p[2], const double a[2], const double b[2], double EPS_degenerate, double *t_out)
{
    double vx=b[0]-a[0], vy=b[1]-a[1];
    double wx=p[0]-a[0], wy=p[1]-a[1];

    double denom = vx*vx + vy*vy;
    if (fabs(denom) < EPS_degenerate) {
        *t_out = 0.0;  return (p[0]-a[0])*(p[0]-a[0]) + (p[1]-a[1])*(p[1]-a[1]);
    } 
    double t = (wx*vx + wy*vy) / denom;
    *t_out = t;

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
    if ((a[0]-b[0])*(a[0]-b[0]) + (a[1]-b[1])*(a[1]-b[1]) <= TOL*TOL) return 1;
    else return 0;
}

static int points_equal_2D(const double a[2], const double b[2])
{
    if (a[0] == b[0] && a[1] == b[1]) return 1;
    else return 0;
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


/* Geometrical tests */
e_geom_test test_point_in_interval_1D(double x, double a, double b, double EPS_degenerate, double TOL)
{
    if (fabs(b - a) < EPS_degenerate) return TEST_ERROR;
    if (fabs(x - a) <= TOL|| fabs(x - b) <= TOL) return TEST_BOUNDARY;
    if (a > b) { double tmp = a; a = b; b = tmp; }
    if (x + TOL>= a && x - TOL<= b) return TEST_IN;
    return TEST_OUT;
}


static e_geom_test test_point_in_triangle_2D_robust_from_orientations(int o1, int o2, int o3)
{
    /* If any non-zero orientation disagrees, p is outside */
    int ref = (o1 != 0) ? o1 : ((o2 != 0) ? o2 : o3);
    if (ref == 0) return TEST_ERROR;
    if ((o1 != 0 && o1 != ref) || (o2 != 0 && o2 != ref) || (o3 != 0 && o3 != ref)) return TEST_OUT;
    if (o1 == 0 || o2 == 0 || o3 == 0) return TEST_BOUNDARY;
    return TEST_IN;
}


e_geom_test test_point_in_triangle_2D(const double a[2], const double b[2], const double c[2], const double p[2], double EPS_degenerate, double TOL)
{
    if (area_triangle_2D(a, b, c) < EPS_degenerate) return TEST_ERROR;

    /* Robust orientation tests of p vs each edge */
    int o1 = orient2d(a, b, p);
    int o2 = orient2d(b, c, p);
    int o3 = orient2d(c, a, p);

    if (TOL == 0) return test_point_in_triangle_2D_robust_from_orientations(o1, o2, o3);

    /* Numerical tolerance branch */
    double TOL2 = TOL*TOL;
    double d2, t;
    /* Check distance to edges */
    d2 = point_segment_dist2_2D(p, a, b, EPS_degenerate, &t);
    if (d2 <= TOL2 && t >= -EPS_degenerate && t <= 1.0 + EPS_degenerate) return TEST_BOUNDARY;

    d2 = point_segment_dist2_2D(p, b, c, EPS_degenerate, &t);
    if (d2 <= TOL2 && t >= -EPS_degenerate && t <= 1.0 + EPS_degenerate) return TEST_BOUNDARY;

    d2 = point_segment_dist2_2D(p, c, a, EPS_degenerate, &t);
    if (d2 <= TOL2 && t >= -EPS_degenerate && t <= 1.0 + EPS_degenerate) return TEST_BOUNDARY;

    /* If not on boundary, check robust orientation to determine TEST_IN/TEST_OUT */
    return test_point_in_triangle_2D_robust_from_orientations(o1, o2, o3);
}


static e_geom_test test_point_in_triangle_3D_robust(const s_point triangle[3], s_point p)
{
    if (orientation_robust(triangle, p) != 0) return TEST_OUT;

    /* Point coplanr. Project to 2D and check */   
    int drop = coord_to_drop_from_plane(triangle);
    double A2[2]; drop_to_2D(triangle[0], drop, A2);
    double B2[2]; drop_to_2D(triangle[1], drop, B2);
    double C2[2]; drop_to_2D(triangle[2], drop, C2);
    double p2[2]; drop_to_2D(p, drop, p2);

    return test_point_in_triangle_2D(A2, B2, C2, p2, 0, 0);
}


e_geom_test test_point_in_triangle_3D(const s_point triangle[3], s_point p, double EPS_degenerate, double TOL) 
{
    if (area_triangle(triangle) < EPS_degenerate) return TEST_ERROR;

    if (TOL == 0) return test_point_in_triangle_3D_robust(triangle, p);

    s_point closest = project_point_to_plane(p, triangle, EPS_degenerate);
    if (!point_is_valid(closest)) return TEST_ERROR;
    if (distance_squared(closest, p) > TOL*TOL) return TEST_OUT;
    
    /* Point is TOL close to plane. Project to 2D and check */
    int drop = coord_to_drop_from_plane(triangle);
    double A2[2]; drop_to_2D(triangle[0], drop, A2);
    double B2[2]; drop_to_2D(triangle[1], drop, B2);
    double C2[2]; drop_to_2D(triangle[2], drop, C2);
    double p2[2]; drop_to_2D(p, drop, p2);

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

    if (e0 == 0 || e1 == 0 || e2 == 0 || e3 == 0) return TEST_ERROR;

    if ((s0 != 0 && s0 * e0 < 0) ||
        (s1 != 0 && s1 * e1 < 0) ||
        (s2 != 0 && s2 * e2 < 0) ||
        (s3 != 0 && s3 * e3 < 0))
        return TEST_OUT;
    
    /* Count how many face-orientation_robusts are exactly zero */
    int zeros = (s0 == 0) + (s1 == 0) + (s2 == 0) + (s3 == 0);
    if (zeros == 0) return TEST_IN;   // strictly inside
    return TEST_BOUNDARY;
}


e_geom_test test_point_in_tetrahedron(const s_point tetra[4], s_point query, double EPS_degenerate, double TOL)
{   
    if (signed_volume_tetra(tetra) < EPS_degenerate) return TEST_ERROR;

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


e_intersect_type test_segment_segment_intersect_2D_robust(const double A1[2], const double A2[2], const double B1[2], const double B2[2])
{
    /* 1) Check case where segment is point */
    int A_is_point = points_equal_2D(A1, A2);
    int B_is_point = points_equal_2D(B1, B2);
    if (A_is_point && B_is_point)
        return points_equal_2D(A1,B1) ? INTERSECT_DEGENERATE : INTERSECT_EMPTY;
    if (A_is_point) {  /* check if A lies on segment B using orient2d and interval test */
        if (orient2d(B1, B2, A1) == 0 &&
            test_point_in_interval_1D(A1[0], B1[0], B2[0], 0, 0) != TEST_OUT &&
            test_point_in_interval_1D(A1[1], B1[1], B2[1], 0, 0) != TEST_OUT)
            return INTERSECT_DEGENERATE;
        return INTERSECT_EMPTY;
    }
    if (B_is_point) {
        if (orient2d(A1, A2, B1) == 0 &&
            test_point_in_interval_1D(B1[0], A1[0], A2[0], 0, 0) != TEST_OUT &&
            test_point_in_interval_1D(B1[1], A1[1], A2[1], 0, 0) != TEST_OUT)
            return INTERSECT_DEGENERATE;
        return INTERSECT_EMPTY;
    }


    /* 2) Segments are non-degenerate */
    double d1 = orient2d(A1, A2, B1);
    double d2 = orient2d(A1, A2, B2);
    double d3 = orient2d(B1, B2, A1);
    double d4 = orient2d(B1, B2, A2);

    /* 2a) Check if both collinear */
    if (d1 == 0 && d2 == 0 && d3 == 0 && d4 == 0) {
        double spanX = fabs(A1[0]-A2[0]) + fabs(B1[0]-B2[0]);
        double spanY = fabs(A1[1]-A2[1]) + fabs(B1[1]-B2[1]);
        int axis = (spanX >= spanY) ? 0 : 1;  /* choose primary axis */

        double a_min = fmin(A1[axis], A2[axis]), a_max = fmax(A1[axis], A2[axis]);
        double b_min = fmin(B1[axis], B2[axis]), b_max = fmax(B1[axis], B2[axis]);
        if (fmax(a_min,b_min) <= fmin(a_max,b_max)) return INTERSECT_DEGENERATE;
        else return INTERSECT_EMPTY;
    }
    
    /* 2b) Check proper intersection via straddling. */
     if ( ((d1>0 && d2<0) || (d1<0 && d2>0)) &&
          ((d3>0 && d4<0) || (d3<0 && d4>0)) )
        return INTERSECT_NONDEGENERATE;
    
    /* 2c) Single end colinear with segment. Check if point lies within the other segment */
    if (d1 == 0 &&
        test_point_in_interval_1D(B1[0], A1[0], A2[0], 0, 0) != TEST_OUT &&
        test_point_in_interval_1D(B1[1], A1[1], A2[1], 0, 0) != TEST_OUT)
        return INTERSECT_DEGENERATE;
    if (d2 == 0 &&
        test_point_in_interval_1D(B2[0], A1[0], A2[0], 0, 0) != TEST_OUT &&
        test_point_in_interval_1D(B2[1], A1[1], A2[1], 0, 0) != TEST_OUT)
        return INTERSECT_DEGENERATE;
    if (d3 == 0 &&
        test_point_in_interval_1D(A1[0], B1[0], B2[0], 0, 0) != TEST_OUT &&
        test_point_in_interval_1D(A1[1], B1[1], B2[1], 0, 0) != TEST_OUT)
        return INTERSECT_DEGENERATE;
    if (d4 == 0 &&
        test_point_in_interval_1D(A2[0], B1[0], B2[0], 0, 0) != TEST_OUT &&
        test_point_in_interval_1D(A2[1], B1[1], B2[1], 0, 0) != TEST_OUT)
        return INTERSECT_DEGENERATE;
    return INTERSECT_EMPTY;
}


e_intersect_type test_segment_triangle_intersect_2D_robust(const double S1[2], const double S2[2], const double A[2], const double B[2], const double C[2])
{
    e_geom_test i1 = test_point_in_triangle_2D(A, B, C, S1, 0, 0);
    e_geom_test i2 = test_point_in_triangle_2D(A, B, C, S2, 0, 0);
    if (i1 == TEST_ERROR || i2 == TEST_ERROR) return INTERSECT_ERROR;
    if (i1 == TEST_IN || i2 == TEST_IN) return INTERSECT_NONDEGENERATE;
    if (i1 == TEST_BOUNDARY && i2 == TEST_BOUNDARY) return INTERSECT_DEGENERATE;  /* Both on boundary of triangle */

    /* Both points are completely outside. Check intersections with triangle edges */
    e_intersect_type t1 = test_segment_segment_intersect_2D_robust(S1, S2, A, B);
    e_intersect_type t2 = test_segment_segment_intersect_2D_robust(S1, S2, B, C);
    e_intersect_type t3 = test_segment_segment_intersect_2D_robust(S1, S2, C, A);

    if (t1 == INTERSECT_NONDEGENERATE || t2 == INTERSECT_NONDEGENERATE || t3 == INTERSECT_NONDEGENERATE) return INTERSECT_NONDEGENERATE;
    if (t1 == INTERSECT_DEGENERATE || t2 == INTERSECT_DEGENERATE || t3 == INTERSECT_DEGENERATE) return INTERSECT_DEGENERATE;
    return INTERSECT_EMPTY;
}


e_intersect_type test_segment_triangle_intersect_3D_robust(const s_point segment[2], const s_point triangle[3])
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
        
        e_intersect_type type2D = test_segment_triangle_intersect_2D_robust(P2, D2, A2, B2, C2);
        if (type2D == INTERSECT_ERROR) return INTERSECT_ERROR;
        if (type2D == INTERSECT_NONDEGENERATE || type2D == INTERSECT_DEGENERATE) return INTERSECT_DEGENERATE;
        return INTERSECT_EMPTY; 
    }
    else if (o1 == 0) {
        e_geom_test test = test_point_in_triangle_3D_robust(triangle, segment[0]);
        if (test == TEST_ERROR) return INTERSECT_ERROR;
        if (test == TEST_IN || test == TEST_BOUNDARY) return INTERSECT_DEGENERATE;
        return INTERSECT_EMPTY;
    }
    else if (o2 == 0) {
        e_geom_test test = test_point_in_triangle_3D_robust(triangle, segment[1]);
        if (test == TEST_ERROR) return INTERSECT_ERROR;
        if (test == TEST_IN || test == TEST_BOUNDARY) return INTERSECT_DEGENERATE;
        return INTERSECT_EMPTY;
    }

    /* Check that segment crosses triangle's plane */
    if (o1 == o2) return INTERSECT_EMPTY;

    /* Segura's algorithm */
    s_point v1 = triangle[0], v2 = triangle[1], v3 = triangle[2];

    int index_q1 = orientation_robust(triangle, segment[0]) == 1 ? 0 : 1;
    int index_q2 = index_q1 == 0 ? 1 : 0;
    s_point q1 = segment[index_q1];
    s_point q2 = segment[index_q2];

    s_point aux[3];
    aux[0] = q2; aux[1] = v1; aux[2] = v2;
    int s1 = orientation_robust(aux, q1);

    aux[0] = q2; aux[1] = v3; aux[2] = v2;
    int s2 = orientation_robust(aux, q1);

    aux[0] = q2; aux[1] = v1; aux[2] = v3;
    int s3 = orientation_robust(aux, q1);

    if (s1 < 0 || s2 < 0 || s3 < 0) return INTERSECT_EMPTY;
    /* All >= 0 */
    if (s1 == 0 || s2 == 0 || s3 == 0) return INTERSECT_DEGENERATE;  
    return INTERSECT_NONDEGENERATE;  
}



/* TEST_INTERSECTIONS */
/* Even though the count should be robust, the point of intersection need not be robust... */
int segment_segment_intersection_2D(const double A1[2], const double A2[2], const double B1[2], const double B2[2], double EPS_degenerate, double TOL, double out[4])
{   /* RETURNS: abs(output) = number of intersections
       If positive, we could compute the intersection respecting EPS_degenerate */
    /* 1) Degeneracy handling */
    int A_is_point = points_close_2D(A1, A2, EPS_degenerate);
    int B_is_point = points_close_2D(B1, B2, EPS_degenerate);
    if (A_is_point && B_is_point) {
        if (fabs(A1[0]-B1[0])<=TOL && fabs(A1[1]-B1[1])<=TOL) {  /* Both are same point */
            if (out) { out[0] = A1[0]; out[1] = A1[1]; } return 1;
        } else return 0;
    } else if (A_is_point) {
        if (orient2d(B1, B2, A1) == 0 &&
            test_point_in_interval_1D(A1[0], B1[0], B2[0], EPS_degenerate, TOL) != TEST_OUT &&
            test_point_in_interval_1D(A1[1], B1[1], B2[1], EPS_degenerate, TOL) != TEST_OUT) {
            if (out) { out[0] = A1[0]; out[1] = A1[1]; } return 1;
        } else return 0;
    } else if (B_is_point) {
        if (orient2d(A1, A2, B1) == 0 &&
            test_point_in_interval_1D(B1[0], A1[0], A2[0], EPS_degenerate, TOL) != TEST_OUT &&
            test_point_in_interval_1D(B1[1], A1[1], A2[1], EPS_degenerate, TOL) != TEST_OUT) {
            if (out) { out[0] = B1[0]; out[1] = B1[1]; } return 1;
        } else return 0;
    }

    /* 2) General (non-degenerate) case */
    double d1 = orient2d(A1, A2, B1);
    double d2 = orient2d(A1, A2, B2);
    double d3 = orient2d(B1, B2, A1);
    double d4 = orient2d(B1, B2, A2);

    /* 2a) Check if collinear. Works for robust and non-robust cases */
    if (fabs(d1) <= EPS_degenerate && fabs(d2) <= EPS_degenerate &&
        fabs(d3) <= EPS_degenerate && fabs(d4) <= EPS_degenerate) {
        /* collinear: overlapping test via projections */
        double min_x = fmax(fmin(A1[0], A2[0]), fmin(B1[0], B2[0]));
        double max_x = fmin(fmax(A1[0], A2[0]), fmax(B1[0], B2[0]));
        double min_y = fmax(fmin(A1[1], A2[1]), fmin(B1[1], B2[1]));
        double max_y = fmin(fmax(A1[1], A2[1]), fmax(B1[1], B2[1]));
        if (min_x <= max_x + TOL && min_y <= max_y + TOL) {
            if (out) {
                out[0] = min_x; out[1] = min_y;
                out[2] = max_x; out[3] = max_y;
            }
            return 2; /* overlapping collinear segments */
        }
        return 0;
    }
    
    /* 2b) Check proper intersection via straddling. Works for robust and non-robust cases. */
     if ( ((d1 > EPS_degenerate && d2 < -EPS_degenerate) || (d1 < -EPS_degenerate && d2 > EPS_degenerate)) &&
          ((d3 > EPS_degenerate && d4 < -EPS_degenerate) || (d3 < -EPS_degenerate && d4 > EPS_degenerate)) ) {
        if (out) {
            double A_a = A2[1] - A1[1];
            double A_b = A1[0] - A2[0];
            double A_c = A_a*A1[0] + A_b*A1[1];
            double B_a = B2[1] - B1[1];
            double B_b = B1[0] - B2[0];
            double B_c = B_a*B1[0] + B_b*B1[1];

            double det = A_a*B_b - B_a*A_b;
            if (fabs(det) <= EPS_degenerate) return -1;  /* nearly parallel */
            out[0] = (B_b*A_c - A_b*B_c)/det;
            out[1] = (A_a*B_c - B_a*A_c)/det;
        }
        return 1;
    }
    
    /* 2c) Handle endpoint touching */
    if (TOL > 0) {  /* If TOL > 0, compare endpoints across segments (A vs B)*/
        const double endpointsA[2][2] = { {A1[0], A1[1]}, {A2[0], A2[1]} };
        const double endpointsB[2][2] = { {B1[0], B1[1]}, {B2[0], B2[1]} };
        for (int i = 0; i < 2; ++i) {
            for (int j = 0; j < 2; ++j) {
                if (points_close_2D(endpointsA[i], endpointsB[j], TOL)) {
                    if (out) { out[0] = endpointsA[i][0]; out[1] = endpointsA[i][1]; }
                    return 1;
                }
            }
        }
    } else {  /* TOL == 0: check if point lies within the other segment */
        if (d1 == 0 &&
            test_point_in_interval_1D(B1[0], A1[0], A2[0], 0, 0) != TEST_OUT &&
            test_point_in_interval_1D(B1[1], A1[1], A2[1], 0, 0) != TEST_OUT) {
            if (out) { out[0] = B1[0]; out[1] = B1[1]; } return 1;
        }
        if (d2 == 0 &&
            test_point_in_interval_1D(B2[0], A1[0], A2[0], 0, 0) != TEST_OUT &&
            test_point_in_interval_1D(B2[1], A1[1], A2[1], 0, 0) != TEST_OUT) {
            if (out) { out[0] = B2[0]; out[1] = B2[1]; } return 1;
        }
        if (d3 == 0 &&
            test_point_in_interval_1D(A1[0], B1[0], B2[0], 0, 0) != TEST_OUT &&
            test_point_in_interval_1D(A1[1], B1[1], B2[1], 0, 0) != TEST_OUT) {
            if (out) { out[0] = A1[0]; out[1] = A1[1]; } return 1;
        }
        if (d4 == 0 &&
            test_point_in_interval_1D(A2[0], B1[0], B2[0], 0, 0) != TEST_OUT &&
            test_point_in_interval_1D(A2[1], B1[1], B2[1], 0, 0) != TEST_OUT) {
            if (out) { out[0] = A2[0]; out[1] = A2[1]; } return 1;
        }
    }

    return 0;
}


static int segment_plane_intersection_robust(const s_point seg[2], const s_point plane[3], double EPS_degenerate, s_point out[2])
{
    int o0 = orientation_robust(plane, seg[0]);
    int o1 = orientation_robust(plane, seg[1]);

    if (o0 == 0 && o1 == 0) { out[0] = seg[0]; out[1] = seg[1]; return 2; } 
    else if (o0 == 0) { out[0] = seg[0]; return 1; } 
    else if (o1 == 0) { out[0] = seg[1]; return 1; } 
    else if (o0 == o1) { return 0; }  /* same side of plane */
    else { 
        /* fallback to numeric interpolation for exact intersection */
        s_point n = cross_prod(subtract_points(plane[1], plane[0]),
                               subtract_points(plane[2], plane[0]));
        n = normalize_vec(n, EPS_degenerate);
        if (!point_is_valid(n)) return -1;
        double d = dot_prod(n, plane[0]);
        double sA = dot_prod(seg[0], n) - d;
        double sB = dot_prod(seg[1], n) - d;
        if (fabs(sA - sB) < EPS_degenerate) return -1;
        double t = sA / (sA - sB);
        out[0] = interpolate_points(seg[0], seg[1], t);
        return 1;
    }
}


int segment_plane_intersection(const s_point seg[2], const s_point plane[3], double EPS_degenerate, double TOL, s_point out[2])
{
    if (TOL == 0) return segment_plane_intersection_robust(seg, plane, EPS_degenerate, out);

    s_point n = cross_prod(subtract_points(plane[1], plane[0]),
                           subtract_points(plane[2], plane[0]));
    n = normalize_vec(n, EPS_degenerate);
    if (!point_is_valid(n)) return 0;
    double d = dot_prod(n, plane[0]);
    double sA = dot_prod(seg[0], n) - d;
    double sB = dot_prod(seg[1], n) - d;

    /* both endpoints on same side (and not near plane) */
    if (sA * sB > TOL*TOL) return 0;

    /* both nearly on plane : whole segment lies in it */
    if (fabs(sA) < TOL && fabs(sB) < TOL) {
        out[0] = project_point_to_plane(seg[0], plane, EPS_degenerate);
        out[1] = project_point_to_plane(seg[1], plane, EPS_degenerate);
        if (!point_is_valid(out[0]) || !point_is_valid(out[1])) return -2;
        return 2;
    }

    /* find parametric intersection */
    double denom = sA - sB;
    if (fabs(denom) < EPS_degenerate) return -1;

    double t = sA / denom;
    if (t < -TOL || t > 1 + TOL) return 0;
    if (t < 0) t = 0;
    if (t > 1) t = 1;
    out[0] = interpolate_points(seg[0], seg[1], t);
    return 1;
}


int segment_triangle_intersection_2D(const double s1[2], const double s2[2], const double a[2], const double b[2], const double c[2], double EPS_degenerate, double TOL, double out[4])
{
    /* Check if inside */
    e_geom_test i1 = test_point_in_triangle_2D(a, b, c, s1, EPS_degenerate, TOL);
    e_geom_test i2 = test_point_in_triangle_2D(a, b, c, s2, EPS_degenerate, TOL);
    if (i1 == TEST_ERROR || i2 == TEST_ERROR) return 0;
    if (i1 != TEST_OUT && i2 != TEST_OUT) {  /* Both inside / on boundary */
        if (out) { out[0] = s1[0]; out[1] = s1[1]; out[2] = s2[0]; out[3] = s2[1]; } return 2;
    } else if (i1 == TEST_BOUNDARY && i2 == TEST_OUT) { 
        if (out) { out[0] = s1[0]; out[1] = s1[1]; } return 1; 
    } else if (i2 == TEST_BOUNDARY && i1 == TEST_OUT) { 
        if (out) { out[0] = s2[0]; out[1] = s2[1]; } return 1;
    }

    /* Check segment_segment intersections in 2D. Already considers robust / non-robust cases! */
    double intersection1[4], intersection2[4], intersection3[4];
    int n1 = segment_segment_intersection_2D(s1, s2, a, b, EPS_degenerate, TOL, intersection1);
    int n2 = segment_segment_intersection_2D(s1, s2, b, c, EPS_degenerate, TOL, intersection2);
    int n3 = segment_segment_intersection_2D(s1, s2, c, a, EPS_degenerate, TOL, intersection3);
    if (!out && (n1>0 || n2>0 || n3>0)) return (int)fmax(n1, fmax(n2, n3)); 

    int h = 0;
    double TOL2 = TOL*TOL;
    /* append with dedupe */
    for (int ii = 0; ii < n1; ++ii) {
        double p[2] = { intersection1[ii*2+0], intersection1[ii*2+1] };
        append_unique_2D_lim2(out, &h, p, TOL2);
    }
    for (int ii = 0; ii < n2 && h < 2; ++ii) {
        double p[2] = { intersection2[ii*2+0], intersection2[ii*2+1] };
        append_unique_2D_lim2(out, &h, p, TOL2);
    }
    for (int ii = 0; ii < n3 && h < 2; ++ii) {
        double p[2] = { intersection3[ii*2+0], intersection3[ii*2+1] };
        append_unique_2D_lim2(out, &h, p, TOL2);
    }
    if (h == 2) return 2;  /* Two total intersections */

    /* If one endpoint is inside and we found at least one intersection, include it */
    if (((i1 != TEST_OUT &&  i2 == TEST_OUT) || (i2 != TEST_OUT &&  i1 == TEST_OUT)) && h == 1) {
        const double *inside = (i1 != TEST_OUT) ? s1 : s2;
        append_unique_2D_lim2(out, &h, inside, TOL2);
        return 2;
    }

    if (h == 1) return 1;       
    else return 0;
}


static int segment_triangle_intersection_3D_robust(const s_point segment[2], const s_point triangle[3], s_point out[2])
{
    int o1 = orientation_robust(triangle, segment[0]);
    int o2 = orientation_robust(triangle, segment[1]);
    if (o1 != 0 && o2 != 0 && o1 == o2) return 0;  /* Both on the same side of plane */

    int drop = coord_to_drop_from_plane(triangle);
    double A2[2]; drop_to_2D(triangle[0], drop, A2);
    double B2[2]; drop_to_2D(triangle[1], drop, B2);
    double C2[2]; drop_to_2D(triangle[2], drop, C2);
    double P2[2]; drop_to_2D(segment[0], drop, P2);
    double D2[2]; drop_to_2D(segment[1], drop, D2);

    double intersections2D[4];
    int Ni = segment_triangle_intersection_2D(P2, D2, A2, B2, C2, 0, 0, intersections2D);
    if (out) {
        for (int ii=0; ii<Ni; ii++) {
            out[ii] = lift_point_from_dropped_2D(triangle, drop, &intersections2D[ii*2], 0);
        }
    }
    return Ni;
}


int segment_triangle_intersection_3D(const s_point segment[2], const s_point triangle[3], double EPS_degenerate, double TOL, s_point out[2])
{
    if (area_triangle(triangle) < EPS_degenerate) return 0;

    if (TOL == 0) return segment_triangle_intersection_3D_robust(segment, triangle, out);

    /* Non-robust branch */
    double TOL2 = TOL*TOL;
    s_point tmp[2];
    int Ni_plane = segment_plane_intersection(segment, triangle, EPS_degenerate, TOL, tmp);
    int count = 0;
    /* Intersection(s) close enought to closest point on triangle? */
    for (int ii=0; ii<Ni_plane; ii++) {
        s_point closest = closest_point_on_triangle(triangle, EPS_degenerate, tmp[ii]);
        if (distance_squared(closest, tmp[ii]) < TOL2) count++;
    }
    if (out) { out[0] = tmp[0]; out[1] = tmp[1]; }
    return count;
}


