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


// int segment_crosses_triangle_robust_3D(const s_point triangle[3], s_point a, s_point b)
// {
//     if (orientation_robust(triangle, a) == orientation_robust(triangle, b)) return 0;
//     s_point aux[3];
//     aux[2] = a; 
//
//     aux[0] = triangle[0];   aux[1] = triangle[1];
//     int s1 = orientation_robust(aux, b);
//
//     aux[0] = triangle[1];   aux[1] = triangle[2];
//     int s2 = orientation_robust(aux, b);
//
//     aux[0] = triangle[2];   aux[1] = triangle[0];
//     int s3 = orientation_robust(aux, b);
//     
//     if (s1 == s2 && s2 == s3 && s3 == s1) return 1;
//     else return 0;
// }


e_geom_test test_point_in_interval_1D(double x, double a, double b, double EPS_degenerate, double TOL)
{
    if (fabs(b - a) < EPS_degenerate) return ERROR;
    if (fabs(x - a) <= TOL|| fabs(x - b) <= TOL) return BOUNDARY;
    if (a > b) { double tmp = a; a = b; b = tmp; }
    if (x + TOL>= a && x - TOL<= b) return IN;
    return OUT;
}


int segments_candegenerate_intersect_2D(const double A1[2], const double A2[2], const double B1[2], const double B2[2], double EPS)
{
    /* 1) Degenerate B (point) -> treat as point-segment test */
    if (fabs(B1[0]-B2[0])<EPS && fabs(B1[1]-B2[1])<EPS) {
        /* If A is also degenerate (a point) -> check coincidence */
        if (fabs(A1[0]-A2[0])<EPS && fabs(A1[1]-A2[1])<EPS)
            return (fabs(B1[0]-A1[0])<EPS && fabs(B1[1]-A1[1])<EPS) ? 1 : 0;
        
        /* X must be collinear with A1-A2 and within bounding-box of A */
        if (orient2d(A1, A2, B1) == 0 &&
            test_point_in_interval_1D(B1[0], A1[0], A2[0], EPS, EPS) != OUT &&
            test_point_in_interval_1D(B1[1], A1[1], A2[1], EPS, EPS) != OUT) 
            return 1;
        else return 0;
    }

    /* 2) Degenerate A (point) -> treat as point-segment test */
    if (fabs(A1[0]-A2[0])<EPS && fabs(A1[1]-A2[1])<EPS) {
        if (orient2d(B1, B2, A1) == 0 &&
            test_point_in_interval_1D(A1[0], B1[0], B2[0], EPS, EPS) != OUT &&
            test_point_in_interval_1D(A1[1], B1[1], B2[1], EPS, EPS) != OUT)
            return 1;
        else return 0;
    }

    /* 3) General case: use two straddling tests  */
    /* 3a) [A1,A2] vs {B1,B2} */
    int o1 = orient2d(A1, A2, B1);
    int o2 = orient2d(A1, A2, B2);
    if (o1 != 0 && o2 != 0 && o1 == o2) return 0; /* both same nonzero side -> no intersection */
    /* 3b) [B1,B2] vs {A1,A2} */
    o1 = orient2d(B1, B2, A1);
    o2 = orient2d(B1, B2, A2);
    if (o1 != 0 && o2 != 0 && o1 == o2) return 0;

    /* If we get here, segments (possibly including collinear/end-touching cases) intersect */
    return 1;
}


static double point_segment_dist2_2D(const double p[2], const double a[2], const double b[2], double EPS_degenerate, double *t_out)
{
    double vx = b[0] - a[0];
    double vy = b[1] - a[1];
    double wx = p[0] - a[0];
    double wy = p[1] - a[1];

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


e_geom_test test_point_in_triangle_2D(const double a[2], const double b[2], const double c[2], const double p[2], double EPS_degenerate, double TOL)
{
    /* Degenerate triangle check */
    if (area_triangle_2D(a, b, c) < EPS_degenerate) return ERROR;

    /* Robust orientation tests of p vs each edge */
    int o1 = orient2d(a, b, p);
    int o2 = orient2d(b, c, p);
    int o3 = orient2d(c, a, p);

    /* Robust branch */
    if (TOL== 0) {
        /* If any non-zero orientation disagrees, p is outside */
        int ref = (o1 != 0) ? o1 : ((o2 != 0) ? o2 : o3);
        if (ref == 0) return ERROR;
        if ((o1 != 0 && o1 != ref) || (o2 != 0 && o2 != ref) || (o3 != 0 && o3 != ref)) return OUT;
        if (o1 == 0 || o2 == 0 || o3 == 0) return BOUNDARY;
        return IN;
    }

    /* Numerical tolerance branch */
    double TOL2 = TOL*TOL;
    double d2, t;
    /* Check distance to edges */
    d2 = point_segment_dist2_2D(p, a, b, EPS_degenerate, &t);
    if (d2 <= TOL2 && t >= -EPS_degenerate && t <= 1.0 + EPS_degenerate) return BOUNDARY;

    d2 = point_segment_dist2_2D(p, b, c, EPS_degenerate, &t);
    if (d2 <= TOL2 && t >= -EPS_degenerate && t <= 1.0 + EPS_degenerate) return BOUNDARY;

    d2 = point_segment_dist2_2D(p, c, a, EPS_degenerate, &t);
    if (d2 <= TOL2 && t >= -EPS_degenerate && t <= 1.0 + EPS_degenerate) return BOUNDARY;

    /* If not on boundary, check robust orientation to determine IN/OUT */
    int ref = (o1 != 0) ? o1 : ((o2 != 0) ? o2 : o3);
    if (ref == 0) return ERROR;
    if ((o1 != 0 && o1 != ref) || (o2 != 0 && o2 != ref) || (o3 != 0 && o3 != ref)) return OUT;

    return IN;
}


e_geom_test test_point_in_triangle_3D(const s_point triangle[3], s_point p, double EPS_degenerate, double TOL) 
{
    if (area_triangle(triangle) < EPS_degenerate) return ERROR;

    s_point n, t1, t2;
    if (!basis_vectors_plane(triangle, EPS_degenerate, &n, &t1, &t2)) return ERROR;

    /* Check if point is coplanar or EPS-coplanar before proceeding */
    if (TOL == 0) {
        if (orientation_robust(triangle, p) != 0) return OUT;
    } else {
        s_point closest = closest_point_on_triangle(triangle, EPS_degenerate, p);
        if (!point_is_valid(closest)) return ERROR;
        if (distance_squared(closest, p) > TOL*TOL) return OUT;
    }

    s_point p_proj = project_point_to_plane(p, triangle, EPS_degenerate);
    if (!point_is_valid(p_proj)) return ERROR;
    double v1[2] = { dot_prod(triangle[0], t1), dot_prod(triangle[0], t2) };
    double v2[2] = { dot_prod(triangle[1], t1), dot_prod(triangle[1], t2) };
    double v3[2] = { dot_prod(triangle[2], t1), dot_prod(triangle[2], t2) };
    double p2d[2] = { dot_prod(p_proj, t1), dot_prod(p_proj, t2) };

    e_geom_test res = test_point_in_triangle_2D(v1, v2, v3, p2d, EPS_degenerate, TOL);
    return res;
}


e_geom_test test_point_in_tetrahedron(const s_point tetra[4], s_point query)
{   // TODO add non-robust branch?
    s_point tmp[3];

    /* First compute reference signs*/
    tmp[0] = tetra[1];   tmp[1] = tetra[2];   tmp[2] = tetra[3];
    int e0 = orientation_robust(tmp, tetra[0]);

    tmp[0] = tetra[0];   tmp[1] = tetra[3];   tmp[2] = tetra[2];
    int e1 = orientation_robust(tmp, tetra[1]);

    tmp[0] = tetra[0];   tmp[1] = tetra[1];   tmp[2] = tetra[3];
    int e2 = orientation_robust(tmp, tetra[2]);

    tmp[0] = tetra[0];   tmp[1] = tetra[2];   tmp[2] = tetra[1];
    int e3 = orientation_robust(tmp, tetra[3]);

    if (e0 == 0 || e1 == 0 || e2 == 0 || e3 == 0) return ERROR;

    /* Compute signs for the query */
    tmp[0] = tetra[1];   tmp[1] = tetra[2];   tmp[2] = tetra[3];
    int s0 = orientation_robust(tmp, query);

    tmp[0] = tetra[0];   tmp[1] = tetra[3];   tmp[2] = tetra[2];
    int s1 = orientation_robust(tmp, query);

    tmp[0] = tetra[0];   tmp[1] = tetra[1];   tmp[2] = tetra[3];
    int s2 = orientation_robust(tmp, query);

    tmp[0] = tetra[0];   tmp[1] = tetra[2];   tmp[2] = tetra[1];
    int s3 = orientation_robust(tmp, query);

    if ((s0 != 0 && s0 * e0 < 0) ||
        (s1 != 0 && s1 * e1 < 0) ||
        (s2 != 0 && s2 * e2 < 0) ||
        (s3 != 0 && s3 * e3 < 0))
        return OUT;
    
    /* Count how many face-orientation_robusts are exactly zero */
    int zeros = (s0 == 0) + (s1 == 0) + (s2 == 0) + (s3 == 0);
    if (zeros == 0) return IN;   // strictly inside
    return BOUNDARY;
}


s_points_test test_points_in_halfspace(const s_point plane_ordered[3], const s_points *points, double EPS_degenerate, double TOL, e_geom_test out_buff[points->N])
{   /* Right hand rule: normal outwards */
    if (out_buff == NULL) out_buff = malloc(points->N * sizeof(e_geom_test));
    int Nin = 0, Nbdy = 0, Nout = 0;

    if (TOL == 0) {  /* Robust branch */
        for (int ii=0; ii<points->N; ii++) {
            int o = orientation_robust(plane_ordered, points->p[ii]);
            if (o == 1) { out_buff[ii] = IN; Nin++; } 
            else if (o == 0) { out_buff[ii] = BOUNDARY; Nbdy++; }
            else { out_buff[ii] = OUT; Nout++; }
        } 
    } else {  /* Non-robust branch */
        s_point n = cross_prod(subtract_points(plane_ordered[1], plane_ordered[0]),
                               subtract_points(plane_ordered[2], plane_ordered[0]));
        n = normalize_vec(n, EPS_degenerate);
        if (!point_is_valid(n)) {
            for (int ii=0; ii<points->N; ii++) out_buff[ii] = ERROR;
            return (s_points_test){ .Nin = 0, .Nbdy = 0, .Nout = 0, .Nerr = points->N, .indicator = out_buff };
        }
        double d = dot_prod(n, plane_ordered[0]);

        for (int ii=0; ii<points->N; ii++) {
            double s = dot_prod(points->p[ii], n) - d;
            if (s < TOL) { out_buff[ii] = IN; Nin++; } 
            else if (fabs(s) <= TOL) { out_buff[ii] = BOUNDARY; Nbdy++; } 
            else { out_buff[ii] = OUT; Nout++; }
        }
    }
    
    return (s_points_test){ .Nin = Nin, .Nbdy = Nbdy, .Nout = Nout, .Nerr = 0, .indicator = out_buff };
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
        if (!point_is_valid(n)) return 0;
        double d = dot_prod(n, plane[0]);
        double sA = dot_prod(seg[0], n) - d;
        double sB = dot_prod(seg[1], n) - d;
        if (fabs(sA - sB) < EPS_degenerate) return 0;
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
        out[0] = project_point_to_plane(seg[0], plane, EPS_degenerate);;
        out[1] = project_point_to_plane(seg[1], plane, EPS_degenerate);;
        return 2;
    }

    /* find parametric intersection */
    double denom = sA - sB;
    if (fabs(denom) < EPS_degenerate) return 0;

    double t = sA / denom;
    if (t < -TOL || t > 1 + TOL) return 0;

    out[0] = interpolate_points(seg[0], seg[1], t);
    return 1;
}


int segment_segment_intersection_2D(const double A1[2], const double A2[2], const double B1[2], const double B2[2], double EPS_degenerate, double TOL, double out[4])
{
    /* 1) Degeneracy handling */
    int A_is_point = (fabs(A1[0]-A2[0]) <= EPS_degenerate && fabs(A1[1]-A2[1]) <= EPS_degenerate);
    int B_is_point = (fabs(B1[0]-B2[0]) <= EPS_degenerate && fabs(B1[1]-B2[1]) <= EPS_degenerate);

    if (A_is_point && B_is_point) {
        if (fabs(A1[0]-B1[0]) <= TOL && fabs(A1[1]-B1[1]) <= TOL) {
            out[0] = A1[0]; out[1] = A1[1];
            return 1;
        }
        return 0;
    }
    if (A_is_point) {
        if (orient2d(B1, B2, A1) == 0 &&
            test_point_in_interval_1D(A1[0], B1[0], B2[0], EPS_degenerate, TOL) != OUT &&
            test_point_in_interval_1D(A1[1], B1[1], B2[1], EPS_degenerate, TOL) != OUT) {
            out[0] = A1[0]; out[1] = A1[1]; return 1;
        }
        return 0;
    }
    if (B_is_point) {
        if (orient2d(A1, A2, B1) == 0 &&
            test_point_in_interval_1D(B1[0], A1[0], A2[0], EPS_degenerate, TOL) != OUT &&
            test_point_in_interval_1D(B1[1], A1[1], A2[1], EPS_degenerate, TOL) != OUT) {
            out[0] = B1[0]; out[1] = B1[1]; return 1;
        }
        return 0;
    }

    /* 2) General (non-degenerate) case */
    double d1 = orient2d(A1, A2, B1);
    double d2 = orient2d(A1, A2, B2);
    double d3 = orient2d(B1, B2, A1);
    double d4 = orient2d(B1, B2, A2);

    /* 2a) Check if collinear */
    if (d1 == 0 && d2 == 0 && d3 == 0 && d4 == 0) {
        /* Project onto X or Y axis depending on largest span */
        double min_x = fmax(fmin(A1[0], A2[0]), fmin(B1[0], B2[0]));
        double max_x = fmin(fmax(A1[0], A2[0]), fmax(B1[0], B2[0]));
        double min_y = fmax(fmin(A1[1], A2[1]), fmin(B1[1], B2[1]));
        double max_y = fmin(fmax(A1[1], A2[1]), fmax(B1[1], B2[1]));
        if (min_x <= max_x + TOL && min_y <= max_y + TOL) {  /* overlapping or touching collinear segments */
            out[0] = min_x; out[1] = min_y;
            out[2] = max_x; out[3] = max_y;
            return 2;
        }
        return 0;
    }

    /* 2b) Check proper intersection via straddling */
    if (((d1 > 0 && d2 < 0) || (d1 < 0 && d2 > 0)) &&
        ((d3 > 0 && d4 < 0) || (d3 < 0 && d4 > 0))) {
        /* Compute intersection using line equations */
        double A_a = A2[1] - A1[1];
        double A_b = A1[0] - A2[0];
        double A_c = A_a*A1[0] + A_b*A1[1];

        double B_a = B2[1] - B1[1];
        double B_b = B1[0] - B2[0];
        double B_c = B_a*B1[0] + B_b*B1[1];

        double det = A_a*B_b - B_a*A_b;
        if (fabs(det) < EPS_degenerate) return 0;  /* nearly parallel */

        out[0] = (B_b*A_c - A_b*B_c)/det;
        out[1] = (A_a*B_c - B_a*A_c)/det;
        return 1;
    }

    if (TOL > 0) {  /* 2c) Handle endpoint touching within tolerance */
        double endpoints[4][2] = { {A1[0],A1[1]}, {A2[0],A2[1]}, {B1[0],B1[1]}, {B2[0],B2[1]} };
        for (int i=0; i<4; i++)
            for (int j=i+1; j<4; j++)
                if (fabs(endpoints[i][0] - endpoints[j][0]) < TOL && fabs(endpoints[i][1] - endpoints[j][1]) < TOL) {
                    out[0] = endpoints[i][0]; out[1] = endpoints[i][1]; return 1;
                }
    }

    return 0;
}


static int append_unique_2D_lim2(double out[4], int *h, const double p[2], double TOL2)
{
    for (int i = 0; i < *h; ++i) {
        double q[2] = { out[i*2+0], out[i*2+1] };
        if ((q[0]-p[0])*(q[0]-p[0]) + (q[1]-p[1])*(q[1]-p[1]) <= TOL2) return 0;
    }
    if (*h >= 2) return 0;  /* already full (safety) */
    out[(*h)*2+0] = p[0];
    out[(*h)*2+1] = p[1];
    (*h)++;
    return 1;
}


int clip_segment_with_triangle_2D(const double s1[2], const double s2[2], const double v1[2], const double v2[2], const double v3[2], double EPS_degenerate, double TOL, double out[4])
{
    /* Check if inside */
    e_geom_test i1 = test_point_in_triangle_2D(v1, v2, v3, s1, EPS_degenerate, TOL);
    e_geom_test i2 = test_point_in_triangle_2D(v1, v2, v3, s2, EPS_degenerate, TOL);
    if (i1 == ERROR || i2 == ERROR) return 0;
    if (i1 != OUT && i2 != OUT) {  /* Both inside / on boundary */
        out[0] = s1[0]; out[1] = s1[1];
        out[2] = s2[0]; out[3] = s2[1];
        return 2;
    }
    if (i1 == BOUNDARY && i2 == OUT) { out[0] = s1[0]; out[1] = s1[1]; return 1; }
    if (i2 == BOUNDARY && i1 == OUT) { out[0] = s2[0]; out[1] = s2[1]; return 1; }

    // Check segment_segment intersections in 2D
    double intersection1[4], intersection2[4], intersection3[4];
    int n1 = segment_segment_intersection_2D(s1, s2, v1, v2, EPS_degenerate, TOL, intersection1);
    int n2 = segment_segment_intersection_2D(s1, s2, v2, v3, EPS_degenerate, TOL, intersection2);
    int n3 = segment_segment_intersection_2D(s1, s2, v3, v1, EPS_degenerate, TOL, intersection3);

    int h = 0;
    double TOL2 = TOL*TOL;
    /* append with dedupe */
    for (int ii = 0; ii < n1; ++ii) {
        double p[2] = { intersection1[ii*2+0], intersection1[ii*2+1] };
        append_unique_2D_lim2(out, &h, p, TOL2);
        if (h == 2) break;
    }
    for (int ii = 0; ii < n2 && h < 2; ++ii) {
        double p[2] = { intersection2[ii*2+0], intersection2[ii*2+1] };
        append_unique_2D_lim2(out, &h, p, TOL2);
        if (h == 2) break;
    }
    for (int ii = 0; ii < n3 && h < 2; ++ii) {
        double p[2] = { intersection3[ii*2+0], intersection3[ii*2+1] };
        append_unique_2D_lim2(out, &h, p, TOL2);
        if (h == 2) break;
    }
    if (h == 2) return 2;  /* Two total intersections */

    /* If one endpoint is inside and we found at least one intersection, include it */
    if (((i1 != OUT &&  i2 == OUT) || (i2 != OUT &&  i1 == OUT)) && h == 1) {
        const double *inside = (i1 != OUT) ? s1 : s2;
        append_unique_2D_lim2(out, &h, inside, TOL2);
        return 2;
    }

    if (h == 1) return 1;       
    else return 0;
}

