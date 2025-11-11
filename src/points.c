#include "points.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>


int point_is_valid(s_point p)
{
    if (isnan(p.x) || isnan(p.y) || isnan(p.z)) return 0;
    return 1;
}


int points_is_valid(const s_points *p)
{
    if (p->p == NULL) return 0;
    else return 1;
}


void free_points(s_points *points)
{
    free(points->p);
    memset(points, 0, sizeof(s_points));
}


void print_points(const s_points *points)
{
    for (int ii=0; ii<points->N; ii++) {
        printf("%f, %f, %f\n", points->p[ii].x, points->p[ii].y, points->p[ii].z);
    }
}


static int count_lines(FILE* file)
{
    const int BUF_SIZE = 2048;
    char buf[BUF_SIZE];
    int counter = 0;
    while(1) {
        int nread = fread(buf, 1, BUF_SIZE, file);
        if (ferror(file)) return -1;

        for(int ii=0; ii<nread; ii++)
            if (buf[ii] == '\n') counter++;

        if (nread < BUF_SIZE && feof(file)) {
            if (nread > 0 && buf[nread - 1] != '\n') counter++;
            break;
        }
    }
    rewind(file);
    return counter;
}


s_points read_points_from_csv(const char *file)
{
    s_points out = {0};

    FILE *f = fopen(file, "r");
    if (!f) goto error;
    
    int Nlines = count_lines(f);
    if (Nlines <= 0) goto error;

    out.p = malloc(sizeof(s_point) * Nlines);
    if (!out.p) goto error;

    for (int ii=0; ii<Nlines; ii++)
        if (fscanf(f, "%lf,%lf,%lf", &out.p[ii].x, &out.p[ii].y, &out.p[ii].z) != 3) 
            goto error;
    
    if (fclose(f) == EOF) { f = NULL; goto error; }

    out.N = Nlines;
    return out;

    error:
        fprintf(stderr, "Error in 'read_points_from_csv'\n");
        if (f) fclose(f);
        if (out.p) free(out.p);
        return points_NAN;
}


int write_points_to_csv(const char *file, const char *f_access_mode, const s_points *points)
{
    FILE *f = fopen(file, f_access_mode);
    if (!f) goto error;

    for (int ii=0; ii<points->N; ii++)
        if (fprintf(f, "%f, %f, %f\n", points->p[ii].x, points->p[ii].y, points->p[ii].z) < 0)
            goto error;

    if (fclose(f) == EOF) { f = NULL;  goto error; }
    return 1;

    error:
        fprintf(stderr, "Error in 'write_points_to_csv'\n");
        if (f) fclose(f);
        return 0;
}


s_points copy_points(const s_points *points)
{
    s_points copy = {.N = points->N,
                     .p = malloc(sizeof(s_point) * points->N)};
    if (!copy.p) goto error;
    memcpy(copy.p, points->p, sizeof(s_point) * points->N);
    return copy;

    error:
        fprintf(stderr, "Error in 'copy_points'\n");
        return points_NAN;
}


s_points copy_points_remove_duplicates(const s_points *points, double TOL)
{
    int *mark_dup = NULL;   s_points out = {0};
    const double TOL2 = TOL * TOL;

    mark_dup = malloc(sizeof(int) * points->N);
    if (!mark_dup) goto error;
    memset(mark_dup, 0, sizeof(int) * points->N);

    for(int ii=0; ii<points->N-1; ii++) {
        for (int jj=ii+1; jj<points->N; jj++) {
            double d2 = distance_squared(points->p[ii], points->p[jj]);
            if (d2 <= TOL2) mark_dup[jj] = 1;
        }
    }

    int count_dup = 0;
    for (int ii=0; ii<points->N; ii++)
        if (mark_dup[ii] == 1) count_dup++;

    out.N = points->N - count_dup, 
    out.p = malloc(sizeof(s_point) * (points->N - count_dup));
    if (!out.p) goto error;

    int jj = 0;
    for (int ii=0; ii<points->N; ii++) 
        if (mark_dup[ii] == 0) out.p[jj++] = points->p[ii];

    free(mark_dup);
    return out;

    error:
        if (mark_dup) free(mark_dup);
        if (out.p) free(out.p);
        fprintf(stderr, "Error in 'copy_points_remove_duplicates'\n");
        return points_NAN;
}


s_point sum_points(s_point u, s_point v)
{
    return (s_point){{{u.x + v.x, u.y + v.y, u.z + v.z}}};
}


s_point subtract_points(s_point u, s_point v)
{
    return (s_point){{{u.x - v.x, u.y - v.y, u.z - v.z}}};
}


s_point spherical_to_cartesian(double r, double theta, double phi) 
{
    s_point cart;
    cart.x = r * sin(theta) * cos(phi);
    cart.y = r * sin(theta) * sin(phi);
    cart.z = r * cos(theta);
    return cart;
}

double dot_prod(s_point u, s_point v)
{
    return u.x*v.x + u.y*v.y + u.z*v.z;
}


s_point cross_prod(s_point u, s_point v)
{
    s_point out;
    out.x = u.y*v.z - u.z*v.y;
    out.y = u.z*v.x - u.x*v.z;
    out.z = u.x*v.y - u.y*v.x;
    return out;
}


double norm_squared(s_point v)
{
    return dot_prod(v, v);
}


double norm(s_point v)
{
    return sqrt(norm_squared(v));
}


double distance_squared(s_point a, s_point b)
{
    s_point diff = subtract_points(a, b);
    return norm_squared(diff);
}


double distance(s_point a, s_point b)
{
    return sqrt(distance_squared(a, b));
}


double max_distance(const s_points *points, s_point query)
{   
    double maxd2 = 0;
    for (int ii=0; ii<points->N; ii++) {
        double d2 = distance_squared(points->p[ii], query);
        if (maxd2 < d2) maxd2 = d2;
    }
    return sqrt(maxd2);
}


s_point scale_point(s_point a, double s)
{
    s_point out = a;
    out.x *= s;
    out.y *= s;
    out.z *= s;
    return out;
}


s_point normalize_vec(s_point v, double EPS)
{
    double n = norm(v);
    if (n < EPS) return point_NAN;
    return scale_point(v, 1.0/n);
}


s_point interpolate_points(s_point a, s_point b, double t)
{
    s_point diff = subtract_points(b, a);
    return sum_points(a, scale_point(diff, t));
}


void bounding_box_points(const s_points *points, s_point *min_out, s_point *max_out)
{
    if (points->N == 0) {
        *min_out = (s_point){0};
        *max_out = (s_point){0};
        return;
    }

    s_point min, max;
    min.x = points->p[0].x;   min.y = points->p[0].y;   min.z = points->p[0].z;
    max.x = points->p[0].x;   max.y = points->p[0].y;   max.z = points->p[0].z;
    if (points->N == 1) {
        *min_out = min;
        *max_out = max;
        return;
    }

    for (int ii=1; ii<points->N; ii++) {
        if (points->p[ii].x < min.x) 
            min.x = points->p[ii].x;
        if (points->p[ii].y < min.y) 
            min.y = points->p[ii].y;
        if (points->p[ii].z < min.z) 
            min.z = points->p[ii].z;

        if (points->p[ii].x > max.x)
            max.x = points->p[ii].x;
        if (points->p[ii].y > max.y) 
            max.y = points->p[ii].y;
        if (points->p[ii].z > max.z) 
            max.z = points->p[ii].z;
    }
    
    *min_out = min;
    *max_out = max;
}


s_point span_points(const s_points *points)
{
    s_point min, max;
    bounding_box_points(points, &min, &max);
    return subtract_points(max, min);
}


double absolute_tolerance_from_scale(const s_points *points, double rel_tol)
{
    assert(rel_tol >= 0);
    if (!points_is_valid(points)) return NAN;

    s_point span = span_points(points);
    double max_range = fmax(span.x, fmax(span.y, span.z));

    // If points are degenerate (all equal), define a minimal scale
    if (max_range == 0) max_range = 1;

    return rel_tol * max_range;
}


s_point point_average(const s_points *points)
{
    s_point out = points->p[0];
    for (int ii=1; ii<points->N; ii++) {
        out.x += points->p[ii].x;
        out.y += points->p[ii].y;
        out.z += points->p[ii].z;
    }
    out.x /= points->N;
    out.y /= points->N;
    out.z /= points->N;
    return out;
}


int coord_with_largest_component_3D(s_point v)
{
    double a = fabs(v.coords[0]);
    double b = fabs(v.coords[1]);
    double c = fabs(v.coords[2]);

    if (a >= b && a >= c) return 0;
    if (b >= a && b >= c) return 1;
    return 2;
}


int coord_with_smallest_component_3D(s_point v)
{
    double a = fabs(v.coords[0]);
    double b = fabs(v.coords[1]);
    double c = fabs(v.coords[2]);

    if (a <= b && a <= c) return 0;
    if (b <= a && b <= c) return 1;
    return 2;
}


s_point random_point_uniform_3D(s_point min, s_point max)
{
    double ux = rand() / ((double) RAND_MAX + 1.0);
    double uy = rand() / ((double) RAND_MAX + 1.0);
    double uz = rand() / ((double) RAND_MAX + 1.0);
    s_point out;
    out.x = min.x + (max.x - min.x) * ux;
    out.y = min.y + (max.y - min.y) * uy;
    out.z = min.z + (max.z - min.z) * uz;
    return out;
}


double area_triangle(const s_point triangle[3])
{
    return norm(cross_prod(subtract_points(triangle[2], triangle[0]), 
                           subtract_points(triangle[1], triangle[0]))) / 2.0;
}


double signed_volume_tetra(const s_point tetra[4])
{
    s_point v0 = subtract_points(tetra[0], tetra[3]);
    s_point v1 = subtract_points(tetra[1], tetra[3]);
    s_point v2 = subtract_points(tetra[2], tetra[3]);
    return 1.0 / 6.0 * (dot_prod(cross_prod(v0, v1), v2));
}


int basis_vectors_plane(const s_point plane[3], double EPS_degenerate, s_point *out_n, s_point *out_t1, s_point *out_t2)
{
    s_point d1 = subtract_points(plane[1], plane[0]);
    s_point d2 = subtract_points(plane[2], plane[0]);
    s_point n = cross_prod(d1, d2);
    *out_n = normalize_vec(n, EPS_degenerate);
    if (!point_is_valid(*out_n)) {
        *out_t1 = (s_point){0};
        *out_t2 = (s_point){0};
        return 0;
    }

    int ref_coord = coord_with_smallest_component_3D(n);
    s_point ref = (ref_coord == 0) ?   (s_point){{{1,0,0}}} :
                  ( (ref_coord == 1) ? (s_point){{{0,1,0}}} :
                                       (s_point){{{0,0,1}}} );
    *out_t1 = normalize_vec(cross_prod(ref, *out_n), EPS_degenerate);
    *out_t2 = normalize_vec(cross_prod(*out_n, *out_t1), EPS_degenerate);
    return 1;
}


int plane_equation_from_points(const s_point plane[3], double EPS_degenerate, s_point *abc_out, double *d_out)
{   /* Plane: x  s.t.  dot(abc_out, x) + d_out = 0; */
    s_point u = subtract_points(plane[1], plane[0]);
    s_point v = subtract_points(plane[2], plane[0]);
    
    s_point n = cross_prod(u, v);

    *abc_out = normalize_vec(n, EPS_degenerate);
    if (!point_is_valid(*abc_out)) { *d_out = NAN; return 0; }
    *d_out = -dot_prod(*abc_out, plane[0]);
    return 1;
}


s_point project_point_to_plane(s_point p, const s_point plane[3], double EPS_degenerate) 
{
    s_point n = cross_prod(subtract_points(plane[1], plane[0]),
                           subtract_points(plane[2], plane[0]));
    n = normalize_vec(n, EPS_degenerate);
    if (!point_is_valid(n)) return point_NAN; 

    double d = dot_prod(n, plane[0]);
    double dist = dot_prod(p, n) - d;
    return subtract_points(p, scale_point(n, dist));
}


s_point closest_point_on_triangle(const s_point triangle[3], double EPS_degenerate, s_point p)
{   /* From "C Ericson. Real-Time Collision Detection" */
    if (area_triangle(triangle) < EPS_degenerate) return point_NAN;

    s_point A = triangle[0], B = triangle[1], C = triangle[2];
    s_point AB = {{{B.x-A.x, B.y-A.y, B.z-A.z}}};
    s_point AC = {{{C.x-A.x, C.y-A.y, C.z-A.z}}};
    s_point AP = {{{p.x-A.x, p.y-A.y, p.z-A.z}}};

    // In vertex A
    double d1 = dot_prod(AB, AP), d2 = dot_prod(AC, AP);
    if (d1 <= 0 && d2 <= 0) return A;

    // In vertex B
    s_point BP = {{{ p.x-B.x, p.y-B.y, p.z-B.z }}};
    double d3 = dot_prod(AB, BP), d4 = dot_prod(AC, BP);
    if (d3 >= 0 && d4 <= d3) return B;

    // In edge AB
    double vc = d1*d4 - d3*d2;
    if (vc <= 0 && d1 >= 0 && d3 <= 0) {
        double v = d1 / (d1 - d3);
        return (s_point){{{A.x + v*AB.x, A.y + v*AB.y, A.z + v*AB.z}}};
    }

    // In vertex C
    s_point CP = {{{ p.x-C.x, p.y-C.y, p.z-C.z }}};
    double d5 = dot_prod(AB, CP), d6 = dot_prod(AC, CP);
    if (d6 >= 0 && d5 <= d6) return C;

    // In edge AC
    double vb = d5*d2 - d1*d6;
    if (vb <= 0 && d2 >= 0 && d6 <= 0) {
        double w = d2 / (d2 - d6);
        return (s_point){{{A.x + w*AC.x, A.y + w*AC.y, A.z + w*AC.z}}};
    }

    // In edge BC
    double va = d3*d6 - d5*d4;
    if (va <= 0 && (d4 - d3) >= 0 && (d5 - d6) >= 0) {
        s_point BC = {{{ C.x-B.x, C.y-B.y, C.z-B.z }}};
        double w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return (s_point){{{B.x + w*BC.x, B.y + w*BC.y, B.z + w*BC.z}}};
    }

    // Inside face region
    // Barycentric coordinates (u,v,w)
    double denom = va + vb + vc;
    if (denom < EPS_degenerate) return point_NAN;
    double v = vb / denom;
    double w = vc / denom;
    return (s_point){{{A.x + AB.x*v + AC.x*w, A.y + AB.y*v + AC.y*w, A.z + AB.z*v + AC.z*w}}};
}


s_point closest_point_on_segment(const s_point segment[2], double EPS_degenerate, s_point p)
{
    s_point AB = subtract_points(segment[1], segment[0]);
    s_point pA = subtract_points(p, segment[0]);
    // Project c onto ab, computing parameterized position d(t) = a + t*(b – a)
    double denom = dot_prod(AB, AB);
    if (denom < EPS_degenerate) return point_NAN;
    double t = dot_prod(pA, AB) / denom;
    // If outside segment, clamp t (and therefore d) to the closest endpoint
    if (t < 0) t = 0;
    if (t > 1) t = 1;
    return (s_point){{{segment[0].x + t*AB.x, segment[0].y + t*AB.y, segment[0].z + t*AB.z}}};
}


