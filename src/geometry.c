#include "geometry.h"
#include "predicates.h"
#include "string.h"
#include "stdlib.h"
#include "stdio.h"
#include <math.h>
#include <assert.h>



s_points copy_points(const s_points *points)
{
    s_points copy = {.N = points->N,
                     .p = malloc(sizeof(s_point) * points->N)};
    memcpy(copy.p, copy.p, sizeof(s_point) * points->N);
    return copy;
}


void free_points(s_points *points)
{
    free(points->p);
    memset(points, 0, sizeof(s_points));
}


int orientation(const s_point p[3], s_point q)
{
    double aux = orient3d(p[0].coords, p[1].coords, p[2].coords, q.coords);
    if (aux > 0) return 1;
    else if (aux < 0) return -1;
    else return 0;
}


int in_sphere(const s_point p[4], s_point q)
{   

    int factor;
    if (orientation(p, p[3]) == 1) factor = 1;
    else factor = -1;

    double aux = insphere(p[0].coords, p[1].coords, p[2].coords, p[3].coords, q.coords);
    
    if (aux > 0) return factor;
    else if (aux < 0) return -factor;
    else return 0;
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


s_point normalize_3d(s_point v)
{
    double n = norm(v);
    return (s_point){{{v.x/n, v.y/n, v.z/n}}};
}


s_point find_center_mass(const s_points *points)
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


int coord_with_largest_component_3d(s_point x)
{
    double *n = x.coords;
    int out = 2;
    if (fabs(n[0]) > fabs(n[1]) && fabs(n[0]) > fabs(n[2])) out = 0;
    else if (fabs(n[1]) > fabs(n[0]) && fabs(n[1]) > fabs(n[2])) out = 1;
    return out;
}


s_point random_point_uniform_3d(s_point min, s_point max)
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


int segment_crosses_triangle_3d(const s_point triangle[3], s_point a, s_point b)
{
    if (orientation(triangle, a) == orientation(triangle, b)) return 0;
    s_point aux[3];
    aux[2] = a; 

    aux[0] = triangle[0];   aux[1] = triangle[1];
    int s1 = orientation(aux, b);

    aux[0] = triangle[1];   aux[1] = triangle[2];
    int s2 = orientation(aux, b);

    aux[0] = triangle[2];   aux[1] = triangle[0];
    int s3 = orientation(aux, b);
    
    if (s1 == s2 && s2 == s3 && s3 == s1) return 1;
    else return 0;
}


int between_1d(double x, double a, double b, double eps)
{
    // returns 1 if x lies between [a,b]
    if (a > b) { double t=a; a=b; b=t; }
    return (x + eps >= a && x - eps <= b);
}


int segments_intersect_2d(const s_point AB[2], const s_point pd[2])
{
    const double EPS = 1e-9;
    double Ax = AB[0].x, Ay = AB[0].y;
    double Bx = AB[1].x, By = AB[1].y;
    double px = pd[0].x, py = pd[0].y;
    double dx = pd[1].x, dy = pd[1].y;

    // 1) Degenerate pd, treat as X
    if (fabs(px - dx) < EPS && fabs(py - dy) < EPS) {
        double X[2] = { px, py };
        if (fabs(Ax - Bx) < EPS && fabs(Ay - By) < EPS) {
            // If A, B also degenerate
            return (fabs(px - Ax) < EPS && fabs(py - Ay) < EPS);
        }
        if (orient2d(AB[0].coords, AB[1].coords, X) == 0 &&
            between_1d(px, Ax, Bx, EPS) &&
            between_1d(py, Ay, By, EPS)) {
            return 1;
        }
        return 0;
    }

    // 2) Degenerate AB, treat as Y
    if (fabs(Ax - Bx) < EPS && fabs(Ay - By) < EPS) {
        double Y[2] = { Ax, Ay };
        if (orient2d(pd[0].coords, pd[1].coords, Y) == 0 &&
            between_1d(Ax, px, dx, EPS) &&
            between_1d(Ay, py, dy, EPS)) {
            return 1;
        }
        return 0;
    }

    // 3) General case: two straddling tests
    // 3a) [A,B] vs {p,d}
    int o1 = orient2d(AB[0].coords, AB[1].coords, pd[0].coords);
    int o2 = orient2d(AB[0].coords, AB[1].coords, pd[1].coords);
    if (o1 != 0 && o2 != 0 && o1 == o2) return 0;

    // 3b) [p,d] vs {A,B}
    o1 = orient2d(pd[0].coords, pd[1].coords, AB[0].coords);
    o2 = orient2d(pd[0].coords, pd[1].coords, AB[1].coords);
    if (o1 != 0 && o2 != 0 && o1 == o2) return 0;

    return 1;
}


s_point closest_point_on_triangle(const s_point triangle[3], s_point p)
{
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
    assert(fabs(denom) > 1e-9  && "Triangle is degenerate?");  // Triangle is degenarate?
    double v = vb / denom;
    double w = vc / denom;
    return (s_point){{{A.x + AB.x*v + AC.x*w, A.y + AB.y*v + AC.y*w, A.z + AB.z*v + AC.z*w}}};
}


s_point closest_point_on_segment(const s_point segment[2], s_point p)
{
    s_point AB = subtract_points(segment[1], segment[0]);
    s_point pA = subtract_points(p, segment[0]);
    // Project c onto ab, computing parameterized position d(t) = a + t*(b – a)
    double denom = dot_prod(AB, AB);
    assert(denom > 1e-9);
    double t = dot_prod(pA, AB) / denom;
    // If outside segment, clamp t (and therefore d) to the closest endpoint
    if (t < 0) t = 0;
    if (t > 1) t = 1;
    return (s_point){{{segment[0].x + t*AB.x, segment[0].y + t*AB.y, segment[0].z + t*AB.z}}};
}


int point_in_triangle_2d(const s_point triangle[3], s_point p)
{
    int o1 = orient2d(triangle[0].coords, triangle[1].coords, p.coords);
    int o2 = orient2d(triangle[1].coords, triangle[2].coords, p.coords);
    int o3 = orient2d(triangle[2].coords, triangle[0].coords, p.coords);
    
    // Find reference sign (non-zero)
    int signs[3] = {o1, o2, o3};
    int ref_sign = 0;
    for (int ii=0; ii<3; ii++) {
        if (signs[ii] != 0) {
            ref_sign = signs[ii];
            break;
        }
    }
    assert(ref_sign != 0);

    for (int ii=0; ii<3; ii++) {
        if (signs[ii] != 0 && signs[ii] != ref_sign) return 0;
    }
    return 1;
}


int point_in_triangle_3d(const s_point triangle[3], s_point p)
{
    s_point aux = closest_point_on_triangle(triangle, p);
    if (distance_squared(aux, p) < 1e-6) return 1;
    else return 0;
}


s_points remove_duplicate_points(const s_points *points, double tol_dist)
{
    double tol2 = tol_dist * tol_dist;

    int *mark_dup = malloc(sizeof(int) * points->N);
    memset(mark_dup, 0, sizeof(int) * points->N);

    for(int ii=0; ii<points->N-1; ii++) {
        for (int jj=ii+1; jj<points->N; jj++) {
            double d2 = distance_squared(points->p[ii], points->p[jj]);
            if (d2 <= tol2) mark_dup[jj] = 1;
        }
    }

    int count_dup = 0;
    for (int ii=0; ii<points->N; ii++) {
        if (mark_dup[ii] == 1) count_dup++;
    }

    s_points out = { .N = points->N - count_dup, 
                     .p = malloc(sizeof(s_point) * (points->N - count_dup)) };
    int jj = 0;
    for (int ii=0; ii<points->N; ii++) {
        if (mark_dup[ii] == 0) out.p[jj++] = points->p[ii];
    }

    free(mark_dup);
    return out;
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
    s_points out = {0, NULL};

    FILE *f = fopen(file, "r");
    if (!f) goto error;
    
    int Nlines = count_lines(f);
    if (Nlines <= 0) goto error;

    out.p = malloc(sizeof(s_point) * Nlines);
    if (!out.p) goto error;

    for (int ii=0; ii<Nlines; ii++) {
        if (fscanf(f, "%lf,%lf,%lf", &out.p[ii].x, &out.p[ii].y, &out.p[ii].z) != 3) 
            goto error;
    }

    if (fclose(f) == EOF) { f = NULL; goto error; }

    out.N = Nlines;
    return out;

    error:
        fprintf(stderr, "Error in 'read_points_from_csv'\n");
        if (f) fclose(f);
        if (out.p) { free(out.p);  out.p = NULL;  out.N = 0; }
        return out;
}


int write_points_to_csv(const char *file, const char *f_access_mode, const s_points *points)
{
    FILE *f = fopen(file, f_access_mode);
    if (!f) goto error;

    for (int ii=0; ii<points->N; ii++) {
        if (fprintf(f, "%f, %f, %f\n", points->p[ii].x, points->p[ii].y, points->p[ii].z) < 0)
            goto error;
    }

    if (fclose(f) == EOF) { f = NULL;  goto error; }
    return EXIT_SUCCESS;

    error:
        fprintf(stderr, "Error in 'write_points_to_csv'\n");
        if (f) fclose(f);
        return EXIT_FAILURE;
}


double volume_tetrahedron_approx(s_point p1, s_point p2, s_point p3, s_point p4)
{
    return fabs(1.0/6.0 * orient3d(p1.coords, p2.coords, p3.coords, p4.coords));
}


