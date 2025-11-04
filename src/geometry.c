#include "geometry.h"
#include "string.h"
#include "stdlib.h"
#include "stdio.h"
#include <math.h>
#include <assert.h>
#include <stdatomic.h>
#include "float.h"



s_points copy_points(const s_points *points)
{
    s_points copy = {.N = points->N,
                     .p = malloc(sizeof(s_point) * points->N)};
    memcpy(copy.p, points->p, sizeof(s_point) * points->N);
    return copy;
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



// --------- PREDICATES -------------
static atomic_flag predicates_init_flag = ATOMIC_FLAG_INIT;
static inline void ensure_predicates_initialized(void)
{
    if (!atomic_flag_test_and_set(&predicates_init_flag)) {
        exactinit();
    }

}


int orientation(const s_point p[3], s_point q)
{
    ensure_predicates_initialized();

    double aux = orient3d(p[0].coords, p[1].coords, p[2].coords, q.coords);
    if (aux > 0) return 1;
    else if (aux < 0) return -1;
    else return 0;
}


int in_sphere(const s_point p[4], s_point q)
{   
    ensure_predicates_initialized();

    int o = orientation(p, p[3]);
    int factor;
    if (o == 0) return 0;
    else if (o == 1) factor = 1;
    else factor = -1;

    double aux = insphere(p[0].coords, p[1].coords, p[2].coords, p[3].coords, q.coords);
    
    if (aux > 0) return factor;
    else if (aux < 0) return -factor;
    else return 0;
}


// ------ POINTS -------
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


s_point normalize_3d(s_point v)
{
    double n = norm(v);
    assert(n > 1e-14 && "Degenerate vector"); 
    return scale_point(v, 1.0/n);
}


void bounding_box_points(const s_points *points, s_point *min_out, s_point *max_out)
{
    s_point min, max;
    min.x = DBL_MAX;   min.y = DBL_MAX;   min.z = DBL_MAX;
    max.x = -DBL_MAX;  max.y = -DBL_MAX;  max.z = -DBL_MAX;

    for (int ii=0; ii<points->N; ii++) {
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


int coord_with_largest_component_3d(s_point v)
{
    double a = fabs(v.coords[0]);
    double b = fabs(v.coords[1]);
    double c = fabs(v.coords[2]);

    if (a >= b && a >= c) return 0;
    if (b >= a && b >= c) return 1;
    return 2;
}


int coord_with_smallest_component_3d(s_point v)
{
    double a = fabs(v.coords[0]);
    double b = fabs(v.coords[1]);
    double c = fabs(v.coords[2]);

    if (a <= b && a <= c) return 0;
    if (b <= a && b <= c) return 1;
    return 2;
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
    const double EPS = 1e-12;
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


int in_triangle_2d(const double a[2], const double b[2], const double c[2], const double p[2])
{
    assert(orient2d(a, b, c) != 0 && "Degenerate 2D triangle");

    double o1 = orient2d(a, b, p);
    double o2 = orient2d(b, c, p);
    double o3 = orient2d(c, a, p);
    int s1 = (o1 > 0) - (o1 < 0);
    int s2 = (o2 > 0) - (o2 < 0);
    int s3 = (o3 > 0) - (o3 < 0);
    
    // Find reference sign (non-zero) (guaranteed to exist because non_degenerate != 0)
    int ref = (s1 != 0) ? s1 : ((s2 != 0) ? s2 : s3);
    assert(ref != 0);
    
    // If any non-zero orientation disagrees, it's outside
    if ((s1 != 0 && s1 != ref) || 
        (s2 != 0 && s2 != ref) || 
        (s3 != 0 && s3 != ref))
        return 0;

    // If any orientation is exactly zero, point is on an edge or vertex 
    // However, we already know the point is NOT outside , so the only possibility is that 
    // it lies in an edge or a vertex
    if (s1 == 0 || s2 == 0 || s3 == 0) return -1;

    return 1;
}


int in_triangle_3d(const s_point triangle[3], s_point p)
{
    // First chack if it is not coplanar
    if (orientation(triangle, p) != 0) return 0;

    // If it is coplanar, change to 2D coordinates in the plane
    s_point d1 = subtract_points(triangle[1], triangle[0]);
    s_point d2 = subtract_points(triangle[2], triangle[0]);
    s_point n = cross_prod(d1, d2);
    if (norm(n) < 1e-14) return 0;  // If too small, return 0
    n = normalize_3d(n);

    int ref_coord = coord_with_smallest_component_3d(n);
    s_point ref = (ref_coord == 0) ?   (s_point){{{1,0,0}}} :
                  ( (ref_coord == 1) ? (s_point){{{0,1,0}}} :
                                       (s_point){{{0,0,1}}} );
    s_point t1 = normalize_3d(cross_prod(ref, n));
    s_point t2 = normalize_3d(cross_prod(n, t1));

    // Build new triangle (lives in 2D)   
    double v1[2] = {dot_prod(triangle[0], t1), dot_prod(triangle[0], t2)};
    double v2[2] = {dot_prod(triangle[1], t1), dot_prod(triangle[1], t2)};
    double v3[2] = {dot_prod(triangle[2], t1), dot_prod(triangle[2], t2)};
    double paux[2] = {dot_prod(p, t1), dot_prod(p, t2)};

    return in_triangle_2d(v1, v2, v3, paux);
}


int in_tetrahedron(const s_point tetra[4], s_point query)
{
    s_point tmp[3];

    // First compute reference signs
    tmp[0] = tetra[1];   tmp[1] = tetra[2];   tmp[2] = tetra[3];
    int e0 = orientation(tmp, tetra[0]);

    tmp[0] = tetra[0];   tmp[1] = tetra[3];   tmp[2] = tetra[2];
    int e1 = orientation(tmp, tetra[1]);

    tmp[0] = tetra[0];   tmp[1] = tetra[1];   tmp[2] = tetra[3];
    int e2 = orientation(tmp, tetra[2]);

    tmp[0] = tetra[0];   tmp[1] = tetra[2];   tmp[2] = tetra[1];
    int e3 = orientation(tmp, tetra[3]);

    if (e0 == 0 || e1 == 0 || e2 == 0 || e3 == 0) {
        fprintf(stderr, "Tetrahedron is degenerate.\n");
        exit(1);
    }

    // Compute signs for the query
    tmp[0] = tetra[1];   tmp[1] = tetra[2];   tmp[2] = tetra[3];
    int s0 = orientation(tmp, query);

    tmp[0] = tetra[0];   tmp[1] = tetra[3];   tmp[2] = tetra[2];
    int s1 = orientation(tmp, query);

    tmp[0] = tetra[0];   tmp[1] = tetra[1];   tmp[2] = tetra[3];
    int s2 = orientation(tmp, query);

    tmp[0] = tetra[0];   tmp[1] = tetra[2];   tmp[2] = tetra[1];
    int s3 = orientation(tmp, query);

    if ((s0 != 0 && s0 * e0 < 0) ||
        (s1 != 0 && s1 * e1 < 0) ||
        (s2 != 0 && s2 * e2 < 0) ||
        (s3 != 0 && s3 * e3 < 0))
        return 0;
    
    // Count how many face-orientations are exactly zero
    int zeros = (s0 == 0) + (s1 == 0) + (s2 == 0) + (s3 == 0);
    if (zeros == 0) return 1;   // strictly inside
    if (zeros == 1) return -1;  // on face interior
    if (zeros == 2) return -2;  // on edge 
    if (zeros == 3) return -3;  // at vertex 

    // should not reach here (zeros can't be 4 because tetra not degenerate)
    fprintf(stderr, "Could not determine if point is inside tetrahedron.");
    exit(1);
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




int points_inside_halfspace(const s_point plane_ordered[3], s_points points, int out[points.N])
{   // Right hand rule: normal outwards
    int count = 0;
    memset(out, 0, sizeof(int) * points.N);

    for (int ii=0; ii<points.N; ii++) {
        int o = orientation(plane_ordered, points.p[ii]);
        if (o == 1) {
            out[ii] = 1;
            count++;
        } else if (o == 0) {
            out[ii] = -1;
            count++;
        }
    }
    return count;
}


void plane_equation_from_points(const s_point plane[3], s_point *abc_out, double *d_out)
{  // Plane: x  s.t.  dot(abc_out, x) + d = 0;
    s_point u = subtract_points(plane[1], plane[0]);
    s_point v = subtract_points(plane[2], plane[0]);
    
    s_point n = cross_prod(u, v);
    double n_len = norm(n);
    assert(n_len > 1e-12 && "Plane degenerate");

    *abc_out = scale_point(n, 1.0/n_len);
    *d_out = -dot_prod(*abc_out, plane[0]);
}


s_point interpolate_points(s_point a, s_point b, double t)
{
    s_point diff = subtract_points(b, a);
    return sum_points(a, scale_point(diff, t));
}


s_point project_point_to_plane(s_point p, const s_point plane[3]) 
{
    s_point n = cross_prod(subtract_points(plane[1], plane[0]),
                           subtract_points(plane[2], plane[0]));
    s_point n_unit = normalize_3d(n);
    double d = dot_prod(n_unit, plane[0]);

    double dist = dot_prod(p, n_unit) - d;
    s_point out = subtract_points(p, scale_point(n_unit, dist));  // snap to plane
    // assert(orientation(plane, out) == 0 && "Point is not coplanar after projecting!");
    return out;
}


int segment_plane_intersection(const s_point segment[2], const s_point plane[3], s_point out[2])
{
    // Check degenerate cases:
    int o1 = orientation(plane, segment[0]);
    int o2 = orientation(plane, segment[1]);
    if (o1 == 0 && o2 == 0) {  // Whole segment in plane
        out[0] = segment[0];
        out[1] = segment[1];
        return 2;
    } else if (o1 == 0 && o2 != 0) {  // One end in plane
        out[0] = segment[0];
        return 1;
    } else if (o2 == 0 && o1 != 0) {  // The other end in plane
        out[0] = segment[1];
        return 1;
    } else if (o1 == o2) {  // Both ends in same side of plane
        return 0;
    }

    // Near-repeated points (not captured by orientation)
    double EPS = 1e-12;
    int s0 = 0, s1 = 0;
    if (norm(subtract_points(plane[0], segment[0])) < EPS ||
        norm(subtract_points(plane[1], segment[0])) < EPS ||
        norm(subtract_points(plane[2], segment[0])) < EPS )
        s0 = 1;
    if (norm(subtract_points(plane[0], segment[1])) < EPS ||
        norm(subtract_points(plane[1], segment[1])) < EPS ||
        norm(subtract_points(plane[2], segment[1])) < EPS )
       s1 = 1;

    if (s0 == 1 && s1 == 0) {
        out[0] = project_point_to_plane(segment[0], plane);
        return 1;
    } else if (s0 == 0 && s1 == 1) {
        out[0] = project_point_to_plane(segment[1], plane);
        return 1;
    } else if (s0 == 1 && s1 == 1) {
        out[0] = project_point_to_plane(segment[0], plane);
        out[1] = project_point_to_plane(segment[1], plane);
        return 2;
    }

    // If we reach here, it means that there is a single intersection
    // Form unit normal and unit d so s=n_unit*x-d_unit is approx signed distance
    s_point n = cross_prod(subtract_points(plane[1], plane[0]), 
                           subtract_points(plane[2], plane[0]));
    s_point n_unit = normalize_3d(n);
    double d_unit = dot_prod(n_unit, plane[0]);

    // Compute signed distances (approx distances)
    double sA = dot_prod(segment[0], n_unit) - d_unit;
    double sB = dot_prod(segment[1], n_unit) - d_unit;
    
    s0 = 0;  s1 = 0;
    if (fabs(sA) < EPS) s0 = 1;   // segment[0] is almost on the plane
    if (fabs(sB) < EPS) s1 = 1;   // segment[1] is almost on the plane
    if (s0 == 1 && s1 == 0) {
        out[0] = project_point_to_plane(segment[0], plane);
        return 1;
    } else if (s0 == 0 && s1 == 1) {
        out[0] = project_point_to_plane(segment[1], plane);
        return 1;
    } else if (s0 == 1 && s1 == 1) {
        out[0] = project_point_to_plane(segment[0], plane);
        out[1] = project_point_to_plane(segment[1], plane);
        return 2;
    }


    double denom = (sA - sB);
    assert(fabs(denom) > EPS && "denom is too small");

    double t = sA / denom;
    assert(t < 1 && t > 0 && "interpolation parameter out of range");

    out[0] = interpolate_points(segment[0], segment[1], t);
    return 1;
}


