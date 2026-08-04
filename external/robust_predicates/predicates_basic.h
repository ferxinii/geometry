#ifndef BASIC_PREDICATES_H
#define BASIC_PREDICATES_H

#ifdef __cplusplus
extern "c" {
#endif

// Orientation predicates
int orient2d(double ax, double ay,
             double bx, double by,
             double cx, double cy);
int orient3d(double ax, double ay, double az,
             double bx, double by, double bz,
             double cx, double cy, double cz,
             double dx, double dy, double dz);

// Incircle / insphere with orient test implicit
int powertest_n2_k3_unweighted(double ax, double ay,
                               double bx, double by,
                               double cx, double cy,
                               double dx, double dy);
int powertest_n3_k4_unweighted(double ax, double ay, double az,
                               double bx, double by, double bz,
                               double cx, double cy, double cz,
                               double dx, double dy, double dz,
                               double ex, double ey, double ez);


// Always gate with orient3d(a,b,c,d)==0 first.
// PRECONDITIONS: (1) a,b,c,d coplanar; (2) a,b,c non-collinear.  If a,b,c are
// collinear (no circumcircle) the predicate returns the sentinel 2.
int incircle3d(double ax, double ay, double az,
               double bx, double by, double bz,
               double cx, double cy, double cz,
               double dx, double dy, double dz);

// orient3d_dd -- sign of det3(a-b, c-d, e-f): 
int orient3d_dd(double ax, double ay, double az,
                double bx, double by, double bz,
                double cx, double cy, double cz,
                double dx, double dy, double dz,
                double ex, double ey, double ez,
                double fx, double fy, double fz);

#ifdef __cplusplus
}
#endif

#endif 
