#ifndef ROBUST_PREDICATES_H
#define ROBUST_PREDICATES_H

#ifdef __cplusplus
extern "C" {
#endif

// Orientation predicates
int orient2d(double ax, double ay,
             double bx, double by,
             double cx, double cy);
int orient3d(double ax, double ay, double az,
             double bx, double by, double bz,
             double cx, double cy, double cz,
             double dx, double dy, double dz);

// Power predicates unweighted (incircle / insphere with orient test implicit)
int powertest_n2_k3_unweighted(double ax, double ay,
                               double bx, double by,
                               double cx, double cy,
                               double dx, double dy);
int powertest_n3_k4_unweighted(double ax, double ay, double az,
                               double bx, double by, double bz,
                               double cx, double cy, double cz,
                               double dx, double dy, double dz,
                               double ex, double ey, double ez);

// Power predicates: Weighted
int powertest_n1_k1(double xa, double wa,
                    double xb, double wb);
int powertest_n1_k2(double xa, double wa,
                    double xb, double wb,
                    double xc, double wc);

int powertest_n2_k1(double ax, double ay, double wa,
                    double bx, double by, double wb);
int powertest_n2_k2(double ax, double ay, double wa,
                    double bx, double by, double wb,
                    double cx, double cy, double wc);
int powertest_n2_k3(double ax, double ay, double wa,
                    double bx, double by, double wb,
                    double cx, double cy, double wc,
                    double dx, double dy, double wd);

int powertest_n3_k1(double ax, double ay, double az, double wa,
                    double bx, double by, double bz, double wb);
int powertest_n3_k2(double ax, double ay, double az, double wa,
                    double bx, double by, double bz, double wb,
                    double cx, double cy, double cz, double wc);
int powertest_n3_k3(double ax, double ay, double az, double wa,
                    double bx, double by, double bz, double wb,
                    double cx, double cy, double cz, double wc,
                    double dx, double dy, double dz, double wd);
int powertest_n3_k4(double ax, double ay, double az, double wa,
                    double bx, double by, double bz, double wb,
                    double cx, double cy, double cz, double wc,
                    double dx, double dy, double dz, double wd,
                    double ex, double ey, double ez, double we);

// Orthohypersphere weight predicates
int orthow_n1_k1(double xa, double wa,
                 double alpha);
int orthow_n1_k2(double xa, double wa,
                 double xb, double wb,
                 double alpha);

int orthow_n2_k1(double ax, double ay, double wa,
                 double alpha);
int orthow_n2_k2(double ax, double ay, double wa,
                 double bx, double by, double wb,
                 double alpha);
int orthow_n2_k3(double ax, double ay, double wa,
                 double bx, double by, double wb,
                 double cx, double cy, double wc,
                 double alpha);

int orthow_n3_k1(double ax, double ay, double az, double wa,
                 double alpha);
int orthow_n3_k2(double ax, double ay, double az, double wa,
                 double bx, double by, double bz, double wb,
                 double alpha);
int orthow_n3_k3(double ax, double ay, double az, double wa,
                 double bx, double by, double bz, double wb,
                 double cx, double cy, double cz, double wc,
                 double alpha);
int orthow_n3_k4(double ax, double ay, double az, double wa,
                 double bx, double by, double bz, double wb,
                 double cx, double cy, double cz, double wc,
                 double dx, double dy, double dz, double wd,
                 double alpha);

#ifdef __cplusplus
}
#endif

#endif 
