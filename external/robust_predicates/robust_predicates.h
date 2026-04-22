#ifndef ROBUST_PREDICATES_H
#define ROBUST_PREDICATES_H

/* All return +1, 0, -1 */

int orient2d(double ax, double ay,
             double bx, double by,
             double cx, double cy);

int orient3d(double ax, double ay, double az,
             double bx, double by, double bz,
             double cx, double cy, double cz,
             double dx, double dy, double dz);

int incircle(double ax, double ay,
             double bx, double by,
             double cx, double cy,
             double dx, double dy);

int insphere(double ax, double ay, double az,
             double bx, double by, double bz,
             double cx, double cy, double cz,
             double dx, double dy, double dz,
             double ex, double ey, double ez);

int powertest1d(double xa, double wa,
                double xb, double wb,
                double xc, double wc);

int powertest2d(double ax, double ay, double wa,
                double bx, double by, double wb,
                double cx, double cy, double wc,
                double dx, double dy, double wd);

int powertest3d(double ax, double ay, double az, double wa,
                double bx, double by, double bz, double wb,
                double cx, double cy, double cz, double wc,
                double dx, double dy, double dz, double wd,
                double ex, double ey, double ez, double we);

#endif 
