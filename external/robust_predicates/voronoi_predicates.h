#ifndef ROBUST_PREDICATES_VORONOI_MIRRORING_H
#define ROBUST_PREDICATES_VORONOI_MIRRORING_H

#ifdef __cplusplus
extern "C" {
#endif

// LP mirroring predicates -- 1D slope sign(D_{TRI*,i})
int lp_D_T0_S(double Ax, double Ay, double Az,
              double Cx, double Cy, double Cz,
              double sx, double sy, double sz,
              double tix, double tiy, double tiz);
int lp_D_T1_S(double Ax, double Ay, double Az,
              double Bx, double By, double Bz,
              double sx, double sy, double sz,
              double tix, double tiy, double tiz);
int lp_D_T2_S(double Bx, double By, double Bz,
              double Cx, double Cy, double Cz,
              double sx, double sy, double sz,
              double tix, double tiy, double tiz);

// LP mirroring predicates
int lp_det2(double Ax, double Ay, double Az,
            double Bx, double By, double Bz,
            double Cx, double Cy, double Cz,
            double sx, double sy, double sz,
            double tjx, double tjy, double tjz,
            double tkx, double tky, double tkz);
int lp_det3(double Ax, double Ay, double Az,
            double Bx, double By, double Bz,
            double Cx, double Cy, double Cz,
            double sx, double sy, double sz,
            double tjx, double tjy, double tjz,
            double tkx, double tky, double tkz,
            double tlx, double tly, double tlz);

int lp_feasible_T0_T1_S(double Ax, double Ay, double Az,
                     double sx, double sy, double sz,
                     double tix, double tiy, double tiz);
int lp_feasible_T1_T2_S(double Bx, double By, double Bz,
                     double sx, double sy, double sz,
                     double tix, double tiy, double tiz);
int lp_feasible_T0_T2_S(double Cx, double Cy, double Cz,
                     double sx, double sy, double sz,
                     double tix, double tiy, double tiz);

int lp_feasible_T0_S_S(double Ax, double Ay, double Az,
                    double Bx, double By, double Bz,
                    double Cx, double Cy, double Cz,
                    double sx, double sy, double sz,
                    double tjx, double tjy, double tjz,
                    double tlx, double tly, double tlz);
int lp_feasible_T1_S_S(double Ax, double Ay, double Az,
                    double Bx, double By, double Bz,
                    double Cx, double Cy, double Cz,
                    double sx, double sy, double sz,
                    double tjx, double tjy, double tjz,
                    double tlx, double tly, double tlz);
int lp_feasible_T2_S_S(double Ax, double Ay, double Az,
                    double Bx, double By, double Bz,
                    double Cx, double Cy, double Cz,
                    double sx, double sy, double sz,
                    double tjx, double tjy, double tjz,
                    double tlx, double tly, double tlz);

int lp_feasible_T0_S_T1(double Ax, double Ay, double Az,
                     double Cx, double Cy, double Cz,
                     double sx, double sy, double sz,
                     double tjx, double tjy, double tjz);
int lp_feasible_T1_S_T0(double Ax, double Ay, double Az,
                     double Bx, double By, double Bz,
                     double sx, double sy, double sz,
                     double tjx, double tjy, double tjz);
int lp_feasible_T0_S_T2(double Ax, double Ay, double Az,
                     double Cx, double Cy, double Cz,
                     double sx, double sy, double sz,
                     double tjx, double tjy, double tjz);
int lp_feasible_T1_S_T2(double Ax, double Ay, double Az,
                     double Bx, double By, double Bz,
                     double sx, double sy, double sz,
                     double tjx, double tjy, double tjz);
int lp_feasible_T2_S_T0(double Bx, double By, double Bz,
                     double Cx, double Cy, double Cz,
                     double sx, double sy, double sz,
                     double tjx, double tjy, double tjz);
int lp_feasible_T2_S_T1(double Bx, double By, double Bz,
                     double Cx, double Cy, double Cz,
                     double sx, double sy, double sz,
                     double tjx, double tjy, double tjz);

int lp_feasible_S_S_T0(double Ax, double Ay, double Az,
                    double Bx, double By, double Bz,
                    double Cx, double Cy, double Cz,
                    double sx, double sy, double sz,
                    double tjx, double tjy, double tjz,
                    double tkx, double tky, double tkz);
int lp_feasible_S_S_T1(double Ax, double Ay, double Az,
                    double Bx, double By, double Bz,
                    double Cx, double Cy, double Cz,
                    double sx, double sy, double sz,
                    double tjx, double tjy, double tjz,
                    double tkx, double tky, double tkz);
int lp_feasible_S_S_T2(double Ax, double Ay, double Az,
                    double Bx, double By, double Bz,
                    double Cx, double Cy, double Cz,
                    double sx, double sy, double sz,
                    double tjx, double tjy, double tjz,
                    double tkx, double tky, double tkz);

#ifdef __cplusplus
}
#endif

#endif
