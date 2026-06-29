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

// CDT certify predicate -- sign of D' (time-derivative of insphere det)
// sign of det | a-Rv  dot(Rv-a, n) |  for rows {a,b,c,L}, n=(pb-pa)x(pc-pa)
//             | b-Rv  dot(Rv-b, n) |
//             | c-Rv  dot(Rv-c, n) |
//             | L-Rv  dot(Rv-L, n) |
int cdt_dprime(double ax,  double ay,  double az,
               double bx,  double by,  double bz,
               double cx,  double cy,  double cz,
               double Lx,  double Ly,  double Lz,
               double Rvx, double Rvy, double Rvz,
               double pax, double pay, double paz,
               double pbx, double pby, double pbz,
               double pcx, double pcy, double pcz);

// CDT heap ordering predicate -- sign of D0_B*D'_A - D0_A*D'_B
// used to compare flip times τ_A < τ_B without floating-point τ values
int cdt_sign_cross_dprime(double aAx,  double aAy,  double aAz,
                      double bAx,  double bAy,  double bAz,
                      double cAx,  double cAy,  double cAz,
                      double LAx,  double LAy,  double LAz,
                      double RvAx, double RvAy, double RvAz,
                      double aBx,  double aBy,  double aBz,
                      double bBx,  double bBy,  double bBz,
                      double cBx,  double cBy,  double cBz,
                      double LBx,  double LBy,  double LBz,
                      double RvBx, double RvBy, double RvBz,
                      double pax,  double pay,  double paz,
                      double pbx,  double pby,  double pbz,
                      double pcx,  double pcy,  double pcz);

// Weight-perturbed variant for SoS cascade.
// kA[5] = {k_aA, k_bA, k_cA, k_LA, k_RvA}, kB[5] = {k_aB, ...}.
// Each k is added to the corresponding vertex's weight in the D0 insphere
// determinant (kA[4]/kB[4] subtracts from all rows as the Rv reference).
// D' is unaffected by weights. Normal call: pass {0,0,0,0,0} for both arrays.
int cdt_sign_cross_dprime_weighted(double aAx,  double aAy,  double aAz,
                      double bAx,  double bAy,  double bAz,
                      double cAx,  double cAy,  double cAz,
                      double LAx,  double LAy,  double LAz,
                      double RvAx, double RvAy, double RvAz,
                      double aBx,  double aBy,  double aBz,
                      double bBx,  double bBy,  double bBz,
                      double cBx,  double cBy,  double cBz,
                      double LBx,  double LBy,  double LBz,
                      double RvBx, double RvBy, double RvBz,
                      double pax,  double pay,  double paz,
                      double pbx,  double pby,  double pbz,
                      double pcx,  double pcy,  double pcz,
                      const double kA[5], const double kB[5]);

#ifdef __cplusplus
}
#endif

#endif
