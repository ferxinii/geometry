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

// Tier B: 3D LP-feasibility predicates (lp3_*). See lp3_predicates.tex.
// lp3_D_TSS -- slope of triple (T, S_i, S_j): sign det[n_T; t1-s; t2-s].
// n_T = (Q-A) x (R-A) for tet face (A,Q,R).
int lp3_D_TSS(double Ax, double Ay, double Az,
              double Qx, double Qy, double Qz,
              double Rx, double Ry, double Rz,
              double sx, double sy, double sz,
              double t1x, double t1y, double t1z,
              double t2x, double t2y, double t2z);

// lp3_feasible_TTS_S -- (T,T,S) vertex [bisector s->t cutting tet edge AB]
// tested against bisector s->u.
int lp3_feasible_TTS_S(double Ax, double Ay, double Az,
                       double Bx, double By, double Bz,
                       double sx, double sy, double sz,
                       double tx, double ty, double tz,
                       double ux, double uy, double uz);

// lp3_slope_TTS -- sign(gA-gB), gX=|X-s|^2-|X-t|^2: 1D-edge slope of bisector
// s->t on tet edge AB (LO bound if >0, HI bound if <0). For the (d) 1D-LP.
int lp3_slope_TTS(double Ax, double Ay, double Az,
                  double Bx, double By, double Bz,
                  double sx, double sy, double sz,
                  double tx, double ty, double tz);

// lp3_det_TTS -- sign(g1A*g2B - g1B*g2A): orders the two bisector crossings of
// s->t1 and s->t2 on edge AB.  For the (d) 1D-LP.
int lp3_det_TTS(double Ax, double Ay, double Az,
                double Bx, double By, double Bz,
                double sx, double sy, double sz,
                double t1x, double t1y, double t1z,
                double t2x, double t2y, double t2z);

// lp3_feasible_TSS_S -- (T,S,S) vertex [Voronoi edge s->t1 ^ s->t2 meeting tet
// face T1=(A,Q,R)] tested against bisector s->u.
int lp3_feasible_TSS_S(double Ax, double Ay, double Az,
                       double Qx, double Qy, double Qz,
                       double Rx, double Ry, double Rz,
                       double sx, double sy, double sz,
                       double t1x, double t1y, double t1z,
                       double t2x, double t2y, double t2z,
                       double ux, double uy, double uz);

// lp3_feasible_TSS_T -- (T,S,S) vertex [Voronoi edge s->t1 ^ s->t2 meeting the
// vertex-defining face T1=(A,E,C)] tested against adjacent query face T2=(A,E,D)
// (shared edge A-E). Returns sign(Delta_TSS)*sign(Gamma) with n_T2=(E-A)x(D-A);
// caller multiplies by the T2 interior-orientation sign.
int lp3_feasible_TSS_T(double Ax, double Ay, double Az,
                       double Ex, double Ey, double Ez,
                       double Cx, double Cy, double Cz,
                       double Dx, double Dy, double Dz,
                       double sx, double sy, double sz,
                       double t1x, double t1y, double t1z,
                       double t2x, double t2y, double t2z);

// lp3_feasible_SSS_T -- (S,S,S) vertex [circumcenter of Delaunay tet
// (s,t1,t2,t3)] tested against tet face T=(A,Q,R). Returns
// sign(Delta_SSS)*sign(Gamma) with n_T=(Q-A)x(R-A); caller multiplies by the
// face interior-orientation sign.
int lp3_feasible_SSS_T(double sx, double sy, double sz,
                       double t1x, double t1y, double t1z,
                       double t2x, double t2y, double t2z,
                       double t3x, double t3y, double t3z,
                       double Ax, double Ay, double Az,
                       double Qx, double Qy, double Qz,
                       double Rx, double Ry, double Rz);

#ifdef __cplusplus
}
#endif

#endif
