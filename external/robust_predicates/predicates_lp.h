#ifndef ROBUST_PREDICATES_VORONOI_MIRRORING_H
#define ROBUST_PREDICATES_VORONOI_MIRRORING_H

#ifdef __cplusplus
extern "C" {
#endif

    
int lp2_det2(double Ax, double Ay, double Az,
            double Bx, double By, double Bz,
            double Cx, double Cy, double Cz,
            double sx, double sy, double sz,
            double tjx, double tjy, double tjz,
            double tkx, double tky, double tkz);
int lp2_det3(double Ax, double Ay, double Az,
            double Bx, double By, double Bz,
            double Cx, double Cy, double Cz,
            double sx, double sy, double sz,
            double tjx, double tjy, double tjz,
            double tkx, double tky, double tkz,
            double tlx, double tly, double tlz);

// 1D slope sign(D_{TRI*,i})
int lp2_D_T0_S(double Ax, double Ay, double Az,
              double Cx, double Cy, double Cz,
              double sx, double sy, double sz,
              double tix, double tiy, double tiz);
int lp2_D_T1_S(double Ax, double Ay, double Az,
              double Bx, double By, double Bz,
              double sx, double sy, double sz,
              double tix, double tiy, double tiz);
int lp2_D_T2_S(double Bx, double By, double Bz,
              double Cx, double Cy, double Cz,
              double sx, double sy, double sz,
              double tix, double tiy, double tiz);

int lp2_feasible_T0_T1_S(double Ax, double Ay, double Az,
                     double sx, double sy, double sz,
                     double tix, double tiy, double tiz);
int lp2_feasible_T1_T2_S(double Bx, double By, double Bz,
                     double sx, double sy, double sz,
                     double tix, double tiy, double tiz);
int lp2_feasible_T0_T2_S(double Cx, double Cy, double Cz,
                     double sx, double sy, double sz,
                     double tix, double tiy, double tiz);

int lp2_feasible_T0_S_S(double Ax, double Ay, double Az,
                    double Bx, double By, double Bz,
                    double Cx, double Cy, double Cz,
                    double sx, double sy, double sz,
                    double tjx, double tjy, double tjz,
                    double tlx, double tly, double tlz);
int lp2_feasible_T1_S_S(double Ax, double Ay, double Az,
                    double Bx, double By, double Bz,
                    double Cx, double Cy, double Cz,
                    double sx, double sy, double sz,
                    double tjx, double tjy, double tjz,
                    double tlx, double tly, double tlz);
int lp2_feasible_T2_S_S(double Ax, double Ay, double Az,
                    double Bx, double By, double Bz,
                    double Cx, double Cy, double Cz,
                    double sx, double sy, double sz,
                    double tjx, double tjy, double tjz,
                    double tlx, double tly, double tlz);

int lp2_feasible_T0_S_T1(double Ax, double Ay, double Az,
                     double Cx, double Cy, double Cz,
                     double sx, double sy, double sz,
                     double tjx, double tjy, double tjz);
int lp2_feasible_T1_S_T0(double Ax, double Ay, double Az,
                     double Bx, double By, double Bz,
                     double sx, double sy, double sz,
                     double tjx, double tjy, double tjz);
int lp2_feasible_T0_S_T2(double Ax, double Ay, double Az,
                     double Cx, double Cy, double Cz,
                     double sx, double sy, double sz,
                     double tjx, double tjy, double tjz);
int lp2_feasible_T1_S_T2(double Ax, double Ay, double Az,
                     double Bx, double By, double Bz,
                     double sx, double sy, double sz,
                     double tjx, double tjy, double tjz);
int lp2_feasible_T2_S_T0(double Bx, double By, double Bz,
                     double Cx, double Cy, double Cz,
                     double sx, double sy, double sz,
                     double tjx, double tjy, double tjz);
int lp2_feasible_T2_S_T1(double Bx, double By, double Bz,
                     double Cx, double Cy, double Cz,
                     double sx, double sy, double sz,
                     double tjx, double tjy, double tjz);

int lp2_feasible_S_S_T0(double Ax, double Ay, double Az,
                    double Bx, double By, double Bz,
                    double Cx, double Cy, double Cz,
                    double sx, double sy, double sz,
                    double tjx, double tjy, double tjz,
                    double tkx, double tky, double tkz);
int lp2_feasible_S_S_T1(double Ax, double Ay, double Az,
                    double Bx, double By, double Bz,
                    double Cx, double Cy, double Cz,
                    double sx, double sy, double sz,
                    double tjx, double tjy, double tjz,
                    double tkx, double tky, double tkz);
int lp2_feasible_S_S_T2(double Ax, double Ay, double Az,
                    double Bx, double By, double Bz,
                    double Cx, double Cy, double Cz,
                    double sx, double sy, double sz,
                    double tjx, double tjy, double tjz,
                    double tkx, double tky, double tkz);


// lp3_D_TSS -- slope of triple (T, S_i, S_j): sign det[n_T; t1-s; t2-s].
// n_T = (Q-A) x (R-A) for tet face (A,Q,R).
int lp3_D_TSS(double Ax, double Ay, double Az,
              double Qx, double Qy, double Qz,
              double Rx, double Ry, double Rz,
              double sx, double sy, double sz,
              double t1x, double t1y, double t1z,
              double t2x, double t2y, double t2z);

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

// Caller applies the face interior-orientation sign sigma_f.
int lp3_feasible_TTS_T(double Ax, double Ay, double Az,
                       double Bx, double By, double Bz,
                       double sx, double sy, double sz,
                       double tx, double ty, double tz,
                       double Pcx, double Pcy, double Pcz,
                       double Qcx, double Qcy, double Qcz,
                       double Rcx, double Rcy, double Rcz);

int lp3_feasible_TSS_S(double Ax, double Ay, double Az,
                       double Qx, double Qy, double Qz,
                       double Rx, double Ry, double Rz,
                       double sx, double sy, double sz,
                       double t1x, double t1y, double t1z,
                       double t2x, double t2y, double t2z,
                       double ux, double uy, double uz);

// caller multiplies by the T2 interior-orientation sign.
int lp3_feasible_TSS_T(double Ax, double Ay, double Az,
                       double Ex, double Ey, double Ez,
                       double Cx, double Cy, double Cz,
                       double Dx, double Dy, double Dz,
                       double sx, double sy, double sz,
                       double t1x, double t1y, double t1z,
                       double t2x, double t2y, double t2z);

// lp3_feasible_TSS_T_gen -- (T,S,S) vertex tested against an ARBITRARY query
// face T2=(P2,Q2,R2) (no shared edge; surface edge ordering). 
int lp3_feasible_TSS_T_gen(double Ax, double Ay, double Az,
                           double Q1x, double Q1y, double Q1z,
                           double R1x, double R1y, double R1z,
                           double P2x, double P2y, double P2z,
                           double Q2x, double Q2y, double Q2z,
                           double R2x, double R2y, double R2z,
                           double sx, double sy, double sz,
                           double t1x, double t1y, double t1z,
                           double t2x, double t2y, double t2z);

// Caller multiplies by the face interior-orientation sign.
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
