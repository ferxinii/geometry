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

// lp3_feasible_TTS_T -- (T,T,S) vertex [bisector s->t cutting tet edge AB] tested
// against a query tet face Tc=(Pc,Qc,Rc): -sign(gA*wB-gB*wA)*sign(gA-gB), with
// gX=|X-s|^2-|X-t|^2 and wX=n_c.(X-Pc), n_c=(Qc-Pc)x(Rc-Pc).  Surface edge
// ordering (lp3_predicates.tex Prop "TTS vs T, signed").  Caller applies the
// face interior-orientation sign sigma_f.
int lp3_feasible_TTS_T(double Ax, double Ay, double Az,
                       double Bx, double By, double Bz,
                       double sx, double sy, double sz,
                       double tx, double ty, double tz,
                       double Pcx, double Pcy, double Pcz,
                       double Qcx, double Qcy, double Qcz,
                       double Rcx, double Rcy, double Rcz);

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

// lp3_feasible_TSS_T_gen -- (T,S,S) vertex [Voronoi edge s->t1 ^ s->t2 meeting
// definer face T1=(A,Q1,R1)] tested against an ARBITRARY query face T2=(P2,Q2,R2)
// (no shared edge; surface edge ordering). sign(Delta_TSS)*sign(Gamma), Gamma =
// -c1*det[t2-s;nT1;nT2] + c2*det[t1-s;nT1;nT2] + (nT2.(P2-A))*Delta_TSS.  Degree 7.
int lp3_feasible_TSS_T_gen(double Ax, double Ay, double Az,
                           double Q1x, double Q1y, double Q1z,
                           double R1x, double R1y, double R1z,
                           double P2x, double P2y, double P2z,
                           double Q2x, double Q2y, double Q2z,
                           double R2x, double R2y, double R2z,
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

// incircle3d -- exact "is d on the circumcircle of triangle (a,b,c)", for four
// points that the caller GUARANTEES are coplanar.  Returns 0 iff d lies on that
// circumcircle, nonzero otherwise (only the zero/nonzero result is meaningful --
// the sign of the nonzero case is unspecified).
//
// It is the 3D insphere determinant of (a,b,c, w, d) with w = a + (b-a)x(c-a)
// formed EXACTLY, reduced under the coplanarity precondition: the full insphere
// determinant equals  deg6expr + g*orient3d(a,b,c,d)  where the dropped term g
// carries the expensive |n|^2 sub-expression; since orient3d(a,b,c,d)==0 by
// precondition it vanishes, leaving a degree-6 expression (same shape/size as
// lp3_feasible_SSS_T's factor2).  Because that term is dropped, the result is
// ONLY valid for coplanar input -- always gate with orient3d(a,b,c,d)==0 first.
//
// PRECONDITIONS: (1) a,b,c,d coplanar; (2) a,b,c non-collinear.  If a,b,c are
// collinear (no circumcircle) the predicate returns the sentinel 2.
int incircle3d(double ax, double ay, double az,
               double bx, double by, double bz,
               double cx, double cy, double cz,
               double dx, double dy, double dz);

#ifdef __cplusplus
}
#endif

#endif
