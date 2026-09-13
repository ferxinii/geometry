#ifndef ROBUST_PREDICATES_LP_POWER_H
#define ROBUST_PREDICATES_LP_POWER_H

/* Weighted (power/Laguerre) analogues of the LP-feasibility predicates in
 * predicates_lp.h.  Derived in geometric_tests.tex, appendix:LP, by the
 * single substitution proved there: a bisector constraint S=(s->t)'s normal
 * a_S is untouched by weights, only its offset changes,
 * c_S -> c_S^w = c_S + 1/2*(w_s - w_t).  Every predicate below is that one
 * substitution threaded through its own already-published formula (see
 * predicates_lp_power.cpp for the per-predicate derivation and which pieces
 * are weight-invariant).
 *
 * Naming/argument convention: same point arguments as the unweighted
 * counterpart (predicates_lp.h), in the same order, with one extra double
 * weight argument appended per SEED point (never per tetrahedron vertex --
 * tet vertices carry no weight), in point-then-weight order matching
 * test_orthosphere's wc[k] alongside c[k].  At all-zero weights every
 * function here must agree in sign with its unweighted counterpart -- this
 * is the regression-equivalence check exercised by the test suite. */

#ifdef __cplusplus
extern "C" {
#endif

/* power_lp2_feasible_T0_T1_S -- weighted analogue of lp2_feasible_T0_T1_S.
 * sign(|ti-A|^2 - |s-A|^2 + ws - wt).  Used generically at arbitrary 3D
 * points by compute_vertex_owners (not only against a T0/T1 tet edge, per
 * the tet-case-table's (T,T,T;S) row: this IS that predicate). */
int power_lp2_feasible_T0_T1_S(double Ax, double Ay, double Az,
                               double sx, double sy, double sz,
                               double tix, double tiy, double tiz,
                               double ws, double wt);

/* power_lp3_feasible_TTS_S -- weighted analogue of lp3_feasible_TTS_S.
 * feasible = -sign(gA^w * hB^w - gB^w * hA^w) * sign(gA - gB), where
 * gX^w = gX - ws + wt, hX^w = hX - ws + wu (leaf-substituted, then expanded:
 * the cross terms are NOT a simple top-level correction, see .cpp).  The
 * slope factor sign(gA-gB) is weight-invariant (a_S untouched) and reused
 * unweighted via lp3_slope_TTS. */
int power_lp3_feasible_TTS_S(double Ax, double Ay, double Az,
                             double Bx, double By, double Bz,
                             double sx, double sy, double sz,
                             double tx, double ty, double tz,
                             double ux, double uy, double uz,
                             double ws, double wt, double wu);

/* power_lp3_det_TTS -- weighted analogue of lp3_det_TTS.  Same weighted
 * product as power_lp3_feasible_TTS_S's factor2 (g1 playing the role of g,
 * g2 the role of h), with t1,t2 in place of t,u. */
int power_lp3_det_TTS(double Ax, double Ay, double Az,
                      double Bx, double By, double Bz,
                      double sx, double sy, double sz,
                      double t1x, double t1y, double t1z,
                      double t2x, double t2y, double t2z,
                      double ws, double wt1, double wt2);

/* power_lp3_feasible_TTS_T -- weighted analogue of lp3_feasible_TTS_T.
 * feasible = -sign(gA^w*wB - gB^w*wA) * sign(gA-gB); gX^w = gX-ws+wt as
 * above, wX (query T-face signed distance) is weight-independent, untouched.
 * Caller applies the query face's interior-orientation sign, same as the
 * unweighted version. */
int power_lp3_feasible_TTS_T(double Ax, double Ay, double Az,
                             double Bx, double By, double Bz,
                             double sx, double sy, double sz,
                             double tx, double ty, double tz,
                             double Pcx, double Pcy, double Pcz,
                             double Qcx, double Qcy, double Qcz,
                             double Rcx, double Rcy, double Rcz,
                             double ws, double wt);

/* power_lp3_feasible_TSS_S -- weighted analogue of lp3_feasible_TSS_S.
 * feasible = sign(Delta_TSS) * sign(Gamma^w), Delta_TSS weight-invariant
 * (reused unweighted via lp3_D_TSS).  Gamma^w is linear in c1,c2,cu: direct
 * substitution ci^w = ci + ws - w_{ti} (resp. cu^w = cu + ws - wu), same
 * geometric minors M1,M2,M3 (pure orientation, untouched). */
int power_lp3_feasible_TSS_S(double Ax, double Ay, double Az,
                             double Qx, double Qy, double Qz,
                             double Rx, double Ry, double Rz,
                             double sx, double sy, double sz,
                             double t1x, double t1y, double t1z,
                             double t2x, double t2y, double t2z,
                             double ux, double uy, double uz,
                             double ws, double wt1, double wt2, double wu);

/* power_lp3_feasible_TSS_T -- weighted analogue of lp3_feasible_TSS_T.
 * Same c1^w, c2^w substitution as TSS_S; the two T-face minors M1, M2 are
 * pure orientation, untouched.  Caller applies the query face's interior
 * orientation sign, same as the unweighted version. */
int power_lp3_feasible_TSS_T(double Ax, double Ay, double Az,
                             double Ex, double Ey, double Ez,
                             double Cx, double Cy, double Cz,
                             double Dx, double Dy, double Dz,
                             double sx, double sy, double sz,
                             double t1x, double t1y, double t1z,
                             double t2x, double t2y, double t2z,
                             double ws, double wt1, double wt2);

/* power_lp3_feasible_TSS_T_gen -- weighted analogue of lp3_feasible_TSS_T_gen
 * (arbitrary, non-edge-sharing query face T2; surface edge ordering).  Same
 * c1^w, c2^w substitution; the T2-vs-A offset term "w" (and its x2 scaling)
 * is weight-independent, untouched. */
int power_lp3_feasible_TSS_T_gen(double Ax, double Ay, double Az,
                                 double Q1x, double Q1y, double Q1z,
                                 double R1x, double R1y, double R1z,
                                 double P2x, double P2y, double P2z,
                                 double Q2x, double Q2y, double Q2z,
                                 double R2x, double R2y, double R2z,
                                 double sx, double sy, double sz,
                                 double t1x, double t1y, double t1z,
                                 double t2x, double t2y, double t2z,
                                 double ws, double wt1, double wt2);

/* power_lp3_feasible_SSS_T -- weighted analogue of lp3_feasible_SSS_T.
 * feasible = sign(Delta_SSS) * sign(Gamma^w).  Delta_SSS = det[t1-s;t2-s;t3-s]
 * is weight-invariant (recomputed locally, unweighted, matching the
 * unweighted function's own convention exactly).  Gamma^w substitutes
 * ci^w = ci + ws - w_{ti} for i=1,2,3; the T-face minors M1,M2,M3 are pure
 * orientation, untouched.  Caller applies the face's interior-orientation
 * sign, same as the unweighted version.
 * NOTE: does NOT cover the (S,S,S;S') row (circumcenter vs. a fifth seed) --
 * that reduces to an orthosphere test (test_orthosphere/powertest_n3_k4),
 * not this LP-feasibility framework; see the weighted vertex-labeling work
 * that replaces test_insphere/points_concyclic. */
int power_lp3_feasible_SSS_T(double sx, double sy, double sz,
                             double t1x, double t1y, double t1z,
                             double t2x, double t2y, double t2z,
                             double t3x, double t3y, double t3z,
                             double Ax, double Ay, double Az,
                             double Qx, double Qy, double Qz,
                             double Rx, double Ry, double Rz,
                             double ws, double wt1, double wt2, double wt3);

/* power_lp2_feasible_S_S_T0/T1/T2 -- weighted analogues of
 * lp2_feasible_S_S_T0/T1/T2 (nface_face_crossings' domain-boundary-face
 * crossing test). factor1 = sign(a_j*b_k-a_k*b_j) is built only from the
 * (weight-invariant) bisector normals t-s and is reused unweighted via
 * lp2_det2. factor2 is LINEAR in c_j=|tj-A|^2-|s-A|^2 and c_k likewise (each
 * multiplied only by a weight-invariant a/b coefficient, never by each
 * other), so the direct substitution c_m^w = c_m + (ws - w_tm) is exact --
 * same linear-in-c class as power_lp3_feasible_TSS_S's Gamma. */
int power_lp2_feasible_S_S_T0(double Ax, double Ay, double Az,
                              double Bx, double By, double Bz,
                              double Cx, double Cy, double Cz,
                              double sx, double sy, double sz,
                              double tjx, double tjy, double tjz,
                              double tkx, double tky, double tkz,
                              double ws, double wtj, double wtk);
int power_lp2_feasible_S_S_T1(double Ax, double Ay, double Az,
                              double Bx, double By, double Bz,
                              double Cx, double Cy, double Cz,
                              double sx, double sy, double sz,
                              double tjx, double tjy, double tjz,
                              double tkx, double tky, double tkz,
                              double ws, double wtj, double wtk);
int power_lp2_feasible_S_S_T2(double Ax, double Ay, double Az,
                              double Bx, double By, double Bz,
                              double Cx, double Cy, double Cz,
                              double sx, double sy, double sz,
                              double tjx, double tjy, double tjz,
                              double tkx, double tky, double tkz,
                              double ws, double wtj, double wtk);

#ifdef __cplusplus
}
#endif

#endif
