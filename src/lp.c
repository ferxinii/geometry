
#include "points.h"
#include "gtests.h"
#include "predicates_lp.h"
#include "predicates_lp_power.h"   /* power_lp2_* / power_lp3_*: weighted dispatch */
#include "lp.h"
#include <stdio.h>
#include <assert.h>

#ifndef GEOMETRY_COMPILE_TESTS_POWER
#error "lp.c needs GEOMETRY_COMPILE_TESTS_POWER (weighted lp3_* use test_orthosphere)"
#endif


// ****************************************************************************
// ********************************** lp3 *************************************
// ****************************************************************************

/* ------ lp3_D : slope factor sign det[n1;n2;n3] ---------------------------- */
/* Dispatch by the type multiset of the triple. Used for slope decisions in a
 * Seidel-style solver; the enumerate-and-hull clip does not call it, but it is
 * provided for completeness and matches the lp3_D row of the catalogue.       */
int lp3_D(const s_point tet[4], s_point s, const s_points *seeds,
          s_cid3 a, s_cid3 b, s_cid3 c)
{
    int Tf[3], nT, nS; s_cid3 Sc[3];
    lp3_split3(a, b, c, Tf, &nT, Sc, &nS);

    if (nT == 0) {                                   /* SSS: orient3d(s,t1,t2,t3) */
        s_point t1 = seeds->p[Sc[0].idx], t2 = seeds->p[Sc[1].idx], t3 = seeds->p[Sc[2].idx];
        return test_orientation((s_point[]){t1, t2, t3}, s);
    }
    if (nT == 1) {                                   /* TSS: lp3_D_TSS */
        int fv[3]; lp3_face_idx(Tf[0], fv);
        s_point A = tet[fv[0]], Q = tet[fv[1]], R = tet[fv[2]];
        s_point t1 = seeds->p[Sc[0].idx], t2 = seeds->p[Sc[1].idx];
        return lp3_D_TSS(A.x,A.y,A.z, Q.x,Q.y,Q.z, R.x,R.y,R.z,
                         s.x,s.y,s.z, t1.x,t1.y,t1.z, t2.x,t2.y,t2.z);
    }
    if (nT == 2) {                                   /* TTS: (t-s).(B-A) on shared edge */
        int e[2]; lp3_shared_edge(Tf[0], Tf[1], e);
        s_point A = tet[e[0]], B = tet[e[1]], t = seeds->p[Sc[0].idx];
        /* sign((t-s).(B-A)) via lp_D_T1_S(A,B,s,t); per-edge orientation sign is
         * applied by the caller when this feeds a 1D-edge lo/hi decision. */
        return lp2_D_T1_S(A.x,A.y,A.z, B.x,B.y,B.z, s.x,s.y,s.z, t.x,t.y,t.z);
    }
    /* nT == 3, TTT: determinant of three face normals -- a fixed per-tet
     * orientation sign. Not used by the enumerate-and-hull clip. */
    return test_orientation(tet, tet[3]);
}

/* ---- lp3_feasible : feasibility of vertex (a^b^c) against query d -------- */
/* Returns >0 (strictly inside d), 0 (on d's boundary / degenerate), <0 (outside).
 * For T queries the raw predicate is multiplied by sigma_f[query face] so that
 * the tet interior is the feasible side.
 * weights: NULL-able, index-aligned with seeds->p[]; s_w = weight of s (see
 * lp.h).  When non-NULL every S-involving branch dispatches to its power_*
 * counterpart; the purely combinatorial branches (TTS vs T, TTT vs T) are
 * weight-independent.                                                         */
int lp3_feasible(const s_point tet[4], const int sigma_f[4], s_point s,
                 const s_points *seeds, const double *weights, double s_w,
                 s_cid3 a, s_cid3 b, s_cid3 c, s_cid3 d)
{
    int Tf[3], nT, nS; s_cid3 Sc[3];
    lp3_split3(a, b, c, Tf, &nT, Sc, &nS);

    if (nT == 0) {                                   /* ---- SSS vertex ---- */
        s_point t1 = seeds->p[Sc[0].idx], t2 = seeds->p[Sc[1].idx], t3 = seeds->p[Sc[2].idx];
        if (d.type == LP3_SEED) {                     /* SSS vs S == insphere */
            /* circumcenter of (s,t1,t2,t3); feasible iff u outside circumsphere.
             * (In a valid DT this never rejects.) Sign validated by
             * tests/clip_lp_validate.c. */
            s_point sph[4] = { s, t1, t2, t3 };
            s_point u = seeds->p[d.idx];
            if (weights) {
                double sph_w[4] = { s_w, weights[Sc[0].idx], weights[Sc[1].idx],
                                    weights[Sc[2].idx] };
                return -test_orthosphere(4, sph, sph_w, u, weights[d.idx]);
            }
            return -test_insphere(sph, u);
        } else {                                     /* SSS vs T */
            int fv[3]; lp3_face_idx(d.idx, fv);
            s_point A = tet[fv[0]], Q = tet[fv[1]], R = tet[fv[2]];
            if (weights) {
                double wt1 = weights[Sc[0].idx], wt2 = weights[Sc[1].idx],
                       wt3 = weights[Sc[2].idx];
                return sigma_f[d.idx] *
                       power_lp3_feasible_SSS_T(s.x,s.y,s.z, t1.x,t1.y,t1.z, t2.x,t2.y,t2.z,
                                                t3.x,t3.y,t3.z, A.x,A.y,A.z, Q.x,Q.y,Q.z, R.x,R.y,R.z,
                                                s_w, wt1, wt2, wt3);
            }
            return sigma_f[d.idx] *
                   lp3_feasible_SSS_T(s.x,s.y,s.z, t1.x,t1.y,t1.z, t2.x,t2.y,t2.z,
                                      t3.x,t3.y,t3.z, A.x,A.y,A.z, Q.x,Q.y,Q.z, R.x,R.y,R.z);
        }
    }

    if (nT == 1) {                                   /* ---- TSS vertex ---- */
        int f = Tf[0];
        s_point t1 = seeds->p[Sc[0].idx], t2 = seeds->p[Sc[1].idx];
        if (d.type == LP3_SEED) {                     /* TSS vs S */
            int fv[3]; lp3_face_idx(f, fv);
            s_point A = tet[fv[0]], Q = tet[fv[1]], R = tet[fv[2]], u = seeds->p[d.idx];
            if (weights) {
                double wt1 = weights[Sc[0].idx], wt2 = weights[Sc[1].idx],
                       wu = weights[d.idx];
                return power_lp3_feasible_TSS_S(A.x,A.y,A.z, Q.x,Q.y,Q.z, R.x,R.y,R.z,
                                                s.x,s.y,s.z, t1.x,t1.y,t1.z, t2.x,t2.y,t2.z,
                                                u.x,u.y,u.z, s_w, wt1, wt2, wu);
            }
            return lp3_feasible_TSS_S(A.x,A.y,A.z, Q.x,Q.y,Q.z, R.x,R.y,R.z,
                                      s.x,s.y,s.z, t1.x,t1.y,t1.z, t2.x,t2.y,t2.z,
                                      u.x,u.y,u.z);
        } else {                                     /* TSS vs T (query face d) */
            int g = d.idx;
            int e[2]; lp3_shared_edge(f, g, e);
            /* T1 (definer) = face f, apex = tet[g];  T2 (query) = face g, apex = tet[f]. */
            s_point A = tet[e[0]], E = tet[e[1]], C = tet[g], D = tet[f];
            /* Orientation of the query face T2 must match the predicate's own
             * n_T2 = (E-A)x(D-A) vertex order, so compute it inline rather than
             * from the ascending-order sigma_f[] (which can differ by a swap). */
            int sq = test_orientation((s_point[]){A, E, D}, tet[g]);
            if (weights) {
                double wt1 = weights[Sc[0].idx], wt2 = weights[Sc[1].idx];
                return sq *
                       power_lp3_feasible_TSS_T(A.x,A.y,A.z, E.x,E.y,E.z, C.x,C.y,C.z, D.x,D.y,D.z,
                                                s.x,s.y,s.z, t1.x,t1.y,t1.z, t2.x,t2.y,t2.z,
                                                s_w, wt1, wt2);
            }
            return sq *
                   lp3_feasible_TSS_T(A.x,A.y,A.z, E.x,E.y,E.z, C.x,C.y,C.z, D.x,D.y,D.z,
                                      s.x,s.y,s.z, t1.x,t1.y,t1.z, t2.x,t2.y,t2.z);
        }
    }

    if (nT == 2) {                                   /* ---- TTS vertex (edge) ---- */
        int e[2]; lp3_shared_edge(Tf[0], Tf[1], e);
        s_point A = tet[e[0]], B = tet[e[1]];
        if (d.type == LP3_SEED) {                     /* TTS vs S */
            s_point t = seeds->p[Sc[0].idx], u = seeds->p[d.idx];
            if (weights) {
                double wt = weights[Sc[0].idx], wu = weights[d.idx];
                return power_lp3_feasible_TTS_S(A.x,A.y,A.z, B.x,B.y,B.z, s.x,s.y,s.z,
                                                t.x,t.y,t.z, u.x,u.y,u.z, s_w, wt, wu);
            }
            return lp3_feasible_TTS_S(A.x,A.y,A.z, B.x,B.y,B.z, s.x,s.y,s.z,
                                      t.x,t.y,t.z, u.x,u.y,u.z);
        } else {                                     /* TTS vs T : automatic */
            /* Crossing lies on tet edge (faces Tf[0],Tf[1]); it is on a query
             * face iff that face is one of the two, else strictly interior. */
            return (d.idx == Tf[0] || d.idx == Tf[1]) ? 0 : 1;
        }
    }

    /* nT == 3, TTT vertex == the tet vertex common to faces Tf[0..2]. */
    {
        int V = lp3_common_vertex(Tf[0], Tf[1], Tf[2]);
        if (d.type == LP3_SEED) {                     /* TTT vs S : bisector value at V */
            s_point Vp = tet[V], t = seeds->p[d.idx];
            /* lp_feasible_T0_T1_S(V,s,t) = sign(|V-t|^2 - |V-s|^2) = feasibility
             * (>0 when V is strictly closer to s). */
            if (weights)
                return power_lp2_feasible_T0_T1_S(Vp.x,Vp.y,Vp.z, s.x,s.y,s.z, t.x,t.y,t.z,
                                                  s_w, weights[d.idx]);
            return lp2_feasible_T0_T1_S(Vp.x,Vp.y,Vp.z, s.x,s.y,s.z, t.x,t.y,t.z);
        } else {                                     /* TTT vs T : combinatorial */
            return (d.idx == V) ? 1 : 0;             /* strictly inside only the opposite face */
        }
    }
}


/* Exact one-sided filter: returns 0 only if some bisector s->nbr[j] has every
 * tet vertex strictly on the neighbour's side (cell ^ tet is then empty);
 * 1 = maybe overlaps.  Sound for power bisectors too, since they are affine
 * half-spaces: weights MUST be passed whenever the cell is weighted, or a
 * cell that overlaps under the power metric could be wrongly rejected. */
int lp3_vcell_maybe_hits_tet(const s_point tet[4], s_point s,
                             const s_points *seeds, const double *weights, double s_w,
                             const int *nbr, int NB)
{
    for (int j = 0; j < NB; j++) {
        s_point t = seeds->p[nbr[j]];
        int on_s_side = 0;
        for (int k = 0; k < 4; k++) {
            int fe = weights
                ? power_lp2_feasible_T0_T1_S(tet[k].x,tet[k].y,tet[k].z,
                                             s.x,s.y,s.z, t.x,t.y,t.z, s_w, weights[nbr[j]])
                : lp2_feasible_T0_T1_S(tet[k].x,tet[k].y,tet[k].z,
                                       s.x,s.y,s.z, t.x,t.y,t.z);
            if (fe >= 0) { on_s_side = 1; break; }
        }
        if (!on_s_side) return 0;               /* this bisector cuts the whole tet away */
    }
    return 1;
}


/* Same filter against a triangle (e.g. a boundary face). */
int lp3_vcell_maybe_hits_tri(const s_point tri[3], s_point s,
                             const s_points *seeds, const double *weights, double s_w,
                             const int *nbr, int NB)
{
    for (int j = 0; j < NB; j++) {
        s_point t = seeds->p[nbr[j]];
        int on_s_side = 0;
        for (int k = 0; k < 3; k++) {
            int fe = weights
                ? power_lp2_feasible_T0_T1_S(tri[k].x,tri[k].y,tri[k].z,
                                             s.x,s.y,s.z, t.x,t.y,t.z, s_w, weights[nbr[j]])
                : lp2_feasible_T0_T1_S(tri[k].x,tri[k].y,tri[k].z,
                                       s.x,s.y,s.z, t.x,t.y,t.z);
            if (fe >= 0) { on_s_side = 1; break; }
        }
        if (!on_s_side) return 0;
    }
    return 1;
}
