
#include "points.h"
#include "gtests.h"
#include "predicates_lp.h"
#include "lp.h"
#include <stdio.h>
#include <assert.h>


int lp2_D(const s_point face[3], s_point s,
          const s_points *seeds, const int *best_idx,
          s_cid2 i, s_cid2 j)
{
    s_point A = face[0], B = face[1], C = face[2];

    /* SEED x SEED */
    if (i.type == LP2_SEED && j.type == LP2_SEED) {
        s_point ti = seeds->p[best_idx[i.idx]];
        s_point tj = seeds->p[best_idx[j.idx]];
        return lp2_det2(A.x,A.y,A.z, B.x,B.y,B.z, C.x,C.y,C.z,
                       s.x,s.y,s.z,
                       ti.x,ti.y,ti.z, tj.x,tj.y,tj.z);
    }

    /* SEED x TRI */
    if (i.type == LP2_SEED && j.type != LP2_SEED) {
        s_point ti = seeds->p[best_idx[i.idx]];
        if (j.type == LP2_TRI0) return -lp2_D_T0_S(A.x,A.y,A.z, C.x,C.y,C.z, s.x,s.y,s.z, ti.x,ti.y,ti.z);
        if (j.type == LP2_TRI1) return -lp2_D_T1_S(A.x,A.y,A.z, B.x,B.y,B.z, s.x,s.y,s.z, ti.x,ti.y,ti.z);
        if (j.type == LP2_TRI2) return -lp2_D_T2_S(B.x,B.y,B.z, C.x,C.y,C.z, s.x,s.y,s.z, ti.x,ti.y,ti.z);
    }

    /* TRI x SEED */
    if (i.type != LP2_SEED && j.type == LP2_SEED) {
        s_point tj = seeds->p[best_idx[j.idx]];
        if (i.type == LP2_TRI0) return lp2_D_T0_S(A.x,A.y,A.z, C.x,C.y,C.z, s.x,s.y,s.z, tj.x,tj.y,tj.z);
        if (i.type == LP2_TRI1) return lp2_D_T1_S(A.x,A.y,A.z, B.x,B.y,B.z, s.x,s.y,s.z, tj.x,tj.y,tj.z);
        if (i.type == LP2_TRI2) return lp2_D_T2_S(B.x,B.y,B.z, C.x,C.y,C.z, s.x,s.y,s.z, tj.x,tj.y,tj.z);
    }

    /* TRI x TRI */
    if (i.type == j.type) return 0;  // Repeated arguments == 0
    else if (i.type == LP2_TRI0 && j.type == LP2_TRI1) return 1;
    else if (i.type == LP2_TRI1 && j.type == LP2_TRI0) return -1;
    else if (i.type == LP2_TRI1 && j.type == LP2_TRI2) return 1;
    else if (i.type == LP2_TRI2 && j.type == LP2_TRI1) return -1;
    else if (i.type == LP2_TRI2 && j.type == LP2_TRI0) return 1;
    else if (i.type == LP2_TRI0 && j.type == LP2_TRI2) return -1;

    assert(0 && "lp2_D: case not implemented");
    return 0;
}


int lp2_feasible(const s_point face[3], s_point s,
                 const s_points *seeds, const int *best_idx,
                 s_cid2 L, s_cid2 M, s_cid2 N)
{
    s_point A = face[0], B = face[1], C = face[2];

    /* duplicate SEED constraints -> degenerate intersection, return 0 */
    if (L.type == LP2_SEED && M.type == LP2_SEED && L.idx == M.idx) return 0;
    if (L.type == LP2_SEED && N.type == LP2_SEED && L.idx == N.idx) return 0;
    if (M.type == LP2_SEED && N.type == LP2_SEED && M.idx == N.idx) return 0;

    /* swap so TRI always precedes SEED in (L,M); sign is unchanged by this swap */
    if (L.type == LP2_SEED && M.type != LP2_SEED) { s_cid2 tmp = L; L = M; M = tmp; }

    /* TRI x TRI x SEED */
    if (L.type != LP2_SEED && M.type != LP2_SEED && N.type == LP2_SEED) {
        s_point ti = seeds->p[best_idx[N.idx]];
        if (L.type == LP2_TRI0 && M.type == LP2_TRI1)
            return lp2_feasible_T0_T1_S(A.x,A.y,A.z, s.x,s.y,s.z, ti.x,ti.y,ti.z);
        if (L.type == LP2_TRI1 && M.type == LP2_TRI2)
            return lp2_feasible_T1_T2_S(B.x,B.y,B.z, s.x,s.y,s.z, ti.x,ti.y,ti.z);
        if (L.type == LP2_TRI0 && M.type == LP2_TRI2)
            return lp2_feasible_T0_T2_S(C.x,C.y,C.z, s.x,s.y,s.z, ti.x,ti.y,ti.z);
    }

    /* TRI x SEED x SEED */
    if (L.type != LP2_SEED && M.type == LP2_SEED && N.type == LP2_SEED) {
        s_point tj = seeds->p[best_idx[M.idx]];
        s_point tl = seeds->p[best_idx[N.idx]];
        if (L.type == LP2_TRI0)
            return lp2_feasible_T0_S_S(A.x,A.y,A.z, B.x,B.y,B.z, C.x,C.y,C.z,
                                      s.x,s.y,s.z, tj.x,tj.y,tj.z, tl.x,tl.y,tl.z);
        if (L.type == LP2_TRI1)
            return lp2_feasible_T1_S_S(A.x,A.y,A.z, B.x,B.y,B.z, C.x,C.y,C.z,
                                      s.x,s.y,s.z, tj.x,tj.y,tj.z, tl.x,tl.y,tl.z);
        if (L.type == LP2_TRI2)
            return lp2_feasible_T2_S_S(A.x,A.y,A.z, B.x,B.y,B.z, C.x,C.y,C.z,
                                      s.x,s.y,s.z, tj.x,tj.y,tj.z, tl.x,tl.y,tl.z);
    }

    /* TRI x SEED x TRI */
    if (L.type != LP2_SEED && M.type == LP2_SEED && N.type != LP2_SEED) {
        if (N.type == L.type) return 0;
        s_point tj = seeds->p[best_idx[M.idx]];
        if (L.type == LP2_TRI0 && N.type == LP2_TRI1)
            return lp2_feasible_T0_S_T1(A.x,A.y,A.z, C.x,C.y,C.z, s.x,s.y,s.z, tj.x,tj.y,tj.z);
        if (L.type == LP2_TRI0 && N.type == LP2_TRI2)
            return lp2_feasible_T0_S_T2(A.x,A.y,A.z, C.x,C.y,C.z, s.x,s.y,s.z, tj.x,tj.y,tj.z);
        if (L.type == LP2_TRI1 && N.type == LP2_TRI0)
            return lp2_feasible_T1_S_T0(A.x,A.y,A.z, B.x,B.y,B.z, s.x,s.y,s.z, tj.x,tj.y,tj.z);
        if (L.type == LP2_TRI1 && N.type == LP2_TRI2)
            return lp2_feasible_T1_S_T2(A.x,A.y,A.z, B.x,B.y,B.z, s.x,s.y,s.z, tj.x,tj.y,tj.z);
        if (L.type == LP2_TRI2 && N.type == LP2_TRI0)
            return lp2_feasible_T2_S_T0(B.x,B.y,B.z, C.x,C.y,C.z, s.x,s.y,s.z, tj.x,tj.y,tj.z);
        if (L.type == LP2_TRI2 && N.type == LP2_TRI1)
            return lp2_feasible_T2_S_T1(B.x,B.y,B.z, C.x,C.y,C.z, s.x,s.y,s.z, tj.x,tj.y,tj.z);
    }

    /* SEED x SEED x TRI */
    if (L.type == LP2_SEED && M.type == LP2_SEED && N.type != LP2_SEED) {
        s_point tj = seeds->p[best_idx[L.idx]];
        s_point tk = seeds->p[best_idx[M.idx]];
        if (N.type == LP2_TRI0)
            return lp2_feasible_S_S_T0(A.x,A.y,A.z, B.x,B.y,B.z, C.x,C.y,C.z,
                                      s.x,s.y,s.z, tj.x,tj.y,tj.z, tk.x,tk.y,tk.z);
        if (N.type == LP2_TRI1)
            return lp2_feasible_S_S_T1(A.x,A.y,A.z, B.x,B.y,B.z, C.x,C.y,C.z,
                                      s.x,s.y,s.z, tj.x,tj.y,tj.z, tk.x,tk.y,tk.z);
        if (N.type == LP2_TRI2)
            return lp2_feasible_S_S_T2(A.x,A.y,A.z, B.x,B.y,B.z, C.x,C.y,C.z,
                                      s.x,s.y,s.z, tj.x,tj.y,tj.z, tk.x,tk.y,tk.z);
    }

    /* TRI x TRI x TRI: intersection is a vertex, always inside unless any two match */
    if (L.type != LP2_SEED && M.type != LP2_SEED && N.type != LP2_SEED)
        return (L.type == M.type || N.type == L.type || N.type == M.type) ? 0 : 1;

    /* SEED x SEED x SEED */
    if (L.type == LP2_SEED && M.type == LP2_SEED && N.type == LP2_SEED) {
        s_point tL = seeds->p[best_idx[L.idx]];
        s_point tM = seeds->p[best_idx[M.idx]];
        s_point tN = seeds->p[best_idx[N.idx]];
        return lp2_det2(A.x,A.y,A.z, B.x,B.y,B.z, C.x,C.y,C.z,
                       s.x,s.y,s.z, tL.x,tL.y,tL.z, tM.x,tM.y,tM.z)
             * lp2_det3(A.x,A.y,A.z, B.x,B.y,B.z, C.x,C.y,C.z,
                       s.x,s.y,s.z, tL.x,tL.y,tL.z, tM.x,tM.y,tM.z,
                       tN.x,tN.y,tN.z);
    }

    fprintf(stderr, "lp2_feasible: unhandled constraint combination\n");
    return -1;
}

/* Update one 1D-LP bound. Returns false if the interval becomes empty or
 * the line is globally infeasible. When concurrent=true, a zero feasibility
 * result while tightening a non-null slot signals a concurrent triangle. */
static bool lp2_process_cj(const s_point face[3], s_point s,
                           const s_points *seeds, const int *seeds_idx,
                           s_cid2 ci, s_cid2 cj,
                           s_cid2 *lo, bool *lo_null,
                           s_cid2 *hi, bool *hi_null)
{
    int d = lp2_D(face, s, seeds, seeds_idx, ci, cj);
    if (d == 0) {
        /* C_j parallel to C_i: use whichever bound is non-null as reference
         * (same sign for any point on C_i when C_j is parallel). */
        s_cid2 ref = *lo_null ? *hi : *lo;
        int v = lp2_feasible(face, s, seeds, seeds_idx, ci, ref, cj);
        if (v < 0) return false;
    } else if (d > 0) {
        if (*lo_null) {
            *lo = cj; *lo_null = false;
        } else {
            int v = lp2_feasible(face, s, seeds, seeds_idx, ci, cj, *lo);
            if (v > 0) *lo = cj;
        }
    } else {
        if (*hi_null) {
            *hi = cj; *hi_null = false;
        } else {
            int v = lp2_feasible(face, s, seeds, seeds_idx, ci, cj, *hi);
            if (v > 0) *hi = cj;
        }
    }
    if (!*lo_null && !*hi_null) {
        int empty = lp2_feasible(face, s, seeds, seeds_idx, ci, *lo, *hi);
        if (empty < 0) return false;
    }
    return true;
}


bool lp2_vcell_hits_face(const s_point face[3], const s_points *seeds,
                         int N, int seeds_idx[N], int id,
                          int (*randint)(void *rctx, int), void *rctx)
{
    s_point s = seeds->p[seeds_idx[id]];
    s_cid2 T0 = {LP2_TRI0, -1}, T1 = {LP2_TRI1, -1}, T2 = {LP2_TRI2, -1};

    /* Step 0: triangle degeneracy. Collinear vertices give orient3d = 0 for
     * every fourth point; at most two axis-offset candidates fit any plane,
     * so all three returning 0 is the exact collinearity condition. */
    if (test_points_colinear(face[0], face[1], face[2])) return false;

    /* Step 1: initial candidate = T0 ^ T1 (a triangle vertex). */
    s_cid2 A = T0, B = T1;

    /* Fisher-Yates shuffle. */
    int order[N];
    for (int i = 0; i < N; i++) order[i] = i;
    for (int i = N-1; i > 0; i--) {
        int j = randint(rctx, i+1);
        int tmp = order[i]; order[i] = order[j]; order[j] = tmp;
    }

    /* gate: seed 123 (x~0.833) on the x=1 face (all verts have x=1) */
    for (int oi = 0; oi < N; oi++) {
        int i = order[oi];
        if (i == id) continue;
        s_cid2 ci = {LP2_SEED, i};

        int outer_feas = lp2_feasible(face, s, seeds, seeds_idx, A, B, ci);
        if (outer_feas >= 0) continue;

        /* 1D LP on C_i's line. lo/hi track interval endpoints; *_null = true
         * means that slot is unbounded (no constraint has closed it yet). */
        s_cid2 lo = T0, hi = T0; /* dummy init to silence uninit warnings */
        bool lo_null, hi_null;

        /* Step 3a: seed the interval from C_T0. */
        int d_i_T0 = lp2_D(face, s, seeds, seeds_idx, ci, T0);
        if (d_i_T0 != 0) {
            /* Generic case: one slot from C_T0, other stays open. */
            if (d_i_T0 > 0) { lo = T0; lo_null = false; hi_null = true; }
            else             { hi = T0; hi_null = false; lo_null = true; }
            if (!lp2_process_cj(face, s, seeds, seeds_idx, ci, T1,
                                &lo, &lo_null, &hi, &hi_null))
                return false;
        } else {
            /* Degenerate case: C_T0 is parallel to C_i.
             * Check whole line feasibility w.r.t. C_T0, then seed from C_T1. */
            int par_check = lp2_feasible(face, s, seeds, seeds_idx, ci, T1, T0);
            if (par_check < 0) return false;
            int d_i_T1 = lp2_D(face, s, seeds, seeds_idx, ci, T1);
            if (d_i_T1 > 0) { lo = T1; lo_null = false; hi_null = true; }
            else             { hi = T1; hi_null = false; lo_null = true; }
        }

        /* Step 3c: C_T2, then all SEED constraints before C_i in random order. */
        if (!lp2_process_cj(face, s, seeds, seeds_idx, ci, T2,
                            &lo, &lo_null, &hi, &hi_null))
            return false;

        for (int oj = 0; oj < oi; oj++) {
            int jj = order[oj];
            if (jj == id) continue;
            s_cid2 cj = {LP2_SEED, jj};
            if (!lp2_process_cj(face, s, seeds, seeds_idx, ci, cj,
                                &lo, &lo_null, &hi, &hi_null))
                return false;
        }

        /* Step 4: new candidate is C_i ^ lo. */
        A = ci;
        B = lo;
    }

    return true;
}


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
 * the tet interior is the feasible side.                                      */
int lp3_feasible(const s_point tet[4], const int sigma_f[4], s_point s,
                 const s_points *seeds,
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
            return -test_insphere(sph, u);
        } else {                                     /* SSS vs T */
            int fv[3]; lp3_face_idx(d.idx, fv);
            s_point A = tet[fv[0]], Q = tet[fv[1]], R = tet[fv[2]];
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
            return lp2_feasible_T0_T1_S(Vp.x,Vp.y,Vp.z, s.x,s.y,s.z, t.x,t.y,t.z);
        } else {                                     /* TTT vs T : combinatorial */
            return (d.idx == V) ? 1 : 0;             /* strictly inside only the opposite face */
        }
    }
}


int lp3_vcell_maybe_hits_tet(const s_point tet[4], s_point s,
                             const s_points *seeds,
                             const int *nbr, int NB)
{
    for (int j = 0; j < NB; j++) {
        s_point t = seeds->p[nbr[j]];
        int on_s_side = 0;
        for (int k = 0; k < 4; k++)
            if (lp2_feasible_T0_T1_S(tet[k].x,tet[k].y,tet[k].z,
                                    s.x,s.y,s.z, t.x,t.y,t.z) >= 0) { on_s_side = 1; break; }
        if (!on_s_side) return 0;               /* this bisector cuts the whole tet away */
    }
    return 1;
}


int lp3_vcell_maybe_hits_tri(const s_point tri[3], s_point s,
                             const s_points *seeds,
                             const int *nbr, int NB)
{
    for (int j = 0; j < NB; j++) {
        s_point t = seeds->p[nbr[j]];
        int on_s_side = 0;
        for (int k = 0; k < 3; k++)
            if (lp2_feasible_T0_T1_S(tri[k].x,tri[k].y,tri[k].z,
                                     s.x,s.y,s.z, t.x,t.y,t.z) >= 0) { on_s_side = 1; break; }
        if (!on_s_side) return 0;
    }
    return 1;
}


