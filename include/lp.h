#ifndef _GEOMETRY_LP_H
#define _GEOMETRY_LP_H

#include <assert.h>
#include "points.h"
#include "gtests.h"

#ifdef GEOMETRY_COMPILE_LP

/* Identifies a constraint.  TRI* idx is unused (set to -1).
 * CSEED idx is the constraint's position in the best_idx[] array. */
typedef enum { LP2_TRI0, LP2_TRI1, LP2_TRI2, LP2_SEED } e_ctype2;
typedef struct { e_ctype2 type; int idx; } s_cid2;

int lp2_D(const s_point face[3], s_point s,
          const s_points *seeds, const int *best_idx,
          s_cid2 i, s_cid2 j);
int lp2_feasible(const s_point face[3], s_point s,
                 const s_points *seeds, const int *best_idx,
                 s_cid2 L, s_cid2 M, s_cid2 N);
bool lp2_vcell_hits_face(const s_point face[3], const s_points *seeds,
                         int N, int seeds_idx[N], int id,
                         int (*randint)(void *rctx, int), void *rctx);


/* Constraint identifier. CTET: idx = face 0..3 (face f = the 3 verts != v[f]).
 * CSEED: idx = partner seed index into seeds->p[]. */
typedef enum { LP3_TET, LP3_SEED } e_ctype3;
typedef struct { e_ctype3 type; int idx; } s_cid3;


int lp3_D(const s_point tet[4], s_point s, const s_points *seeds,
          s_cid3 a, s_cid3 b, s_cid3 c);
int lp3_feasible(const s_point tet[4], const int sigma_f[4], s_point s,
                 const s_points *seeds,
                 s_cid3 a, s_cid3 b, s_cid3 c, s_cid3 d);
int lp3_vcell_maybe_hits_tet(const s_point tet[4], s_point s,
                            const s_points *seeds,
                            const int *nbr, int nb);
int lp3_vcell_maybe_hits_tri(const s_point tri[3], s_point s,
                            const s_points *seeds,
                            const int *nbr, int NB);

// ***************************** HELPERS ***************************************
static inline void lp3_face_idx(int f, int out[3])
{   /* The three vertex indices of face f (the ones != f), ascending. */
    int k = 0;
    for (int i = 0; i < 4; i++) if (i != f) out[k++] = i;
}

static inline void lp3_shared_edge(int f, int g, int edge[2])
{  /* The two shared-edge vertex indices of faces f and g (the ones != f and != g). */
    int k = 0;
    for (int i = 0; i < 4; i++) if (i != f && i != g) edge[k++] = i;
}

static inline int lp3_common_vertex(int f, int g, int h)
{   /* The single vertex common to faces f,g,h == the remaining index (0+1+2+3=6). */
    return 6 - f - g - h;
}

static inline void lp3_tet_face_orient(const s_point tet[4], int sigma_f[4])
{  /* Per-tet face interior-orientation signs. sigma_f[f] = orient3d(face f, v[f]) */
    for (int f = 0; f < 4; f++) {
        int v[3]; lp3_face_idx(f, v);
        s_point A = tet[v[0]], B = tet[v[1]], C = tet[v[2]], W = tet[f];
        sigma_f[f] = test_orientation((s_point[]){A, B, C}, W);
    }
}

static inline void lp3_split3(s_cid3 a, s_cid3 b, s_cid3 c,
                       int Tf[3], int *nT, s_cid3 Sc[3], int *nS)
{  /* Split (a,b,c) into the tet-face ids (T) and seed cids (S). */
    s_cid3 in[3] = { a, b, c };
    *nT = 0; *nS = 0;
    for (int i = 0; i < 3; i++) {
        if (in[i].type == LP3_TET) Tf[(*nT)++] = in[i].idx;
        else                       Sc[(*nS)++] = in[i];
    }
}

#else  // GEOMETRY_COMPILE_LP
static_assert(0, "Need to turn ON CMake option 'GEOMETRY_COMPILE_LP' (or flag GEOMETRY_COMPILE_LP)");
#endif  // GEOMETRY_COMPILE_LP


#endif
