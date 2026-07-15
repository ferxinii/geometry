/* Degenerate-face repair for closed manifold trimeshes.  See trimesh.h for
 * the contract and the rationale (collapse-only: calibration found no cap
 * faces in the motivating datasets, only blades with microscopic base edges).
 *
 * The repair operates on raw point/face arrays across rounds and is rebuilt
 * through trimesh_from_arrays() once at the end, so the result passes exactly
 * the validation (closed manifold, adjacency, consistent winding, outward
 * orientation) that every constructed trimesh passes.  Manifoldness during
 * the rounds is maintained by the link condition on every collapse. */

#include "trimesh.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <stdbool.h>

/* ---- small helpers ----------------------------------------------------- */

static double dist(s_point a, s_point b)
{
    double dx = a.x - b.x, dy = a.y - b.y, dz = a.z - b.z;
    return sqrt(dx*dx + dy*dy + dz*dz);
}

static s_point face_normal(const s_point *P, const int *f)
{
    s_point u = { .x = P[f[1]].x - P[f[0]].x,
                  .y = P[f[1]].y - P[f[0]].y,
                  .z = P[f[1]].z - P[f[0]].z };
    s_point v = { .x = P[f[2]].x - P[f[0]].x,
                  .y = P[f[2]].y - P[f[0]].y,
                  .z = P[f[2]].z - P[f[0]].z };
    return (s_point){ .x = u.y*v.z - u.z*v.y,
                      .y = u.z*v.x - u.x*v.z,
                      .z = u.x*v.y - u.y*v.x };
}

/* thinness q = 2*area/Lmax^2 and the local index (0..2) of the SHORTEST edge.
 * Edge j connects f[(j+1)%3] and f[(j+2)%3] (opposite vertex j). */
static double face_q(const s_point *P, const int *f, int *shortest_edge)
{
    double L[3];
    for (int j = 0; j < 3; j++)
        L[j] = dist(P[f[(j+1)%3]], P[f[(j+2)%3]]);
    int jmin = 0, jmax = 0;
    for (int j = 1; j < 3; j++) {
        if (L[j] < L[jmin]) jmin = j;
        if (L[j] > L[jmax]) jmax = j;
    }
    if (shortest_edge) *shortest_edge = jmin;
    s_point n = face_normal(P, f);
    double area2 = sqrt(n.x*n.x + n.y*n.y + n.z*n.z);   /* == 2*area */
    return (L[jmax] > 0.0) ? area2 / (L[jmax] * L[jmax]) : 0.0;
}

/* CSR vertex adjacency (unique undirected neighbours), rows sorted. */
typedef struct { int *off, *adj, nv; } s_vadj;

static void vadj_free(s_vadj *A) { free(A->off); free(A->adj); }

static int vadj_build(const int *faces, int nf, int nv, s_vadj *A)
{
    memset(A, 0, sizeof(*A));
    A->nv = nv;
    A->off = calloc((size_t)nv + 1, sizeof(int));
    int *deg = calloc((size_t)nv, sizeof(int));
    if (!A->off || !deg) { free(deg); vadj_free(A); return 0; }

    for (int i = 0; i < nf; i++)
        for (int j = 0; j < 3; j++) {           /* each face edge twice total */
            deg[faces[i*3 + j]] += 2;           /* over the two incident faces, */
        }                                       /* dedup below                  */
    for (int v = 0; v < nv; v++) A->off[v+1] = A->off[v] + deg[v];
    A->adj = malloc((size_t)A->off[nv] * sizeof(int));
    if (!A->adj) { free(deg); vadj_free(A); return 0; }

    memset(deg, 0, (size_t)nv * sizeof(int));
    for (int i = 0; i < nf; i++)
        for (int j = 0; j < 3; j++) {
            int a = faces[i*3 + j];
            int b = faces[i*3 + (j+1)%3];
            A->adj[A->off[a] + deg[a]++] = b;
            A->adj[A->off[b] + deg[b]++] = a;
        }
    free(deg);

    /* sort + dedup each row in place (rows are tiny: insertion sort) */
    for (int v = 0; v < nv; v++) {
        int *row = A->adj + A->off[v];
        int  n   = A->off[v+1] - A->off[v];
        for (int i = 1; i < n; i++) {
            int key = row[i], k = i - 1;
            while (k >= 0 && row[k] > key) { row[k+1] = row[k]; k--; }
            row[k+1] = key;
        }
        int m = 0;
        for (int i = 0; i < n; i++)
            if (i == 0 || row[i] != row[m-1]) row[m++] = row[i];
        /* store the deduped length in the row itself via a sentinel scheme is
         * overkill; keep a second offsets pass instead */
        for (int i = m; i < n; i++) row[i] = -1;   /* padding */
    }
    return 1;
}

/* Is x a neighbour of v?  (row is sorted, -1 padding at the end) */
static bool vadj_has(const s_vadj *A, int v, int x)
{
    for (int i = A->off[v]; i < A->off[v+1]; i++) {
        int y = A->adj[i];
        if (y < 0 || y > x) return false;
        if (y == x) return true;
    }
    return false;
}

/* Link condition for collapsing edge (a,b): the common neighbours of a and b
 * must be EXACTLY the two vertices opposite the edge (c1, c2). */
static bool link_condition_ok(const s_vadj *A, int a, int b, int c1, int c2)
{
    for (int i = A->off[a]; i < A->off[a+1]; i++) {
        int x = A->adj[i];
        if (x < 0) break;
        if (x == b || x == c1 || x == c2) continue;
        if (vadj_has(A, b, x)) return false;   /* extra common neighbour */
    }
    /* c1 and c2 must actually be common neighbours (sanity) */
    return vadj_has(A, a, c1) && vadj_has(A, b, c1) &&
           vadj_has(A, a, c2) && vadj_has(A, b, c2);
}

/* ---- the repair --------------------------------------------------------- */

typedef struct { double q; int a, b; } s_cand;

static int cand_cmp(const void *x, const void *y)
{
    const s_cand *p = (const s_cand *)x, *q = (const s_cand *)y;
    if (p->q != q->q) return (p->q < q->q) ? -1 : 1;
    if (p->a != q->a) return p->a - q->a;
    return p->b - q->b;
}

s_trimesh trimesh_repaired(const s_trimesh *in,
                           const s_trimesh_repair_opts *opts_in,
                           s_trimesh_repair_stats *stats)
{
    s_trimesh_repair_stats st = {0};
    if (stats) *stats = st;
    if (!trimesh_is_valid(in)) return trimesh_NAN;

    const s_trimesh_repair_opts opts =
        opts_in ? *opts_in : TRIMESH_REPAIR_DEFAULTS;

    /* bbox diagonal of the ORIGINAL mesh: fixed scale for all thresholds */
    s_point lo = in->points.p[0], hi = in->points.p[0];
    for (int i = 1; i < in->points.N; i++) {
        s_point p = in->points.p[i];
        if (p.x < lo.x) lo.x = p.x;  if (p.x > hi.x) hi.x = p.x;
        if (p.y < lo.y) lo.y = p.y;  if (p.y > hi.y) hi.y = p.y;
        if (p.z < lo.z) lo.z = p.z;  if (p.z > hi.z) hi.z = p.z;
    }
    const double diag = dist(lo, hi);
    const double max_len = opts.max_collapse_len * diag;

    /* working copies */
    int nv = in->points.N, nf = in->Nf;
    s_point *P = malloc((size_t)nv * sizeof(s_point));
    int     *F = malloc((size_t)nf * 3 * sizeof(int));
    if (!P || !F) { free(P); free(F); return trimesh_NAN; }
    memcpy(P, in->points.p, (size_t)nv * sizeof(s_point));
    memcpy(F, in->faces, (size_t)nf * 3 * sizeof(int));

    s_cand *cands   = NULL;
    char   *dirty   = NULL;
    int     failed  = 0;

    int round = 0;
    for (round = 0; opts.q_min > 0.0 && round < opts.max_rounds; round++) {
        /* candidates on the current arrays */
        free(cands);
        cands = malloc((size_t)nf * sizeof(s_cand));
        if (!cands) { failed = 1; break; }
        int nc = 0;
        for (int i = 0; i < nf; i++) {
            int se;
            double q = face_q(P, &F[i*3], &se);
            if (q >= opts.q_min) continue;
            int a = F[i*3 + (se+1)%3], b = F[i*3 + (se+2)%3];
            if (dist(P[a], P[b]) > max_len) continue;   /* displacement bound */
            if (a > b) { int t = a; a = b; b = t; }
            cands[nc++] = (s_cand){ .q = q, .a = a, .b = b };
        }
        if (nc == 0) break;
        qsort(cands, (size_t)nc, sizeof(s_cand), cand_cmp);

        s_vadj A;
        if (!vadj_build(F, nf, nv, &A)) { failed = 1; break; }
        free(dirty);
        dirty = calloc((size_t)nv, 1);
        if (!dirty) { vadj_free(&A); failed = 1; break; }

        int edits = 0;
        for (int k = 0; k < nc; k++) {
            int a = cands[k].a, b = cands[k].b;
            if (k > 0 && a == cands[k-1].a && b == cands[k-1].b) continue;
            if (dirty[a] || dirty[b]) continue;

            /* the two faces sharing edge (a,b), and their opposite vertices */
            int c[2], nfab = 0;
            bool nbr_dirty = false;
            for (int i = 0; i < nf && nfab <= 2; i++) {
                const int *f = &F[i*3];
                int ia = -1, ib = -1;
                for (int j = 0; j < 3; j++) {
                    if (f[j] == a) ia = j;
                    if (f[j] == b) ib = j;
                }
                if (ia < 0 || ib < 0) continue;
                if (nfab < 2) c[nfab] = f[3 - ia - ib];
                nfab++;
            }
            if (nfab != 2) continue;               /* should not happen        */

            /* independence: no vertex of the 1-rings may be dirty */
            for (int i = A.off[a]; i < A.off[a+1] && !nbr_dirty; i++)
                if (A.adj[i] >= 0 && dirty[A.adj[i]]) nbr_dirty = true;
            for (int i = A.off[b]; i < A.off[b+1] && !nbr_dirty; i++)
                if (A.adj[i] >= 0 && dirty[A.adj[i]]) nbr_dirty = true;
            if (nbr_dirty) continue;

            if (!link_condition_ok(&A, a, b, c[0], c[1])) { st.n_unfixable++; continue; }

            /* flip guard: no surviving incident face may reverse */
            s_point m = { .x = 0.5 * (P[a].x + P[b].x),
                          .y = 0.5 * (P[a].y + P[b].y),
                          .z = 0.5 * (P[a].z + P[b].z) };
            bool flip = false;
            for (int i = 0; i < nf && !flip; i++) {
                const int *f = &F[i*3];
                bool has_a = (f[0]==a || f[1]==a || f[2]==a);
                bool has_b = (f[0]==b || f[1]==b || f[2]==b);
                if (has_a == has_b) continue;      /* untouched, or dies       */
                s_point n0 = face_normal(P, f);
                int g[3];
                for (int j = 0; j < 3; j++) g[j] = (f[j]==a || f[j]==b) ? -1 : f[j];
                s_point pv[3];
                for (int j = 0; j < 3; j++) pv[j] = (g[j] < 0) ? m : P[g[j]];
                s_point n1 = { .x = (pv[1].y-pv[0].y)*(pv[2].z-pv[0].z) - (pv[1].z-pv[0].z)*(pv[2].y-pv[0].y),
                               .y = (pv[1].z-pv[0].z)*(pv[2].x-pv[0].x) - (pv[1].x-pv[0].x)*(pv[2].z-pv[0].z),
                               .z = (pv[1].x-pv[0].x)*(pv[2].y-pv[0].y) - (pv[1].y-pv[0].y)*(pv[2].x-pv[0].x) };
                if (n0.x*n1.x + n0.y*n1.y + n0.z*n1.z <= 0.0) flip = true;
            }
            if (flip) { st.n_unfixable++; continue; }

            /* apply: a <- midpoint, b -> a, drop the two shared faces */
            double disp = 0.5 * dist(P[a], P[b]);
            if (disp > st.max_displacement) st.max_displacement = disp;
            P[a] = m;
            int w = 0;
            for (int i = 0; i < nf; i++) {
                int *f = &F[i*3];
                bool has_a = (f[0]==a || f[1]==a || f[2]==a);
                bool has_b = (f[0]==b || f[1]==b || f[2]==b);
                if (has_a && has_b) continue;      /* the collapsed pair dies  */
                for (int j = 0; j < 3; j++)
                    F[w*3 + j] = (f[j] == b) ? a : f[j];
                w++;
            }
            nf = w;

            /* dirty the whole neighbourhood (adjacency rows are now stale
             * there; remaining candidates touching it wait for next round) */
            dirty[a] = dirty[b] = 1;
            for (int i = A.off[a]; i < A.off[a+1]; i++)
                if (A.adj[i] >= 0) dirty[A.adj[i]] = 1;
            for (int i = A.off[b]; i < A.off[b+1]; i++)
                if (A.adj[i] >= 0) dirty[A.adj[i]] = 1;

            st.n_collapses++;
            edits++;
        }
        vadj_free(&A);

        if (edits == 0) break;                    /* fixed point (or all guarded) */
        st.n_unfixable = 0;                       /* only the LAST round's count  */
    }
    free(cands);
    free(dirty);

    if (failed) { free(P); free(F); return trimesh_NAN; }

    if (round >= opts.max_rounds && opts.q_min > 0.0)
        fprintf(stderr, "trimesh_repaired: round cap (%d) reached; "
                        "%d faces may remain degenerate\n",
                opts.max_rounds, st.n_unfixable);
    st.n_rounds = round;

    /* compact unreferenced vertices, rebuild + revalidate via the constructor */
    int *remap = malloc((size_t)nv * sizeof(int));
    if (!remap) { free(P); free(F); return trimesh_NAN; }
    for (int i = 0; i < nv; i++) remap[i] = -1;
    int nv2 = 0;
    for (int i = 0; i < nf * 3; i++)
        if (remap[F[i]] < 0) remap[F[i]] = nv2++;
    s_point *P2 = malloc((size_t)nv2 * sizeof(s_point));
    if (!P2) { free(remap); free(P); free(F); return trimesh_NAN; }
    for (int i = 0; i < nv; i++)
        if (remap[i] >= 0) P2[remap[i]] = P[i];
    for (int i = 0; i < nf * 3; i++) F[i] = remap[F[i]];
    free(remap); free(P);

    s_trimesh out = trimesh_from_arrays(P2, nv2, F, nf, 0.0);
    free(P2); free(F);
    if (!trimesh_is_valid(&out)) return trimesh_NAN;

    st.volume_drift = fabs(volume_trimesh(&out) - volume_trimesh(in));
    if (stats) *stats = st;
    return out;
}
