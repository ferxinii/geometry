#include "points.h"
#include "gtests.h"
#include "dynarray.h"
#include "convh.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include <assert.h>
#include <stdbool.h>


/* ===== Geometry helpers ===== */

static bool set_true_if_not_already(bool *arr, int id)
{
    if (!arr[id]) { arr[id] = true; return 1; }
    return 0;
}

static void vertices_face(const s_points *points, const int *faces, int face_id, s_point out[3])
{
    out[0] = points->p[faces[face_id*3+0]];
    out[1] = points->p[faces[face_id*3+1]];
    out[2] = points->p[faces[face_id*3+2]];
}

static void flip_face(int *faces, int face_id)
{
    int tmp = faces[face_id*3+1];
    faces[face_id*3+1] = faces[face_id*3+2];
    faces[face_id*3+2] = tmp;
}

static int next_vid_isused_notinface(int Np, bool isused[Np], int face_ids[3], int current_id)
{
    for (int ii = current_id + 1; ii < Np; ++ii)
        if (isused[ii] && ii != face_ids[0] && ii != face_ids[1] && ii != face_ids[2])
            return ii;
    return -1;
}

static int orient_face_if_needed(const s_points *points, bool isused[points->N],
                                  int Nfaces, int faces[3*Nfaces], int face_id)
{
    int face_vids[3] = {faces[face_id*3+0], faces[face_id*3+1], faces[face_id*3+2]};
    s_point face_vertices[3];
    vertices_face(points, faces, face_id, face_vertices);

    int p = -1;
    p = next_vid_isused_notinface(points->N, isused, face_vids, p);
    if (p == -1) return -1;

    int o = test_orientation(face_vertices, points->p[p]);
    while (o == 0) {
        p = next_vid_isused_notinface(points->N, isused, face_vids, p);
        if (p == -1) return -1;
        o = test_orientation(face_vertices, points->p[p]);
    }

    if (o == -1) { flip_face(faces, face_id); return 1; }
    return 0;
}

/* ===== Initial tetrahedron ===== */

/* Pick 4 non-coplanar seed vertices for the initial tetrahedron.
 * Strategy (from Quickhull): find the 6 AABB extreme points (one per
 * min/max face of the bounding box), pick the most-separated pair (v0,v1),
 * then the point farthest from that edge (v2), then the point farthest
 * from that triangle (v3). */
static int initial_tetra_vids(const s_points *points, double EPS_degenerate, int out[4])
{
    int N = points->N;
    assert(N >= 4);

    /* Phase 1: collect the 6 axis-aligned extreme-point indices (min/max for x, y, z).
     * cands[0,2,4] = indices of min-x, min-y, min-z points.
     * cands[1,3,5] = indices of max-x, max-y, max-z points. */
    int     cands[6] = {-1,-1,-1,-1,-1,-1};
    s_point bb_min = {.x =  DBL_MAX, .y =  DBL_MAX, .z =  DBL_MAX};
    s_point bb_max = {.x = -DBL_MAX, .y = -DBL_MAX, .z = -DBL_MAX};
    for (int i = 0; i < N; i++) {
        s_point pt = points->p[i];
        if (pt.x < bb_min.x) { bb_min.x = pt.x; cands[0] = i; }
        if (pt.x > bb_max.x) { bb_max.x = pt.x; cands[1] = i; }
        if (pt.y < bb_min.y) { bb_min.y = pt.y; cands[2] = i; }
        if (pt.y > bb_max.y) { bb_max.y = pt.y; cands[3] = i; }
        if (pt.z < bb_min.z) { bb_min.z = pt.z; cands[4] = i; }
        if (pt.z > bb_max.z) { bb_max.z = pt.z; cands[5] = i; }
    }

    /* Phase 2: v0, v1  -- the most-separated pair among the 6 candidates.
     * Skip same-index pairs (one point can be extreme on multiple axes). */
    int v0 = -1, v1 = -1;
    double best = -1.0;
    for (int a = 0; a < 6; a++) {
        for (int b = a+1; b < 6; b++) {
            if (cands[b] == cands[a]) continue;
            double d = distance_squared(points->p[cands[a]], points->p[cands[b]]);
            if (d > best) { best = d; v0 = cands[a]; v1 = cands[b]; }
        }
    }
    if (v0 < 0) { fprintf(stderr, "ch_quickhull3D: all points coincide.\n"); return 0; }

    /* Phase 3: v2  -- farthest point from the line (v0, v1). */
    int v2 = -1;
    double best2 = -1.0;
    s_point line[2] = {points->p[v0], points->p[v1]};
    for (int i = 0; i < N; i++) {
        if (i == v0 || i == v1) continue;
        double d = distance_sq_point_line(line, EPS_degenerate, points->p[i]);
        if (d > best2) { best2 = d; v2 = i; }
    }
    if (v2 < 0) { fprintf(stderr, "ch_quickhull3D: all points are collinear.\n"); return 0; }

    /* Phase 4: v3  -- farthest point from the plane (v0, v1, v2).
     * If signed_distance_point_to_plane returns NaN (degenerate triangle  --
     * v0,v1,v2 collinear), no v3 is found and v3<0 triggers the error below. */
    s_point tri[3] = {points->p[v0], points->p[v1], points->p[v2]};
    int v3 = -1;
    double best3 = -1.0;
    for (int i = 0; i < N; i++) {
        if (i == v0 || i == v1 || i == v2) continue;
        double d = fabs(signed_distance_point_to_plane(points->p[i], tri, EPS_degenerate));
        if (d > best3) { best3 = d; v3 = i; }
    }

    /* Robust coplanarity check  -- catches collinear inputs too. */
    if (v3 < 0 || test_orientation(tri, points->p[v3]) == 0) {
        fprintf(stderr, "ch_quickhull3D: all points are coplanar.\n");
        return 0;
    }

    out[0] = v0; out[1] = v1; out[2] = v2; out[3] = v3;
    return 1;
}

static int initial_tetrahedron(const s_points *points, bool *isused,
                                double EPS_degenerate, int faces[4*3], int adj[4*3])
{
    int vertex_ids[4];
    if (!initial_tetra_vids(points, EPS_degenerate, vertex_ids)) return 0;

    memset(isused, 0, points->N * sizeof(bool));
    isused[vertex_ids[0]] = isused[vertex_ids[1]] = isused[vertex_ids[2]] = isused[vertex_ids[3]] = true;

    for (int i=0; i<4; i++)
        for (int j=0, k=0; j<4; j++)
            if (vertex_ids[j] != vertex_ids[i]) faces[i*3+k++] = vertex_ids[j];

    for (int k=0; k<4; k++)
        if (orient_face_if_needed(points, isused, 4, faces, k) == -1) {
            fprintf(stderr, "ch_quickhull3D: initial tetrahedron orientation failed.\n");
            return 0;
        }

    /* adj[i*3+k] = the face that does NOT contain vertex faces[i*3+k] */
    for (int i = 0; i < 4; i++)
        for (int k = 0; k < 3; k++) {
            int v = faces[i*3+k];
            for (int j = 0; j < 4; j++) {
                if (j == i) continue;
                int *fj = faces + j*3;
                if (fj[0] != v && fj[1] != v && fj[2] != v) { adj[i*3+k] = j; break; }
            }
        }

    return 1;
}


/* ===== Visibility ===== */

static int coplanar_visible(const int *adj,
                             s_dynarray *cop_fids,
                             s_dynarray *out_indicator)
{
    int N_visible = 0;

    /* Propagate visibility from anchor-visible faces (test_orientation < 0)
     * into the coplanar group via adjacency.  A coplanar face becomes visible
     * if any of its adj neighbours is already visible (anchor or coplanar). */
    bool changed = true;
    while (changed) {
        changed = false;
        for (int ii = 0; ii < (int)cop_fids->N; ii++) {
            int fid; dynarray_get_value(cop_fids, ii, &fid);
            if (((bool *)out_indicator->items)[fid]) continue;
            for (int k = 0; k < 3; k++) {
                int nbr = adj[fid*3+k];
                if (nbr >= 0 && ((bool *)out_indicator->items)[nbr]) {
                    if (set_true_if_not_already(out_indicator->items, fid)) {
                        N_visible++;
                        changed = true;
                    }
                    break;
                }
            }
        }
    }

    return N_visible;
}


static int visible_faces_from_point(const s_points *points, int Nfaces, int faces[Nfaces*3],
                                     const int *adj,
                                     s_point p, double EPS,
                                     s_dynarray *coplanar_fids,
                                     s_dynarray *out_indicator)
{
    if (!dynarray_ensure_capacity(out_indicator, Nfaces)) {
        fprintf(stderr, "ch_quickhull3D: dynarray error in visible_faces_from_point.\n");
        return -1;
    }
    dynarray_memset0(out_indicator);
    out_indicator->N = Nfaces;
    coplanar_fids->N = 0;
    int N_visible = 0;

    for (int j = 0; j < Nfaces; ++j) {
        s_point face_pts[3]; vertices_face(points, faces, j, face_pts);
        int o = test_orientation(face_pts, p);
        if (o < 0) {
            if (set_true_if_not_already(out_indicator->items, j)) N_visible++;
        } else if (o == 0) {
            /* Exact duplicate of a hull vertex: already on the hull, skip. */
            int vi0=faces[j*3+0], vi1=faces[j*3+1], vi2=faces[j*3+2];
            if ((p.x==points->p[vi0].x && p.y==points->p[vi0].y && p.z==points->p[vi0].z) ||
                (p.x==points->p[vi1].x && p.y==points->p[vi1].y && p.z==points->p[vi1].z) ||
                (p.x==points->p[vi2].x && p.y==points->p[vi2].y && p.z==points->p[vi2].z))
                return 0;
            e_geom_test test = test_point_in_triangle_3D(face_pts, p, EPS, 0.0);
            if (test == TEST_BOUNDARY || test == TEST_IN)
                return 0;
            if ((test == TEST_OUT || test == TEST_DEGENERATE) &&
                !dynarray_push(coplanar_fids, &j)) return 0;
        }
    }

    if (coplanar_fids->N > 0)
        N_visible += coplanar_visible(adj, coplanar_fids, out_indicator);

    return N_visible;
}


/* ===== Horizon ===== */

static int extract_horizon_dfs(const int *faces, const int *adj,
                                 bool *is_visible, bool *dfs_visited, int *dfs_stack,
                                 int start_face, int *out_Nhorizon, int *out_N_vf,
                                 s_dynarray *horizon, s_dynarray *horizon_nbr)
{
    horizon->N = 0;
    horizon_nbr->N = 0;
    if (start_face < 0) { *out_Nhorizon = 0; return 1; }

    int top = 0;
    dfs_stack[top++] = start_face;
    dfs_visited[start_face] = true;

    while (top > 0) {
        int f = dfs_stack[--top];
        for (int k = 0; k < 3; k++) {
            int nbr = adj[f*3+k];
            if (nbr < 0) continue;

            if (is_visible[nbr]) {
                if (!dfs_visited[nbr]) { dfs_visited[nbr] = true; dfs_stack[top++] = nbr; }
                continue;
            }

            /* Non-visible neighbour: candidate horizon edge (canonical vertex order) */
            int vA = faces[f*3+(k+1)%3], vB = faces[f*3+(k+2)%3];
            int ca = vA < vB ? vA : vB, cb = vA < vB ? vB : vA;

            bool dup = false;
            int dup_hi = -1;
            int *h = horizon->items;
            for (int hi = 0; hi < (int)horizon->N; hi += 2)
                if (h[hi] == ca && h[hi+1] == cb) { dup = true; dup_hi = hi; break; }

            if (dup) {
                /* Sandwiched face: borders two visible arcs  -- mark it visible */
                is_visible[nbr] = true;
                (*out_N_vf)++;
                if (!dfs_visited[nbr]) { dfs_visited[nbr] = true; dfs_stack[top++] = nbr; }
                /* Remove duplicate horizon entry by swapping with last */
                int last2 = (int)horizon->N - 2;
                h[dup_hi] = h[last2]; h[dup_hi+1] = h[last2+1];
                ((int*)horizon_nbr->items)[dup_hi/2] = ((int*)horizon_nbr->items)[horizon_nbr->N-1];
                horizon->N -= 2; horizon_nbr->N -= 1;
            } else {
                if (!dynarray_push(horizon, &ca)) return 0;
                if (!dynarray_push(horizon, &cb)) return 0;
                if (!dynarray_push(horizon_nbr, &nbr)) return 0;
            }
        }
    }

    *out_Nhorizon = (int)(horizon->N / 2);
    return 1;
}


static void sort_horizon_loop(int *horizon, int *horizon_nbr, int Nhorizon)
{
    for (int j = 0; j < Nhorizon - 1; j++) {
        int tail = horizon[j*2+1];
        for (int k = j+1; k < Nhorizon; k++) {
            if (horizon[k*2+0] != tail && horizon[k*2+1] != tail) continue;
            if (horizon[k*2+1] == tail) {
                int t = horizon[k*2+0]; horizon[k*2+0] = horizon[k*2+1]; horizon[k*2+1] = t;
            }
            if (k != j+1) {
                int te;
                te = horizon[(j+1)*2+0]; horizon[(j+1)*2+0] = horizon[k*2+0]; horizon[k*2+0] = te;
                te = horizon[(j+1)*2+1]; horizon[(j+1)*2+1] = horizon[k*2+1]; horizon[k*2+1] = te;
                int tn = horizon_nbr[j+1]; horizon_nbr[j+1] = horizon_nbr[k]; horizon_nbr[k] = tn;
            }
            break;
        }
    }
}


/* ===== Face management ===== */

static void compact_faces_adj(int Nfaces, int *faces, int *adj, const int *buff_pmap)
{
    for (int j = 0; j < Nfaces; j++) {
        if (buff_pmap[j] < 0) continue;
        int nj = buff_pmap[j];
        faces[nj*3+0] = faces[j*3+0]; faces[nj*3+1] = faces[j*3+1]; faces[nj*3+2] = faces[j*3+2];
        adj[nj*3+0] = adj[j*3+0] >= 0 ? buff_pmap[adj[j*3+0]] : -1;
        adj[nj*3+1] = adj[j*3+1] >= 0 ? buff_pmap[adj[j*3+1]] : -1;
        adj[nj*3+2] = adj[j*3+2] >= 0 ? buff_pmap[adj[j*3+2]] : -1;
    }
}

static void delete_visible_faces(int Nfaces, int *faces, int *adj,
                                   const bool *faces_isvisible, bool *isused,
                                   int *buff_pmap)
{
    int l = 0;
    for (int j = 0; j < Nfaces; j++)
        buff_pmap[j] = faces_isvisible[j] ? -1 : l++;

    for (int j = 0; j < Nfaces; j++)
        if (faces_isvisible[j])
            isused[faces[j*3+0]] = isused[faces[j*3+1]] = isused[faces[j*3+2]] = false;
    compact_faces_adj(Nfaces, faces, adj, buff_pmap);
}

static int add_faces_from_horizon(const s_points *points, bool *isused,
                                   int Nfaces, s_dynarray *faces, s_dynarray *adj,
                                   int Nhorizon, const int *horizon, const int *horizon_nbr,
                                   int query_pid)
{
    int N_realloc = Nfaces + Nhorizon;
    if (N_realloc >= CH_MAX_NUM_FACES) {
        fprintf(stderr, "ch_quickhull3D: max faces exceeded (%d + %d).\n", Nfaces, Nhorizon);
        return 0;
    }
    faces->N = adj->N = N_realloc * 3;
    if (!dynarray_ensure_capacity(faces, N_realloc*3) ||
        !dynarray_ensure_capacity(adj,   N_realloc*3)) {
        fprintf(stderr, "ch_quickhull3D: dynarray error in add_faces_from_horizon.\n");
        return 0;
    }

    int start = Nfaces;
    int *f = faces->items, *a = adj->items;

    /* Phase A: create new faces, init adj to -1 */
    for (int j = 0; j < Nhorizon; j++) {
        f[Nfaces*3+0] = horizon[j*2+0]; f[Nfaces*3+1] = horizon[j*2+1]; f[Nfaces*3+2] = query_pid;
        a[Nfaces*3+0] = a[Nfaces*3+1] = a[Nfaces*3+2] = -1;
        isused[horizon[j*2+0]] = isused[horizon[j*2+1]] = true;
        Nfaces++;
    }

    /* Phase B: orient new faces */
    for (int k = start; k < Nfaces; k++) {
#ifndef NDEBUG
        int o = orient_face_if_needed(points, isused, Nfaces, f, k);
        assert(o != -1 && "ch_quickhull3D: face orientation failed in add_faces_from_horizon.");
#else
        orient_face_if_needed(points, isused, Nfaces, f, k);
#endif
    }

    /* Phase C: wire each new face to its old non-visible neighbour. */
    for (int j = 0; j < Nhorizon; j++) {
        int fi  = start + j;
        int nbr = horizon_nbr[j];
        int hA  = horizon[j*2+0], hB = horizon[j*2+1];

        assert(nbr >= 0);
        assert((f[nbr*3+0]==hA||f[nbr*3+1]==hA||f[nbr*3+2]==hA) &&
               (f[nbr*3+0]==hB||f[nbr*3+1]==hB||f[nbr*3+2]==hB));

        int k_q = -1;
        for (int k = 0; k < 3; k++) if (f[fi*3+k] == query_pid) { k_q = k; break; }
        assert(k_q >= 0);
        a[fi*3+k_q] = nbr;

        for (int k = 0; k < 3; k++) {
            if (a[nbr*3+k] != -1) continue;
            int ea = f[nbr*3+(k+1)%3], eb = f[nbr*3+(k+2)%3];
            if ((ea==hA&&eb==hB)||(ea==hB&&eb==hA)) { a[nbr*3+k] = fi; break; }
        }
    }

    /* Phase D: wire adj between consecutive new faces along the horizon loop. */
    for (int j = 0; j < Nhorizon; j++) {
        int fi  = start + j;
        int fi1 = start + (j + 1) % Nhorizon;
        if (a[fi*3+0] != -1 && a[fi*3+1] != -1 && a[fi*3+2] != -1) continue;
        for (int ka = 0; ka < 3; ka++) {
            if (f[fi*3+ka] == query_pid || a[fi*3+ka] != -1) continue;
            int ea = f[fi*3+(ka+1)%3], eb = f[fi*3+(ka+2)%3];
            for (int kb = 0; kb < 3; kb++) {
                if (f[fi1*3+kb] == query_pid || a[fi1*3+kb] != -1) continue;
                int ea2 = f[fi1*3+(kb+1)%3], eb2 = f[fi1*3+(kb+2)%3];
                if ((ea==ea2&&eb==eb2)||(ea==eb2&&eb==ea2)) {
                    a[fi*3+ka] = fi1; a[fi1*3+kb] = fi; goto NEXT_J;
                }
            }
        }
        NEXT_J:;
    }

    return 1;
}


#ifndef NDEBUG
// The following are only compiled / used if in debug mode
static int edge_pair_cmp(const void *pa, const void *pb)
{
    const int *a = pa, *b = pb;
    if (a[0] != b[0]) return a[0] < b[0] ? -1 : 1;
    if (a[1] != b[1]) return a[1] < b[1] ? -1 : 1;
    return 0;
}

static int topology_is_valid(int Nfaces, int *faces, s_dynarray *out_edges)
{
    out_edges->N = 0;
    for (int f=0; f<Nfaces; f++) {
        int a0=faces[f*3+0], a1=faces[f*3+1], a2=faces[f*3+2];
        if (!dynarray_push(out_edges, &(int){a0<a1?a0:a1})) return 0;
        if (!dynarray_push(out_edges, &(int){a0<a1?a1:a0})) return 0;
        if (!dynarray_push(out_edges, &(int){a1<a2?a1:a2})) return 0;
        if (!dynarray_push(out_edges, &(int){a1<a2?a2:a1})) return 0;
        if (!dynarray_push(out_edges, &(int){a2<a0?a2:a0})) return 0;
        if (!dynarray_push(out_edges, &(int){a2<a0?a0:a2})) return 0;
    }
    qsort(out_edges->items, Nfaces*3, sizeof(int)*2, edge_pair_cmp);

    int i = 0, seen_boundary = 0;
    while (i < Nfaces*3) {
        int *start = (int*)dynarray_get_ptr(out_edges, 2*i);
        int run_end = i + 1;
        while (run_end < Nfaces*3) {
            int *curr = (int*)dynarray_get_ptr(out_edges, 2*run_end);
            if (!curr) return -1;
            if (curr[0]!=start[0] || curr[1]!=start[1]) break;
            run_end++;
        }
        int run_len = run_end - i;
        if (run_len < 2) {
            seen_boundary = 1;
        } else if (run_len > 2) {
            fprintf(stderr, "ch_quickhull3D: non-manifold edge (%d,%d) shared by %d faces.\n",
                    start[0], start[1], run_len);
            return 1;
        }
        i = run_end;
    }
    return seen_boundary ? 1 : 0;
}
#endif /* NDEBUG */


/* ===== Conflict lists ===== */

typedef struct { int pid; double dist; int next; } s_cl_node;

static int farthest_visible_face(const s_points *points, const int *faces,
                                  int from, int to, s_point p, double EPS_degenerate,
                                  double *out_dist)
{
    int best_f = -1;
    double best_d = 0.0;
    for (int f = from; f < to; f++) {
        s_point fv[3]; vertices_face(points, faces, f, fv);
        double d = signed_distance_point_to_plane(p, fv, EPS_degenerate);
        if (d > best_d) { best_d = d; best_f = f; }
    }
    if (out_dist) *out_dist = best_d;
    return best_f;
}

static int cl_alloc(s_dynarray *cl_arena, int *cl_freelist)
{
    if (*cl_freelist >= 0) {
        int node = *cl_freelist;
        *cl_freelist = ((s_cl_node*)cl_arena->items)[node].next;
        return node;
    }
    s_cl_node z = {0, 0.0, -1};
    if (!dynarray_push(cl_arena, &z)) return -1;
    return (int)cl_arena->N - 1;
}

static void cl_insert(s_cl_node *arena, int *cl_head, int node, int face)
{
    double dist = arena[node].dist;
    int *prev = &cl_head[face];
    while (*prev >= 0 && arena[*prev].dist >= dist) prev = &arena[*prev].next;
    arena[node].next = *prev;
    *prev = node;
}

static void build_conflict_lists(const s_points *points,
                                   const bool *buff_isused, const int *faces, int Nfaces,
                                   s_dynarray *cl_arena, int *cl_head, int *cl_freelist,
                                   double EPS_degenerate)
{
    for (int i = 0; i < points->N; i++) {
        if (buff_isused[i]) continue;
        double best_d;
        int best_f = farthest_visible_face(points, faces, 0, Nfaces,
                                            points->p[i], EPS_degenerate, &best_d);
        if (best_f < 0) continue;
        int node = cl_alloc(cl_arena, cl_freelist);
        if (node < 0) return;
        ((s_cl_node*)cl_arena->items)[node].pid  = i;
        ((s_cl_node*)cl_arena->items)[node].dist = best_d;
        cl_insert(cl_arena->items, cl_head, node, best_f);
    }
}

static int find_eye_point(s_cl_node *arena, int *cl_head, int Nfaces,
                           const bool *buff_isused,
                           int *cl_freelist, int *out_f_eye)
{
    int    best_pid  = -1;
    double best_dist = -1.0;
    *out_f_eye = -1;
    for (int f = 0; f < Nfaces; f++) {
        /* Drain stale heads: point already added to hull */
        while (cl_head[f] >= 0 && buff_isused[arena[cl_head[f]].pid]) {
            int old = cl_head[f];
            cl_head[f] = arena[old].next;
            arena[old].next = *cl_freelist; *cl_freelist = old;
        }
        if (cl_head[f] < 0) continue;
        if (arena[cl_head[f]].dist > best_dist) {
            best_dist = arena[cl_head[f]].dist;
            best_pid  = arena[cl_head[f]].pid;
            *out_f_eye = f;
        }
    }
    return best_pid;
}

static int collect_pending_nodes(s_dynarray *cl_arena, int Nfaces, const bool *faces_isvisible,
                                   const bool *buff_isused,
                                   int *cl_head, int *cl_freelist)
{
    s_cl_node *ar = cl_arena->items;
    int cl_pending = -1;
    for (int f = 0; f < Nfaces; f++) {
        if (!faces_isvisible[f]) continue;
        int node = cl_head[f]; cl_head[f] = -1;
        while (node >= 0) {
            int nx = ar[node].next;
            if (buff_isused[ar[node].pid]) {
                ar[node].next = *cl_freelist; *cl_freelist = node;
            } else {
                ar[node].next = cl_pending; cl_pending = node;
            }
            node = nx;
        }
    }
    return cl_pending;
}

static void repartition_points(int pending, const s_points *points,
                                const int *faces, int first_new_face, int Nfaces,
                                s_dynarray *cl_arena, int *cl_head, int *cl_freelist,
                                const bool *buff_isused, double EPS_degenerate)
{
    while (pending >= 0) {
        s_cl_node *ar = cl_arena->items;
        int next = ar[pending].next;
        int pid  = ar[pending].pid;

        if (!buff_isused[pid]) {
            double best_d;
            int best_f = farthest_visible_face(points, faces, first_new_face, Nfaces,
                                                points->p[pid], EPS_degenerate, &best_d);
            if (best_f < 0)
                best_f = farthest_visible_face(points, faces, 0, first_new_face,
                                               points->p[pid], EPS_degenerate, &best_d);
            if (best_f >= 0) {
                ar = cl_arena->items;
                ar[pending].pid = pid; ar[pending].dist = best_d;
                cl_insert(ar, cl_head, pending, best_f);
                pending = next;
                continue;
            }
        }
        /* Interior or already on hull: free the node */
        ar[pending].next = *cl_freelist; *cl_freelist = pending;
        pending = next;
    }
}



/* ===== Main algorithm ===== */

int quickhull_3d(const s_points *in_vertices, double EPS_degenerate,
                 bool buff_isused[in_vertices->N],
                 int **out_faces, int *N_out_faces)
{
    if (!in_vertices || !buff_isused || !out_faces || !N_out_faces) return -1;
    *out_faces = NULL; *N_out_faces = 0;

    s_dynarray faces = {0}, adj = {0}, faces_isvisible = {0},
               horizon = {0}, AUX_horizon_nbr = {0},
               AUX_cop_fids = {0},
               dfs_vis = {0}, dfs_stk = {0}, cl_arena = {0};
    int  *buff_pmap          = NULL;
    int  *cl_head = NULL;
    int   cl_freelist = -1;
    int   Nfaces = 0;
    int   ret;

    /* Helper Macro */
    #define ALLOC_FAIL(msg) do { fprintf(stderr, "ch_quickhull3D: " msg "\n"); goto error; } while(0)

    if (in_vertices->N <= 3) goto error_degenerate;

    Nfaces = 4;
    faces = dynarray_initialize(sizeof(int), 12);
    adj   = dynarray_initialize(sizeof(int), 12);
    if (!faces.items || !adj.items) ALLOC_FAIL("malloc failed.");
    adj.N = 12;
    if (!initial_tetrahedron(in_vertices, buff_isused, EPS_degenerate,
                              faces.items, adj.items))
        goto error_degenerate;
    if (in_vertices->N == 4) {
        *out_faces = faces.items; *N_out_faces = 4;
        dynarray_free(&adj);
        return 1;
    }

    buff_pmap          = malloc(CH_MAX_NUM_FACES * sizeof(int));   if (!buff_pmap)          ALLOC_FAIL("malloc failed.");
    dfs_vis = dynarray_initialize(sizeof(bool), 64);        if (!dfs_vis.items)  ALLOC_FAIL("malloc failed.");
    dfs_stk = dynarray_initialize(sizeof(int),  64);        if (!dfs_stk.items)  ALLOC_FAIL("malloc failed.");
    cl_arena = dynarray_initialize(sizeof(s_cl_node), 64);  if (!cl_arena.items) ALLOC_FAIL("malloc failed.");
    cl_head = malloc(CH_MAX_NUM_FACES * sizeof(int));       if (!cl_head) ALLOC_FAIL("malloc failed.");
    for (int i = 0; i < CH_MAX_NUM_FACES; i++) cl_head[i] = -1;

    build_conflict_lists(in_vertices, buff_isused,
                         faces.items, Nfaces, &cl_arena, cl_head, &cl_freelist,
                         EPS_degenerate);

    faces_isvisible = dynarray_initialize(sizeof(bool), 0);
    AUX_horizon_nbr = dynarray_initialize(sizeof(int), 0);
    horizon         = dynarray_initialize(sizeof(int), 0);
    AUX_cop_fids    = dynarray_initialize(sizeof(int), 0);
    if (!faces_isvisible.items || !AUX_horizon_nbr.items || !horizon.items ||
        !AUX_cop_fids.items)
        ALLOC_FAIL("malloc failed.");

    for (;;) {
        /* 1. Select the eye point: farthest from any conflict-list face */
        int f_eye = -1;
        int current_id = find_eye_point(cl_arena.items, cl_head, Nfaces,
                                        buff_isused, &cl_freelist, &f_eye);
        if (current_id < 0) break;

        /* Pop the eye-point node from its conflict list */
        { s_cl_node *ar = cl_arena.items;
          int popped = cl_head[f_eye];
          cl_head[f_eye] = ar[popped].next;
          ar[popped].next = cl_freelist; cl_freelist = popped; }

        /* 2. Classify faces visible from current_id */
        int N_vf = visible_faces_from_point(
            in_vertices, Nfaces, faces.items, adj.items,
            in_vertices->p[current_id], EPS_degenerate,
            &AUX_cop_fids, &faces_isvisible);
        if (N_vf == -1) { fprintf(stderr, "ch_quickhull3D: visible_faces error.\n"); goto error; }
        if (N_vf == Nfaces) { fprintf(stderr, "ch_quickhull3D: point sees all faces.\n"); goto error; }
        if (N_vf == 0) { continue; }
        buff_isused[current_id] = true;

        /* 3. Collect pending nodes from visible faces' conflict lists */
        int cl_pending = collect_pending_nodes(&cl_arena, Nfaces, faces_isvisible.items,
                                               buff_isused, cl_head, &cl_freelist);

        /* 4. Extract horizon */
        int Nhorizon;
        {
            int f_start = -1;
            bool *vis_ = faces_isvisible.items;
            for (int fi = 0; fi < Nfaces && f_start < 0; fi++) if (vis_[fi]) f_start = fi;
            if (!dynarray_ensure_capacity(&dfs_vis, Nfaces) ||
                !dynarray_ensure_capacity(&dfs_stk, Nfaces))
                { fprintf(stderr, "ch_quickhull3D: DFS scratch alloc error.\n"); goto error; }
            memset(dfs_vis.items, 0, Nfaces * sizeof(bool));
            if (!extract_horizon_dfs(faces.items, adj.items, faces_isvisible.items,
                                     dfs_vis.items, dfs_stk.items, f_start,
                                     &Nhorizon, &N_vf, &horizon, &AUX_horizon_nbr))
                { fprintf(stderr, "ch_quickhull3D: horizon extraction error.\n"); goto error; }
        }
        sort_horizon_loop(horizon.items, AUX_horizon_nbr.items, Nhorizon);


        /* 5. Delete visible faces and compact cl_head / horizon_nbr */
        int old_Nfaces = Nfaces;
        delete_visible_faces(Nfaces, faces.items, adj.items, faces_isvisible.items,
                             buff_isused, buff_pmap);
        for (int j = 0; j < old_Nfaces; j++) {
            if (buff_pmap[j] < 0) continue;
            cl_head[buff_pmap[j]] = cl_head[j];
        }
        for (int j = 0; j < Nhorizon; j++)
            ((int*)AUX_horizon_nbr.items)[j] = buff_pmap[((int*)AUX_horizon_nbr.items)[j]];
        Nfaces -= N_vf;
        faces.N = Nfaces * 3;
        /* Restore isused for vertices present in surviving faces */
        { int *fp = faces.items;
          for (int fi=0; fi<Nfaces; fi++) {
              buff_isused[fp[fi*3+0]] = buff_isused[fp[fi*3+1]] = buff_isused[fp[fi*3+2]] = true;
          }
        }

        /* 6. Add new faces connecting horizon to current_id */
        int first_new_face = Nfaces;
        if (!add_faces_from_horizon(in_vertices, buff_isused, Nfaces, &faces, &adj,
                                    Nhorizon, horizon.items, AUX_horizon_nbr.items, current_id))
            { fprintf(stderr, "ch_quickhull3D: add_faces error.\n"); goto error; }
        Nfaces += Nhorizon;
        faces.N = adj.N = Nfaces * 3;

        /* 7. Repartition pending points into new faces' conflict lists */
        for (int f = first_new_face; f < Nfaces; f++) cl_head[f] = -1;
        repartition_points(cl_pending, in_vertices, faces.items,
                           first_new_face, Nfaces, &cl_arena, cl_head, &cl_freelist,
                           buff_isused, EPS_degenerate);
    }

#ifndef NDEBUG
    /* Topology check */
    { s_dynarray edges = dynarray_initialize(sizeof(int), 0);
      if (!edges.items) goto error;
      if (topology_is_valid(Nfaces, faces.items, &edges)) {
          write_points_to_csv("problematic.csv", "w", in_vertices);
          fprintf(stderr, "ch_quickhull3D: topology invalid  -- problematic points written to problematic.csv"
                          " (EPS=%g)\n", EPS_degenerate);
          dynarray_free(&edges); goto error;
      }
      dynarray_free(&edges);
    }
#endif

    #undef ALLOC_FAIL

    free(buff_pmap); free(cl_head);
    dynarray_free(&cl_arena); dynarray_free(&adj);
    dynarray_free(&faces_isvisible); dynarray_free(&horizon); dynarray_free(&AUX_horizon_nbr);
    dynarray_free(&AUX_cop_fids);
    dynarray_free(&dfs_vis); dynarray_free(&dfs_stk);
    *out_faces = realloc(faces.items, Nfaces * 3 * sizeof(int));
    *N_out_faces = Nfaces;
    return 1;

error_degenerate:
    ret = 0;
    goto cleanup;
error:
    fprintf(stderr, "ch_quickhull3D: error. Faces N = %d / %d\n", Nfaces, CH_MAX_NUM_FACES);
    ret = -1;
cleanup:
    if (faces.items)         dynarray_free(&faces);
    if (adj.items)           dynarray_free(&adj);
    if (buff_pmap)           free(buff_pmap);
    if (cl_head)             free(cl_head);
    if (cl_arena.items)      dynarray_free(&cl_arena);
    if (dfs_vis.items)       dynarray_free(&dfs_vis);
    if (dfs_stk.items)       dynarray_free(&dfs_stk);
    if (faces_isvisible.items)   dynarray_free(&faces_isvisible);
    if (horizon.items)           dynarray_free(&horizon);
    if (AUX_horizon_nbr.items)   dynarray_free(&AUX_horizon_nbr);
    if (AUX_cop_fids.items)      dynarray_free(&AUX_cop_fids);
    *out_faces = NULL; *N_out_faces = 0;
    return ret;
}

