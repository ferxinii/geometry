#include "trimesh.h"
#include "gtests.h"
#include "hash.h"
#include "dynarray.h"
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>


/* -----------------------------------------------------------------------
 * Internal hash helpers (edge keys are sorted int[2] pairs)
 * ----------------------------------------------------------------------- */

static size_t mix64(uint64_t x)
{
    x ^= x >> 33;
    x *= 0xff51afd7ed558ccdULL;
    x ^= x >> 33;
    x *= 0xc4ceb9fe1a85ec53ULL;
    x ^= x >> 33;
    return (size_t)x;
}

static size_t hash_fun_edge(const void *key)
{
    const int *v = key;
    uint64_t h = ((uint64_t)(unsigned)v[0] * 73856093u) ^
                 ((uint64_t)(unsigned)v[1] * 19349663u);
    return mix64(h);
}

static bool hash_eq_edge(const void *k1, const void *k2)
{
    const int *a = k1, *b = k2;
    return a[0] == b[0] && a[1] == b[1];
}


/* -----------------------------------------------------------------------
 * trimesh_is_valid / free / copy
 * ----------------------------------------------------------------------- */

int trimesh_is_valid(const s_trimesh *m)
{
    return m->points.p != NULL && m->points.N > 0 &&
           m->Nf > 0 &&
           m->faces     != NULL &&
           m->fnormals  != NULL &&
           m->adjacency != NULL;
}

void free_trimesh(s_trimesh *m)
{
    free_points(&m->points);
    free(m->faces);
    free(m->fnormals);
    free(m->adjacency);
    memset(m, 0, sizeof(s_trimesh));
}

s_trimesh copy_trimesh(const s_trimesh *m)
{
    if (!trimesh_is_valid(m)) return trimesh_NAN;

    s_trimesh out = trimesh_NAN;
    out.points = copy_points(&m->points);
    if (!out.points.p) return trimesh_NAN;

    out.Nf = m->Nf;

    out.faces = malloc(sizeof(int) * m->Nf * 3);
    if (!out.faces) goto fail;
    memcpy(out.faces, m->faces, sizeof(int) * m->Nf * 3);

    out.fnormals = malloc(sizeof(s_point) * m->Nf);
    if (!out.fnormals) goto fail;
    memcpy(out.fnormals, m->fnormals, sizeof(s_point) * m->Nf);

    out.adjacency = malloc(sizeof(int) * m->Nf * 3);
    if (!out.adjacency) goto fail;
    memcpy(out.adjacency, m->adjacency, sizeof(int) * m->Nf * 3);

    return out;

fail:
    free_trimesh(&out);
    return trimesh_NAN;
}


/* -----------------------------------------------------------------------
 * trimesh_volume
 * ----------------------------------------------------------------------- */

double volume_trimesh(const s_trimesh *m)
{
    s_point c = point_average(&m->points);
    double V = 0.0;
    for (int i = 0; i < m->Nf; i++) {
        s_point v0 = subtract_points(m->points.p[m->faces[i*3+0]], c);
        s_point v1 = subtract_points(m->points.p[m->faces[i*3+1]], c);
        s_point v2 = subtract_points(m->points.p[m->faces[i*3+2]], c);
        V += dot_prod(v0, cross_prod(v1, v2));
    }
    return V / 6.0;
}


/* -----------------------------------------------------------------------
 * trimesh_from_arrays
 * ----------------------------------------------------------------------- */

typedef struct {
    int face;
    int local_edge;
    int count;   /* 0 = empty, 1 = one face seen, 2 = both faces seen */
} s_edge_entry;

s_trimesh trimesh_from_arrays(const s_point *points, int N_points,
                              const int *faces, int N_faces,
                              double EPS_DEG)
{
    (void)EPS_DEG;

    if (!points || N_points <= 0 || !faces || N_faces <= 0)
        return trimesh_NAN;

    s_trimesh m = trimesh_NAN;

    /* 1. Copy points */
    m.points.N = N_points;
    m.points.p = malloc(sizeof(s_point) * N_points);
    if (!m.points.p) return trimesh_NAN;
    memcpy(m.points.p, points, sizeof(s_point) * N_points);

    /* 2. Copy faces */
    m.Nf = N_faces;
    m.faces = malloc(sizeof(int) * N_faces * 3);
    if (!m.faces) goto fail;
    memcpy(m.faces, faces, sizeof(int) * N_faces * 3);

    /* 3. Per-face normals: n = cross(B-A, C-A) */
    m.fnormals = malloc(sizeof(s_point) * N_faces);
    if (!m.fnormals) goto fail;
    for (int i = 0; i < N_faces; i++) {
        s_point A = m.points.p[m.faces[i*3+0]];
        s_point B = m.points.p[m.faces[i*3+1]];
        s_point C = m.points.p[m.faces[i*3+2]];
        m.fnormals[i] = cross_prod(subtract_points(B, A), subtract_points(C, A));
    }

    /* 4. Build adjacency via edge hash */
    m.adjacency = malloc(sizeof(int) * N_faces * 3);
    if (!m.adjacency) goto fail;
    for (int i = 0; i < N_faces * 3; i++) m.adjacency[i] = -1;

    s_hash_table ht;
    if (!hash_init(&ht, sizeof(int)*2, sizeof(s_edge_entry),
                   (size_t)N_faces * 4, (size_t)N_faces * 3,
                   hash_fun_edge, hash_eq_edge, NULL))
        goto fail;

    bool adj_ok = true;
    for (int i = 0; i < N_faces && adj_ok; i++) {
        for (int j = 0; j < 3 && adj_ok; j++) {
            /* local edge j connects (j+1)%3 and (j+2)%3 */
            int va = m.faces[i*3 + (j+1)%3];
            int vb = m.faces[i*3 + (j+2)%3];
            int key[2] = { va < vb ? va : vb, va < vb ? vb : va };

            s_edge_entry *e = hash_get_or_create(&ht, key);
            if (!e) { adj_ok = false; break; }

            if (e->count == 0) {
                e->face = i; e->local_edge = j; e->count = 1;
            } else if (e->count == 1) {
                m.adjacency[i*3 + j]                        = e->face;
                m.adjacency[e->face*3 + e->local_edge]      = i;
                e->count = 2;
            } else {
                fprintf(stderr, "trimesh_from_arrays: non-manifold edge (%d,%d).\n", va, vb);
                adj_ok = false;
            }
        }
    }
    hash_free(&ht);
    if (!adj_ok) goto fail;

    /* Validate: every edge must have exactly 2 incident faces (closed mesh) */
    for (int i = 0; i < N_faces * 3; i++) {
        if (m.adjacency[i] == -1) {
            fprintf(stderr, "trimesh_from_arrays: boundary edge on face %d local %d "
                    "(mesh not closed).\n", i / 3, i % 3);
            goto fail;
        }
    }

    /* 5. BFS to propagate consistent winding across the mesh.
     *
     * For a consistently oriented mesh adjacent faces always traverse their
     * shared edge in OPPOSITE directions.  If we find a neighbour that
     * traverses the shared edge in the SAME direction as the current face,
     * we flip it (swap vertex positions 1<->2, negate normal, swap adjacency
     * entries 1<->2 to keep the adjacency table correct after the flip).
     */
    {
        bool *visited = calloc((size_t)N_faces, sizeof(bool));
        int  *queue   = malloc(sizeof(int) * (size_t)N_faces);
        if (!visited || !queue) { free(visited); free(queue); goto fail; }

        for (int start = 0; start < N_faces; start++) {
            if (visited[start]) continue;

            int qhead = 0, qtail = 0;
            visited[start] = true;
            queue[qtail++] = start;

            while (qhead < qtail) {
                int fi = queue[qhead++];
                for (int j = 0; j < 3; j++) {
                    int nbr = m.adjacency[fi*3 + j];
                    if (visited[nbr]) continue;
                    visited[nbr] = true;

                    /* Directed edge of fi at local j: va -> vb */
                    int va = m.faces[fi*3 + (j+1)%3];
                    int vb = m.faces[fi*3 + (j+2)%3];

                    /* Find which local edge of nbr is this shared edge */
                    int k = -1;
                    for (int kk = 0; kk < 3; kk++) {
                        int ea = m.faces[nbr*3 + (kk+1)%3];
                        int eb = m.faces[nbr*3 + (kk+2)%3];
                        if ((ea == va && eb == vb) || (ea == vb && eb == va))
                            { k = kk; break; }
                    }
                    /* k == -1 would mean corrupted adjacency; can't happen here */

                    int ea = m.faces[nbr*3 + (k+1)%3];
                    int eb = m.faces[nbr*3 + (k+2)%3];

                    if (ea == va && eb == vb) {
                        /* Same direction -> flip neighbour */
                        int tmp;
                        tmp = m.faces[nbr*3+1];
                        m.faces[nbr*3+1] = m.faces[nbr*3+2];
                        m.faces[nbr*3+2] = tmp;

                        m.fnormals[nbr] = scale_point(m.fnormals[nbr], -1.0);

                        /* Swapping vertex positions 1<->2 also swaps local edges 1<->2
                         * (local edge k is opposite vertex at position k), so the
                         * adjacency entries for positions 1 and 2 must be swapped. */
                        tmp = m.adjacency[nbr*3+1];
                        m.adjacency[nbr*3+1] = m.adjacency[nbr*3+2];
                        m.adjacency[nbr*3+2] = tmp;
                    }

                    queue[qtail++] = nbr;
                }
            }
        }
        free(visited);
        free(queue);
    }

    /* 6. Global orientation via divergence-theorem signed volume.
     *
     * V = (1/6) sum_i dot(v0_i, cross(v1_i, v2_i))
     *
     * Vertices are shifted by the mesh centroid before summing to reduce
     * floating-point cancellation for meshes far from the origin.
     * V < 0 means normals point inward; flip the entire mesh.
     */
    {
        s_point c = point_average(&m.points);
        double V = 0.0;
        for (int i = 0; i < N_faces; i++) {
            s_point v0 = subtract_points(m.points.p[m.faces[i*3+0]], c);
            s_point v1 = subtract_points(m.points.p[m.faces[i*3+1]], c);
            s_point v2 = subtract_points(m.points.p[m.faces[i*3+2]], c);
            V += dot_prod(v0, cross_prod(v1, v2));
        }
        V /= 6.0;

        if (V < 0.0) {
            for (int i = 0; i < N_faces; i++) {
                int tmp = m.faces[i*3+1];
                m.faces[i*3+1] = m.faces[i*3+2];
                m.faces[i*3+2] = tmp;

                m.fnormals[i] = scale_point(m.fnormals[i], -1.0);

                tmp = m.adjacency[i*3+1];
                m.adjacency[i*3+1] = m.adjacency[i*3+2];
                m.adjacency[i*3+2] = tmp;
            }
        }
    }

    return m;

fail:
    free_trimesh(&m);
    return trimesh_NAN;
}


/* -----------------------------------------------------------------------
 * point_in_trimesh  (Phase 2 - ray casting)
 * ----------------------------------------------------------------------- */

int point_in_trimesh(const s_trimesh *m, s_point p, double EPS_DEG, int max_retries)
{
    if (!trimesh_is_valid(m)) return -1;

    /* TODO (post-MVP): replace arbitrary perturbation coefficients with a
     * pre-computed set of well-spread directions (e.g. icosahedron vertices)
     * so retries have a geometric coverage guarantee. */

    /* TODO (post-MVP): ray_len should be 2x the diagonal of the mesh AABB,
     * computed from m->points, and exposed as a parameter so callers with
     * large meshes are not silently under-covered by the fixed 1e9. */

    double dx = 1.0, dy = 0.0, dz = 0.0;
    double ray_len = 1e9;

    for (int attempt = 0; attempt <= max_retries; attempt++) {
        s_point ray_end;
        ray_end.x = p.x + ray_len * dx;
        ray_end.y = p.y + ray_len * dy;
        ray_end.z = p.z + ray_len * dz;
        s_point segment[2] = { p, ray_end };

        int count = 0;
        bool retry = false;

        for (int i = 0; i < m->Nf; i++) {
            s_point tri[3] = {
                m->points.p[m->faces[i*3+0]],
                m->points.p[m->faces[i*3+1]],
                m->points.p[m->faces[i*3+2]]
            };
            e_intersect_type r = test_segment_triangle_intersect_3D(
                    segment, tri, EPS_DEG, 0.0);
            if      (r == INTERSECT_NONDEGENERATE) count++;
            else if (r == INTERSECT_DEGENERATE || r == INTERSECT_ERROR)
                { retry = true; break; }
        }

        if (!retry) return count % 2;   /* odd = inside */

        /* Perturb direction slightly and retry */
        dx += 0.1 * (attempt + 1);
        dy += 0.07 * (attempt + 1);
        dz += 0.13 * (attempt + 1);
        double len = sqrt(dx*dx + dy*dy + dz*dz);
        dx /= len; dy /= len; dz /= len;
    }

    fprintf(stderr, "point_in_trimesh: degenerate intersection persists after %d retries.\n",
            max_retries);
    return -1;
}


/* -----------------------------------------------------------------------
 * trimesh_from_obj
 * ----------------------------------------------------------------------- */

s_trimesh trimesh_from_obj(const char *fname, double EPS_DEG)
{
    FILE *f = fopen(fname, "r");
    if (!f) {
        fprintf(stderr, "trimesh_from_obj: cannot open '%s'\n", fname);
        return trimesh_NAN;
    }

    s_dynarray dverts = dynarray_initialize(sizeof(s_point),    256);
    s_dynarray dfaces = dynarray_initialize(3 * sizeof(int),    512);
    if (!dverts.items || !dfaces.items) goto fail;

    char line[1024];
    int lineno = 0;
    while (fgets(line, sizeof(line), f)) {
        lineno++;

        if (line[0] == 'v' && (line[1] == ' ' || line[1] == '\t')) {
            /* Vertex position */
            double x, y, z;
            if (sscanf(line + 1, "%lf %lf %lf", &x, &y, &z) != 3) {
                fprintf(stderr, "trimesh_from_obj: bad vertex on line %d\n", lineno);
                goto fail;
            }
            s_point v = {.x=x, .y=y, .z=z};
            if (!dynarray_push(&dverts, &v)) goto fail;

        } else if (line[0] == 'f' && (line[1] == ' ' || line[1] == '\t')) {
            /* Face: parse vertex indices, handling v, v/vt, v/vt/vn, v//vn */
            const char *p = line + 1;
            int idx[4], n = 0;
            while (n < 4) {
                while (*p == ' ' || *p == '\t') p++;
                if (*p == '\n' || *p == '\r' || *p == '\0') break;
                int vi;
                if (sscanf(p, "%d", &vi) != 1) break;
                idx[n++] = vi > 0 ? vi - 1 : (int)dverts.N + vi;   /* 1-indexed, negative OK */
                while (*p && *p != ' ' && *p != '\t' && *p != '\n' && *p != '\r') p++;
            }
            if (n < 3) {
                fprintf(stderr, "trimesh_from_obj: degenerate face on line %d\n", lineno);
                goto fail;
            }
            /* Fan triangulation: (0,1,2), (0,2,3), ... */
            for (int t = 1; t < n - 1; t++) {
                int tri[3] = {idx[0], idx[t], idx[t+1]};
                if (!dynarray_push(&dfaces, tri)) goto fail;
            }
        }
        /* All other lines (vn, vt, #, g, s, mtllib, usemtl, ...) are ignored */
    }
    fclose(f);
    f = NULL;

    if (dverts.N == 0 || dfaces.N == 0) {
        fprintf(stderr, "trimesh_from_obj: no geometry found in '%s'\n", fname);
        goto fail;
    }

    s_trimesh m = trimesh_from_arrays(dverts.items, (int)dverts.N,
                                      dfaces.items,  (int)dfaces.N, EPS_DEG);
    dynarray_free(&dverts);
    dynarray_free(&dfaces);
    return m;

fail:
    if (f) fclose(f);
    dynarray_free(&dverts);
    dynarray_free(&dfaces);
    return trimesh_NAN;
}


/* -----------------------------------------------------------------------
 * generalized winding number
 * ----------------------------------------------------------------------- */

double trimesh_winding_number(const s_trimesh *m, s_point p)
{
    double w = 0.0;
    for (int i = 0; i < m->Nf; i++) {
        const int *f = &m->faces[i*3];
        const double ax = m->points.p[f[0]].x - p.x;
        const double ay = m->points.p[f[0]].y - p.y;
        const double az = m->points.p[f[0]].z - p.z;
        const double bx = m->points.p[f[1]].x - p.x;
        const double by = m->points.p[f[1]].y - p.y;
        const double bz = m->points.p[f[1]].z - p.z;
        const double cx = m->points.p[f[2]].x - p.x;
        const double cy = m->points.p[f[2]].y - p.y;
        const double cz = m->points.p[f[2]].z - p.z;

        const double la = sqrt(ax*ax + ay*ay + az*az);
        const double lb = sqrt(bx*bx + by*by + bz*bz);
        const double lc = sqrt(cx*cx + cy*cy + cz*cz);

        const double det = ax*(by*cz - bz*cy)
                         - ay*(bx*cz - bz*cx)
                         + az*(bx*cy - by*cx);
        const double den = la*lb*lc
                         + (ax*bx + ay*by + az*bz) * lc
                         + (bx*cx + by*cy + bz*cz) * la
                         + (cx*ax + cy*ay + cz*az) * lb;

        w += atan2(det, den);          /* Omega/2 for this face */
    }
    return w / (2.0 * M_PI);           /* sum(Omega) / 4pi */
}

int point_in_trimesh_winding(const s_trimesh *m, s_point p)
{
    return fabs(trimesh_winding_number(m, p)) > 0.5;
}
