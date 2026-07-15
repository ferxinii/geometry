#ifndef GEOMETRY_TRIMESH_H
#define GEOMETRY_TRIMESH_H

#include "points.h"

/*
 * Closed, oriented triangle surface mesh.
 * Structurally identical to s_convh but with no convexity constraint,
 * plus face-face adjacency.
 *
 * Convention: local edge j of face i is the edge OPPOSITE vertex j,
 * i.e. it connects faces[i*3 + (j+1)%3] and faces[i*3 + (j+2)%3].
 */
typedef struct trimesh {
    s_points  points;      /* mesh vertices */
    int       Nf;          /* number of triangular faces */
    int      *faces;       /* flat: [v0,v1,v2, v0',v1',v2', ...], length Nf*3 */
    s_point  *fnormals;    /* one unnormalized outward normal per face, length Nf */
    int      *adjacency;   /* adjacency[i*3 + j] = index of face adjacent to face i
                              across local edge j; -1 if none */
} s_trimesh;

#define trimesh_NAN (s_trimesh){0}

int       trimesh_is_valid(const s_trimesh *m);
void      free_trimesh(s_trimesh *m);
s_trimesh copy_trimesh(const s_trimesh *m);
double    volume_trimesh(const s_trimesh *m);

/*
 * Build from raw arrays. Computes normals and adjacency.
 * Validates: closed manifold (every edge has exactly 2 incident faces).
 * Winding is made consistent via BFS; global outward orientation via
 * the divergence-theorem signed volume.
 * Returns trimesh_NAN on error.
 */
s_trimesh trimesh_from_arrays(const s_point *points, int N_points,
                              const int *faces, int N_faces,
                              double EPS_DEG);

/*
 * Load a closed triangle mesh from a Wavefront OBJ file.
 * Handles v/vt, v/vt/vn, and v//vn face formats; polygonal faces are
 * fan-triangulated.  Normals/adjacency are computed; winding and outward
 * orientation are fixed automatically (same as trimesh_from_arrays).
 * Returns trimesh_NAN on error.
 */
s_trimesh trimesh_from_obj(const char *fname, double EPS_DEG);

/*
 * Returns 1 if p is strictly inside the closed trimesh, 0 if outside,
 * -1 on error. Uses ray casting in the +x direction. On degenerate
 * intersection retries with a perturbed direction up to max_retries times.
 */
int point_in_trimesh(const s_trimesh *m, s_point p, double EPS_DEG, int max_retries);

/*
 * Generalized winding number of p wrt the mesh: the total signed solid angle
 * subtended by the faces, / 4pi (van Oosterom-Strackee per triangle).  For a
 * CLOSED, consistently OUTWARD-oriented mesh (the constructors' invariant)
 * this is ~1 inside and ~0 outside, fractional only for points within
 * numerical noise of the surface (or inside cavities thinner than that
 * noise).  Pure per-query summation, O(Nf): no rays, no retries, no special
 * positions, no construction step that can fail -- robust for arbitrarily
 * degenerate valid meshes.
 */
double trimesh_winding_number(const s_trimesh *m, s_point p);

/* Inside test via winding number (> 0.5).  Prefer this over the ray-casting
 * point_in_trimesh when the mesh may contain degenerate geometry. */
int point_in_trimesh_winding(const s_trimesh *m, s_point p);

/* ---- degenerate-face repair --------------------------------------------
 *
 * Collapse-based repair of near-degenerate faces on a CLOSED MANIFOLD mesh.
 * A face is degenerate when its thinness q = 2*area / Lmax^2 falls below
 * q_min; the repair collapses the face's SHORTEST edge to its midpoint,
 * guarded by
 *   - the link condition (link(a) ^ link(b) == the two vertices opposite the
 *     edge), which preserves the closed-manifold invariant,
 *   - a normal flip guard (no incident face may reverse orientation),
 *   - a hard displacement bound: edges longer than
 *     max_collapse_len * bbox_diagonal are never collapsed.
 * Rounds repeat until a fixed point (or max_rounds, with a stderr warning).
 * Anything unfixable under the guards is left in place and counted in
 * stats->n_unfixable.
 *
 * NOTE: a vertex-onto-edge "weld" operation (for cap faces whose apex hovers
 * over a non-incident edge) was considered and deliberately left out:
 * calibration over the motivating datasets (Voronoi-cut lung cells) found no
 * cap faces in any input -- all degeneracies are blades with a microscopic
 * base edge, which the collapse removes.  Add the weld only if a real input
 * ever needs it.
 *
 * Geometric side effects are bounded and reported: every displaced vertex
 * moves by at most max_collapse_len * bbox_diag / 2 (stats->max_displacement
 * has the realized maximum) and stats->volume_drift the realized |dV|.
 *
 * Intended use: consumers that need a well-conditioned APPROXIMATION of the
 * surface (e.g. medial-axis computation, tetrahedralization oracles).  Do NOT
 * feed repaired surfaces to consumers that require the exact input geometry
 * (volume accounting, watertight cell partitions): independently repaired
 * neighbouring cells are NOT cross-cell consistent.
 */
typedef struct {
    double q_min;            /* face thinness threshold (0 disables repair)  */
    double max_collapse_len; /* max collapsible edge length, relative to the
                                bbox diagonal                                */
    int    max_rounds;
} s_trimesh_repair_opts;

#define TRIMESH_REPAIR_DEFAULTS     (s_trimesh_repair_opts){ .q_min = 1e-3, .max_collapse_len = 1e-2,                              .max_rounds = 50 }

typedef struct {
    int    n_collapses;
    int    n_rounds;
    int    n_unfixable;       /* degenerate faces left (guards refused)      */
    double max_displacement;  /* absolute                                    */
    double volume_drift;      /* |V_out - V_in|                              */
} s_trimesh_repair_stats;

/* Returns a NEW mesh (input untouched) satisfying the same closed-manifold
 * postcondition as the constructors.  opts == NULL uses the defaults;
 * stats == NULL is allowed.  Returns trimesh_NAN on error (invalid input or
 * allocation failure). */
s_trimesh trimesh_repaired(const s_trimesh *in,
                           const s_trimesh_repair_opts *opts,
                           s_trimesh_repair_stats *stats);

#endif
