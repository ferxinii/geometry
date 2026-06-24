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

#endif
