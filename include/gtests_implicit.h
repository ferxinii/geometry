#ifndef _GTESTS_IMPLICIT_H
#define _GTESTS_IMPLICIT_H

#include "points.h"

#ifdef GEOMETRY_COMPILE_TESTS_IMPLICIT
/* Indirect (implicit-point) geometric predicates, exact arithmetic over
 * points that may be explicit or an implicit linear combination (LNC) of
 * two prior points 
 *
 * Backed by external/implicit_predicates (Attene, "Indirect predicates for
 * geometric constructions"); see that directory's README for the underlying
 * paper. All state lives in a caller-owned s_implicit_pred_reg -- no globals,
 * so independent builds (e.g. one per thread) can each own their own.
 */

typedef struct implicit_predicate_registry s_implicit_pred_reg;

/* Must be called once per process before any predicate (sets FPU rounding
 * mode). Idempotent. */
void itest_init(void);

/* One registry per DT/CDT build; no state is shared between registries. */
s_implicit_pred_reg *itest_registry_create(void);
void itest_registry_destroy(s_implicit_pred_reg *reg);

/* Register an explicit point at registry index id. */
void itest_point_set_explicit(s_implicit_pred_reg *reg, int id, s_point p);

/* Register a Steiner point at registry index id.
 * Represents t*V[v1] + (1-t)*V[v2] (paper convention).
 * v1 and v2 must already be registered (explicit or LNC) in this registry. */
void itest_point_set_lnc(s_implicit_pred_reg *reg, int id, int v1, int v2, double t);

/* orient3D -- returns -1, 0, or +1.
 * Arguments are registry point indices; any may be explicit or LNC. */
int itest_orient3d(s_implicit_pred_reg *reg, int a, int b, int c, int d);

/* inSphere -- returns -1 (inside), 0 (on), +1 (outside).
 * Arguments are registry point indices; any may be explicit or LNC. */
int itest_insphere(s_implicit_pred_reg *reg, int a, int b, int c, int d, int e);

/* SoS-wrapped inSphere for gift-wrapping (Algorithm 1 from the paper).
 * i1..i4 are the tet vertices, i5 is the query point.
 * Returns -1 if i5 is inside or on circumsphere, +1 otherwise.
 * Never returns 0 (ties broken by global index order). */
int itest_perturbed_insphere(s_implicit_pred_reg *reg, int i1, int i2, int i3, int i4, int i5);

/* TRUE if point p is inside or on the boundary of tet <t0,t1,t2,t3>. */
int itest_point_in_tet(s_implicit_pred_reg *reg, int p, int t0, int t1, int t2, int t3);

/* TRUE if segment <s1,s2> intersects triangle <v0,v1,v2>
 * (not coplanar; interior of segment only). */
int itest_segment_crosses_triangle(s_implicit_pred_reg *reg, int s1, int s2, int v0, int v1, int v2);

/* TRUE if interior of s1-s2 crosses interior of <v0,v1,v2> at a single point.
 * Returns false if either endpoint is coplanar with the triangle. */
int itest_inner_seg_crosses_inner_tri(s_implicit_pred_reg *reg, int s1, int s2, int v0, int v1, int v2);

/* TRUE if interior of coplanar segments a-b and p-q cross at a single point.
 * Points are assumed coplanar; tries all three axis projections. */
int itest_inner_segs_cross(s_implicit_pred_reg *reg, int a, int b, int p, int q);

/* TRUE if point p is strictly inside (not on boundary of) triangle <a,b,c>.
 * Points are assumed coplanar; tries all three axis projections. */
int itest_point_in_inner_tri(s_implicit_pred_reg *reg, int p, int a, int b, int c);

/* Point p vs triangle <a,b,c>: -1 OUTSIDE, 0 on BOUNDARY, +1 strictly INSIDE.
 * Mirrors test_point_in_triangle_3D(.,0,0) (OUT if not coplanar). */
int itest_point_in_triangle(s_implicit_pred_reg *reg, int p, int a, int b, int c);

/* TRUE if the CLOSED segments <s1,s2> and <p,q> cross (coplanar inputs). */
int itest_segments_cross(s_implicit_pred_reg *reg, int s1, int s2, int p, int q);

/* Get approximate float coordinates of any registered point (explicit or
 * LNC). For DT insertion of Steiner points. Returns 0 if id is out of range. */
int itest_get_approx_coords(s_implicit_pred_reg *reg, int id, double *out);

#else  // GEOMETRY_COMPILE_TESTS_IMPLICIT
static_assert(0, "Need to turn ON CMake option 'GEOMETRY_COMPILE_TESTS_IMPLICIT' (or flag GEOMETRY_COMPILE_TESTS_IMPLICIT)");
#endif  // GEOMETRY_COMPILE_TESTS_IMPLICIT

#endif
