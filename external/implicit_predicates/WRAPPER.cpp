#include "indirect_predicates.h"
#include "dynarray.h"
#include "points.h"

/* point_entry maps a registry-visible point id to where its backing object
 * actually lives.
 *   type/exact_local_id : the EXACT representation used for direct predicate
 *     calls on this id -- explicit_pts[exact_local_id] if type==EXPLICIT3D,
 *     lnc_pts[exact_local_id] if type==LNC.
 *   shadow_local_id : explicit_pts[shadow_local_id] is ALWAYS a valid,
 *     rounded-double explicitPoint3D for this id, regardless of type. The
 *     underlying library's LNC can only combine two EXPLICIT points, so a
 *     nested Steiner (one whose endpoint is itself a prior Steiner) is built
 *     from its parent's shadow, not the parent's exact LNC. For an
 *     EXPLICIT3D id, shadow_local_id == exact_local_id (same slot). */
typedef struct {
    Point_Type type;
    int exact_local_id;
    int shadow_local_id;
} point_entry;

struct implicit_predicate_registry {
    s_dynarray explicit_pts; /* s_dynarray of explicitPoint3D */
    s_dynarray lnc_pts;      /* s_dynarray of implicitPoint3D_LNC */
    s_dynarray entries;      /* s_dynarray of point_entry, indexed by registry id */
};

static genericPoint &resolve(implicit_predicate_registry *reg, int id)
{
    point_entry *e = (point_entry *)dynarray_get_ptr(&reg->entries, (size_t)id);
    if (e->type == EXPLICIT3D)
        return *(genericPoint *)dynarray_get_ptr(&reg->explicit_pts, (size_t)e->exact_local_id);
    return *(genericPoint *)dynarray_get_ptr(&reg->lnc_pts, (size_t)e->exact_local_id);
}

static void ensure_entry(implicit_predicate_registry *reg, int id)
{
    static const point_entry blank = {UNDEF, -1, -1};
    while ((int)reg->entries.N <= id) dynarray_push(&reg->entries, &blank);
}

extern "C" {

void itest_init()
{
    initFPU();
}

implicit_predicate_registry *itest_registry_create()
{
    implicit_predicate_registry *reg = new implicit_predicate_registry();
    reg->explicit_pts = dynarray_initialize(sizeof(explicitPoint3D), 0);
    reg->lnc_pts      = dynarray_initialize(sizeof(implicitPoint3D_LNC), 0);
    reg->entries      = dynarray_initialize(sizeof(point_entry), 0);
    return reg;
}

void itest_registry_destroy(implicit_predicate_registry *reg)
{
    if (!reg) return;
    dynarray_free(&reg->explicit_pts);
    dynarray_free(&reg->lnc_pts);
    dynarray_free(&reg->entries);
    delete reg;
}

void itest_point_set_explicit(implicit_predicate_registry *reg, int id, s_point p)
{
    ensure_entry(reg, id);
    explicitPoint3D ep(p.x, p.y, p.z);
    int local_id = (int)reg->explicit_pts.N;
    dynarray_push(&reg->explicit_pts, &ep);
    point_entry e = {EXPLICIT3D, local_id, local_id}; /* shadow == exact: same slot */
    dynarray_change_entry(&reg->entries, (unsigned)id, &e);
}

void itest_point_set_lnc(implicit_predicate_registry *reg, int id, int v1, int v2, double t)
{
    /* paper: t*V[v1] + (1-t)*V[v2]
     * library: (1-t_lib)*P + t_lib*Q  ->  P=V[v1], Q=V[v2], t_lib = 1.0-t */
    ensure_entry(reg, id);

    point_entry *e1 = (point_entry *)dynarray_get_ptr(&reg->entries, (size_t)v1);
    point_entry *e2 = (point_entry *)dynarray_get_ptr(&reg->entries, (size_t)v2);
    /* The library's LNC combines two EXPLICIT points, so a nested Steiner --
     * one placed on a sub-segment whose endpoint is a prior Steiner -- must
     * reference that Steiner's explicit shadow, not its exact LNC. */
    int v1_shadow = e1->shadow_local_id;
    int v2_shadow = e2->shadow_local_id;

    const explicitPoint3D &vp1 = *(const explicitPoint3D *)dynarray_get_ptr(&reg->explicit_pts, (size_t)v1_shadow);
    const explicitPoint3D &vp2 = *(const explicitPoint3D *)dynarray_get_ptr(&reg->explicit_pts, (size_t)v2_shadow);

    implicitPoint3D_LNC lnc(&reg->explicit_pts, v1_shadow, v2_shadow, 1.0 - t);
    int lnc_local_id = (int)reg->lnc_pts.N;
    dynarray_push(&reg->lnc_pts, &lnc);

    /* Also store this point's rounded coords as an explicit shadow, computed
     * with the SAME t*V[v1]+(1-t)*V[v2] formula the builder uses for M
     * (points.p), so this shadow == M exactly and nested LNCs stay
     * consistent with the rounded geometry. A FUTURE LNC that references id
     * as an endpoint resolves through entries[id].shadow_local_id, never
     * through the exact LNC above. */
    explicitPoint3D shadow(t * vp1.X() + (1.0 - t) * vp2.X(),
                           t * vp1.Y() + (1.0 - t) * vp2.Y(),
                           t * vp1.Z() + (1.0 - t) * vp2.Z());
    int shadow_local_id = (int)reg->explicit_pts.N;
    dynarray_push(&reg->explicit_pts, &shadow);

    point_entry e = {LNC, lnc_local_id, shadow_local_id};
    dynarray_change_entry(&reg->entries, (unsigned)id, &e);
}

int itest_orient3d(implicit_predicate_registry *reg, int a, int b, int c, int d)
{
    return genericPoint::orient3D(resolve(reg, a), resolve(reg, b), resolve(reg, c), resolve(reg, d));
}

int itest_insphere(implicit_predicate_registry *reg, int a, int b, int c, int d, int e)
{
    return genericPoint::inSphere(resolve(reg, a), resolve(reg, b), resolve(reg, c), resolve(reg, d), resolve(reg, e));
}

int itest_perturbed_insphere(implicit_predicate_registry *reg, int i1, int i2, int i3, int i4, int i5)
{
    int r = genericPoint::inSphere(resolve(reg, i1), resolve(reg, i2), resolve(reg, i3),
                                    resolve(reg, i4), resolve(reg, i5));
    if (r != 0) return r;

    /* SoS: sort [i1..i5] ascending by insertion sort; track swap parity n. */
    int ids[5] = {i1, i2, i3, i4, i5};
    int n = 0;
    for (int i = 1; i < 5; i++) {
        int key = ids[i], j = i - 1;
        while (j >= 0 && ids[j] > key) { ids[j + 1] = ids[j]; j--; n ^= 1; }
        ids[j + 1] = key;
    }

    /* Try cofactor orient3D calls, skipping each row in turn.
     * Sign for skip k is (-1)^(k+n). */
    for (int skip = 0; skip < 5; skip++) {
        int sub[4], k = 0;
        for (int j = 0; j < 5; j++) if (j != skip) sub[k++] = ids[j];
        int o = genericPoint::orient3D(resolve(reg, sub[0]), resolve(reg, sub[1]),
                                        resolve(reg, sub[2]), resolve(reg, sub[3]));
        if (o != 0) {
            return ((skip ^ n) & 1) ? -o : o;
        }
    }
    return 0; /* all five coincident -- degenerate input */
}

int itest_point_in_tet(implicit_predicate_registry *reg, int p, int t0, int t1, int t2, int t3)
{
    /* Algorithm 3 from the paper: p is inside iff all four orient3D tests are >= 0. */
    if (genericPoint::orient3D(resolve(reg, t0), resolve(reg, t1), resolve(reg, t2), resolve(reg, p))  < 0) return 0;
    if (genericPoint::orient3D(resolve(reg, t0), resolve(reg, t1), resolve(reg, p),  resolve(reg, t3)) < 0) return 0;
    if (genericPoint::orient3D(resolve(reg, t0), resolve(reg, p),  resolve(reg, t2), resolve(reg, t3)) < 0) return 0;
    if (genericPoint::orient3D(resolve(reg, p),  resolve(reg, t1), resolve(reg, t2), resolve(reg, t3)) < 0) return 0;
    return 1;
}

int itest_segment_crosses_triangle(implicit_predicate_registry *reg, int s1, int s2, int v0, int v1, int v2)
{
    /* Algorithm 4 from the paper: interior of s1-s2 crosses interior of <v0,v1,v2>. */
    int o1 = genericPoint::orient3D(resolve(reg, v0), resolve(reg, v1), resolve(reg, v2), resolve(reg, s1));
    int o2 = genericPoint::orient3D(resolve(reg, v0), resolve(reg, v1), resolve(reg, v2), resolve(reg, s2));
    if (o1 == o2) return 0; /* both on same side (or both on plane) */
    int a = genericPoint::orient3D(resolve(reg, v0), resolve(reg, v1), resolve(reg, s1), resolve(reg, s2));
    int b = genericPoint::orient3D(resolve(reg, v1), resolve(reg, v2), resolve(reg, s1), resolve(reg, s2));
    int c = genericPoint::orient3D(resolve(reg, v2), resolve(reg, v0), resolve(reg, s1), resolve(reg, s2));
    if (a * b < 0 || b * c < 0) return 0;
    return 1;
}

int itest_inner_seg_crosses_inner_tri(implicit_predicate_registry *reg, int s1, int s2, int v0, int v1, int v2)
{
    return genericPoint::innerSegmentCrossesInnerTriangle(
        resolve(reg, s1), resolve(reg, s2), resolve(reg, v0), resolve(reg, v1), resolve(reg, v2)) ? 1 : 0;
}

int itest_inner_segs_cross(implicit_predicate_registry *reg, int a, int b, int p, int q)
{
    return genericPoint::innerSegmentsCross(
        resolve(reg, a), resolve(reg, b), resolve(reg, p), resolve(reg, q)) ? 1 : 0;
}

int itest_point_in_inner_tri(implicit_predicate_registry *reg, int p, int a, int b, int c)
{
    return genericPoint::pointInInnerTriangle(
        resolve(reg, p), resolve(reg, a), resolve(reg, b), resolve(reg, c)) ? 1 : 0;
}

int itest_point_in_triangle(implicit_predicate_registry *reg, int p, int a, int b, int c)
{
    /* OUT if not coplanar with the triangle plane. */
    if (genericPoint::orient3D(resolve(reg, a), resolve(reg, b), resolve(reg, c), resolve(reg, p)) != 0) return -1;
    /* Coplanar: pointInTriangle returns in-or-on plus the three edge signs. */
    int o1, o2, o3;
    bool on = genericPoint::pointInTriangle(resolve(reg, p), resolve(reg, a), resolve(reg, b), resolve(reg, c), o1, o2, o3);
    if (!on) return -1;
    return (o1 == 0 || o2 == 0 || o3 == 0) ? 0 : 1;
}

int itest_segments_cross(implicit_predicate_registry *reg, int s1, int s2, int p, int q)
{
    return genericPoint::segmentsCross(
        resolve(reg, s1), resolve(reg, s2), resolve(reg, p), resolve(reg, q)) ? 1 : 0;
}

int itest_get_approx_coords(implicit_predicate_registry *reg, int id, double *out)
{
    if (!reg || id < 0 || id >= (int)reg->entries.N) return 0;
    point_entry *e = (point_entry *)dynarray_get_ptr(&reg->entries, (size_t)id);
    if (e->type == UNDEF) return 0;
    return resolve(reg, id).getApproxXYZCoordinates(out[0], out[1], out[2]) ? 1 : 0;
}

} /* extern "C" */
