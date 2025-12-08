#include "points.h"
#include "gtests.h"
#include "convh.h"
#include "ch_intersect.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>


typedef struct point_list {
    s_point *list;
    int Nmax;
} s_point_list;

static s_point_list initialize_point_list(int Nmax) {
    if (Nmax <= 0) Nmax = CH_N_INIT_POINT_LIST;
    return (s_point_list) { .Nmax = Nmax,
                            .list = malloc(Nmax * sizeof(s_point)) };
}

static int increase_memory_point_list(s_point_list *point_list, int N_needed)
{
    while (N_needed >= point_list->Nmax) {
        s_point *tmp = realloc(point_list->list, 2 * point_list->Nmax * sizeof(s_point));
        if (!tmp) return 0;
        point_list->list = tmp;
        point_list->Nmax *= 2;
    }
    return 1;
}

static void free_point_list(s_point_list *point_list)
{
    free(point_list->list);
    memset(point_list, 0, sizeof(s_point_list));
}


/* Edges */
typedef struct edge {
    int v[2];
} s_edge;

static int edge_pair_cmp(const void *pa, const void *pb) 
{   /* Edge vertices must be ordered!! v[0] < v[1] */
    const s_edge *a = pa;
    const s_edge *b = pb;
    if (a->v[0] < b->v[0]) return -1;
    if (a->v[0] > b->v[0]) return 1;
    if (a->v[1] < b->v[1]) return -1;
    if (a->v[1] > b->v[1]) return 1;
    return 0;
}

typedef struct edge_list {
    s_edge *list;
    int Nmax;
} s_edge_list;

static s_edge_list initialize_edge_list(int Nmax) {
    if (Nmax <= 0) Nmax = CH_N_INIT_POINT_LIST;
    return (s_edge_list) { .Nmax = Nmax,
                           .list = malloc(Nmax * sizeof(s_edge)) };
}

static int increase_memory_edge_list(s_edge_list *edge_list, int N_needed)
{
    while (N_needed >= edge_list->Nmax) {
        s_edge *tmp = realloc(edge_list->list, 2 * edge_list->Nmax * sizeof(s_edge));
        if (!tmp) return 0;
        edge_list->list = tmp;
        edge_list->Nmax *= 2;
    }
    return 1;
}

static void free_edge_list(s_edge_list *edge_list)
{
    free(edge_list->list);
    memset(edge_list, 0, sizeof(s_edge_list));
}

static void list_edges_convhull(const s_convh *C, int *out_Nedges, s_edge_list *out_edges)
{
    increase_memory_edge_list(out_edges, 3 * C->Nf);

    int wrote = 0;
    for (int f=0; f<C->Nf; f++) {
        int a0 = C->faces[f*3+0], a1 = C->faces[f*3+1], a2 = C->faces[f*3+2];

        s_edge e;
        e.v[0] = (a0 < a1) ? a0 : a1;
        e.v[1] = (a0 < a1) ? a1 : a0;
        out_edges->list[wrote++] = e;

        e.v[0] = (a1 < a2) ? a1 : a2;
        e.v[1] = (a1 < a2) ? a2 : a1;
        out_edges->list[wrote++] = e;

        e.v[0] = (a2 < a0) ? a2 : a0;
        e.v[1] = (a2 < a0) ? a0 : a2;
        out_edges->list[wrote++] = e;
    }

    qsort(out_edges->list, wrote, sizeof(s_edge), edge_pair_cmp);

    /* unique in place */
    int uniq = 0;
    for (int i=0; i<wrote; i++) {
        if (uniq == 0) out_edges->list[uniq++] = out_edges->list[i];
        else if (out_edges->list[i].v[0] != out_edges->list[uniq-1].v[0] ||
                 out_edges->list[i].v[1] != out_edges->list[uniq-1].v[1]) 
            out_edges->list[uniq++] = out_edges->list[i];
    }

    *out_Nedges = uniq;
}



/* Clip segment with convhull */
static void append_if_nonexistent_lim2(int *Nout, s_point out[2], s_point p, double TOL)
{
    if (*Nout >= 2) return;
    for (int ii=0; ii<*Nout; ii++) {
        if (distance_squared(out[ii], p) <= TOL*TOL) return;
    }
    out[(*Nout)++] = p;
}

static e_intersect_type core_segment_convhull_surface_intersect(const s_convh *C, const s_point segment[2], double EPS_degenerate, double TOL, int *Nout, s_point out[2])
{
    e_geom_test i0 = test_point_in_convhull(C, segment[0], EPS_degenerate, TOL);
    e_geom_test i1 = test_point_in_convhull(C, segment[1], EPS_degenerate, TOL);

    if (i0 == TEST_IN && i1 == TEST_IN) { 
        if (Nout) *Nout = 0;
        return INTERSECT_EMPTY; 
    };
    if (i0 == TEST_BOUNDARY && i1 == TEST_BOUNDARY) { 
        if (Nout && out) { out[0] = segment[0]; out[1] = segment[1]; *Nout = 2; }
        return INTERSECT_DEGENERATE;
    }
    if (i0 == TEST_BOUNDARY || i1 == TEST_BOUNDARY) {
        if (Nout && out) { out[0] = (i0 == TEST_BOUNDARY) ? segment[0] : segment[1]; *Nout = 1; }
        return INTERSECT_DEGENERATE;
    }


    /* Intersect segment with all faces of the convex hull */
    int h = 0;
    s_point tmp[2];
    for (int ii=0; ii<C->Nf; ii++) {
        if (h == 2) break;
        s_point face[3];
        convh_get_face(C, ii, face);
        
        int o0 = orientation_robust(face, segment[0]);
        int o1 = orientation_robust(face, segment[1]);
        if (o0 != 0 && o1 != 0 && o0 == o1) continue;  /* Skip if on same side of plane */

        s_segment_intersect fi = segment_triangle_intersect_3D(segment, face, EPS_degenerate, TOL);
        for (int jj=0; jj<fi.N; jj++) append_if_nonexistent_lim2(&h, tmp, fi.coords[jj], TOL);
    }

    if (h == 0) {
        if (Nout) *Nout = 0;
        return INTERSECT_EMPTY;
    }
    if (h == 1) { 
        if (Nout && out) { out[0] = tmp[0]; *Nout = 1; } 
        return INTERSECT_NONDEGENERATE;
    }
    if (h == 2) { 
        if (Nout && out) { out[0] = tmp[0]; out[1] = tmp[1]; *Nout = 2; }
        return INTERSECT_NONDEGENERATE;
    }
    fprintf(stderr, "intersect_segment_convhull_surface: Should never reach this!\n");
    exit(1);
}

e_intersect_type test_segment_convhull_surface_intersect(const s_convh *C, const s_point segment[2], double EPS_degenerate, double TOL)
{
    return core_segment_convhull_surface_intersect(C, segment, EPS_degenerate, TOL, NULL, NULL);
}

s_segment_intersect segment_convhull_surface_intersect(const s_convh *C, const s_point segment[2], double EPS_degenerate, double TOL)
{
    s_segment_intersect out;
    out.type = core_segment_convhull_surface_intersect(C, segment, EPS_degenerate, TOL, &out.N, out.coords);
    return out;
}


/* Intersect two convhulls */
int intersection_convhulls(const s_convh *A, const s_convh *B, double EPS_degenerate, double TOL, s_convh *out)
{   /* Returns  1 if intersection non-empty,
                0 if intersection empty,
                -1 if error or intersection empty 
    */
    e_geom_test *buff = malloc((int)fmax(A->points.N, B->points.N) * sizeof(e_geom_test));
    s_edge_list elist = initialize_edge_list((int)fmax(3*A->Nf, 3*B->Nf));
    s_point_list pI = initialize_point_list(0);
    if (!buff || !elist.list || !pI.list) goto error;

    /* First add points of A nonstrictly inside B */
    s_points_test AinB = test_points_in_convhull(B, &A->points, EPS_degenerate, 0, buff);
    if (AinB.Nerr > 0) goto error;
    int NpI = 0;
    for (int ii=0; ii<A->points.N; ii++) 
        if (AinB.indicator[ii] == TEST_IN || AinB.indicator[ii] == TEST_BOUNDARY) {
            if (!increase_memory_point_list(&pI, NpI+1)) goto error;
            pI.list[NpI++] = A->points.p[ii];
        }

    /* Intersect all edges of A with B */
    int N_eA;
    list_edges_convhull(A, &N_eA, &elist);
    for (int ii=0; ii<N_eA; ii++) {
        s_point segment[2] = {A->points.p[elist.list[ii].v[0]], A->points.p[elist.list[ii].v[1]]}; 
        s_segment_intersect clips = segment_convhull_surface_intersect(B, segment, EPS_degenerate, TOL);
        for (int jj=0; jj<clips.N; jj++) {
            if (!increase_memory_point_list(&pI, NpI+1)) goto error;
            pI.list[NpI++] = clips.coords[jj];
        }
    }

    /* Add points of B nonstrictly inside A */
    s_points_test BinA = test_points_in_convhull(A, &B->points, EPS_degenerate, 0, buff);
    if (AinB.Nerr > 0) goto error;
    for (int ii=0; ii<B->points.N; ii++) 
        if (BinA.indicator[ii] == TEST_IN || BinA.indicator[ii] == TEST_BOUNDARY) {
            if (!increase_memory_point_list(&pI, NpI+1)) goto error;
            pI.list[NpI++] = B->points.p[ii];
        }

    /* Intersect all edges of B with A */
    int N_eB;
    list_edges_convhull(B, &N_eB, &elist);
    for (int ii=0; ii<N_eB; ii++) {
        s_point segment[2] = {B->points.p[elist.list[ii].v[0]], B->points.p[elist.list[ii].v[1]]}; 
        s_segment_intersect clips = segment_convhull_surface_intersect(A, segment, EPS_degenerate, TOL);
        for (int jj=0; jj<clips.N; jj++) {
            if (!increase_memory_point_list(&pI, NpI+1)) goto error;
            pI.list[NpI++] = clips.coords[jj];
        }
    }


    int i = convhull_from_points(&(s_points){NpI, pI.list}, EPS_degenerate, TOL, out);
    if (i == -1) goto error;
    if (i == 0) {  /* Could not initialize I. It is empty / degenerate. */
        free(buff);
        free_edge_list(&elist);
        free_point_list(&pI);
        *out = convhull_NAN;
        return 0;
    }
    if (!convhull_is_valid(out)) goto error;  /* Should not happen, but sanity check */
    free(buff);
    free_edge_list(&elist);
    free_point_list(&pI);
    return 1;

    error:
        free(buff);
        free_edge_list(&elist);
        free_point_list(&pI);
        *out = convhull_NAN;
        return -1;
}


static void bisector_plane(const s_convh *A, const s_convh *B, const s_convh *I, double EPS_degenerate, s_point out[3])
{   /* Normal points from A to B */  // TODO handle degeneracies
    s_point cA = convhull_volume_centroid(A, EPS_degenerate);
    s_point cB = convhull_volume_centroid(B, EPS_degenerate);
    s_point n = normalize_vec(subtract_points(cB, cA), EPS_degenerate);
    
    double vA = volume_convhull(A), vB = volume_convhull(B);

    double w = vA / (vA + vB);
    s_point candidate = sum_points(cA, scale_point(subtract_points(cB, cA), w));
    
    // Clamp so plane passes through or touches intersection
    double tmin =  dot_prod(I->points.p[0], n), tmax = dot_prod(I->points.p[0], n);
    for (int i=1; i<I->points.N; i++) {
        double t = dot_prod(I->points.p[i], n);
        if (t < tmin) tmin = t;
        if (t > tmax) tmax = t;
    }
    double tplane = dot_prod(candidate, n);
    if (tplane < tmin) tplane = tmin;
    if (tplane > tmax) tplane = tmax;

    s_point c_plane = scale_point(n, tplane);

    // Orthonormal basis for plane
    int ref_coord = coord_with_smallest_component_3D(n);
    s_point ref = (ref_coord == 0) ?   (s_point){{{1,0,0}}} :
                  ( (ref_coord == 1) ? (s_point){{{0,1,0}}} :
                                       (s_point){{{0,0,1}}} );
    s_point u = normalize_vec(cross_prod(n, ref), EPS_degenerate);
    s_point v = normalize_vec(cross_prod(n, u), EPS_degenerate);

    out[0] = c_plane;
    out[1] = sum_points(c_plane, u);
    out[2] = sum_points(c_plane, v);
}


int clip_convhull_halfspace(const s_convh *C, s_point plane[3], double EPS_degenerate, double TOL, s_convh *out) 
{   /* Returns  1 if clipped
                0 if not clipped
                -1 if error
    */
    s_points_test ptest = {0};  s_point_list p = {0};  s_edge_list elist = {0};
    /* Select points non-strictly inside halfspace */
    ptest = test_points_in_halfspace(plane, &C->points, EPS_degenerate, TOL, NULL);
    if (!ptest.indicator) goto error;
    if (ptest.Nerr > 0) { *out = convhull_NAN; return -1; }
    if (ptest.Nin + ptest.Nbdy == 0) { *out = convhull_NAN; return 0; }

    p = initialize_point_list(ptest.Nin + ptest.Nbdy);
    if (!p.list) goto error;
    int Np = 0;
    for (int ii=0; ii<C->points.N; ii++)
        if (ptest.indicator[ii] == TEST_IN || ptest.indicator[ii] == TEST_BOUNDARY) p.list[Np++] = C->points.p[ii];
    
    /* Clip edges that cross the plane */
    int Nedges;
    elist = initialize_edge_list(3 * C->Nf);
    if (!elist.list) goto error;
    list_edges_convhull(C, &Nedges, &elist);
    for (int ii=0; ii<Nedges; ii++) {
        if (ptest.indicator[elist.list[ii].v[0]] != TEST_ERROR && 
            ptest.indicator[elist.list[ii].v[1]] != TEST_ERROR &&
            ptest.indicator[elist.list[ii].v[0]] != ptest.indicator[elist.list[ii].v[1]]) {  /* Boundary points already considered */
            s_point segment[2] = {C->points.p[elist.list[ii].v[0]], C->points.p[elist.list[ii].v[1]]};
            s_segment_intersect intersections = segment_plane_intersect(segment, plane, EPS_degenerate, TOL);
            for (int jj=0; jj<intersections.N; jj++) {
                if (!increase_memory_point_list(&p, Np+1)) goto error;
                p.list[Np++] = intersections.coords[jj];
            }
        }
    }
    
    int i = convhull_from_points(&(s_points){Np, p.list}, EPS_degenerate, TOL, out);
    if (i == -1) goto error;
    if (i == 0) {
        free_edge_list(&elist);
        free_point_list(&p);
        free(ptest.indicator);
        return 0;
    }
    if (!convhull_is_valid(out)) goto error;
    free_edge_list(&elist);
    free_point_list(&p);
    free(ptest.indicator);
    return 1;

    error:
        if (p.list) free_point_list(&p);
        if (elist.list) free_edge_list(&elist);
        if (ptest.indicator) free(ptest.indicator);
        *out = convhull_NAN;
        return -1;
}


int remove_intersection_convhulls(s_convh *A, s_convh *B, double EPS_degenerate, double TOL, double min_vol_I)
{   /* TODO: What if one is completely contained in the other? */
    /* Returns  1 if OK,
     *          0 if no intersection to remove,
     *          -1 if error
     */
    s_convh I = convhull_NAN, IA = convhull_NAN, IB = convhull_NAN; 
    s_points p_newA = {0}, p_newB = {0};
    e_geom_test *buff = NULL;

    int i = intersection_convhulls(A, B, EPS_degenerate, TOL, &I);
    if (i == 0) return 0;
    if (i == -1) return -1;
    if (volume_convhull(&I) < min_vol_I) { free_convhull(&I); return 0; }

    buff = malloc((int)fmax(A->points.N, B->points.N) * sizeof(e_geom_test));
    if (!buff) goto error;

    s_point plane[3];
    bisector_plane(A, B, &I, EPS_degenerate, plane);  /* Plane normal from A to B*/


    /* Points for reduced A */
    i = clip_convhull_halfspace(&I, plane, EPS_degenerate, TOL, &IA);  /* Divide I with the plane */
    if (i == -1) goto error;
    s_points_test Atest = test_points_in_halfspace(plane, &A->points, EPS_degenerate, TOL, buff);
    p_newA.N = Atest.Nin + Atest.Nbdy + IA.points.N;
    assert(IA.points.N > 0);
    p_newA.p = malloc(sizeof(s_point) * p_newA.N);
    if (!p_newA.p) goto error;
    int jj=0;
    for (int ii=0; ii<A->points.N; ii++)
        if (Atest.indicator[ii] == TEST_IN || Atest.indicator[ii] == TEST_BOUNDARY) p_newA.p[jj++] = A->points.p[ii];  /* Keep points of A nonstriclty inside plane */
    for (int ii=0; ii<IA.points.N; ii++) 
        p_newA.p[jj++] = IA.points.p[ii];


    /* Points for reduced B */
    s_point tmp = plane[0]; plane[0] = plane[1];  plane[1] = tmp;  /* Flip plane normal */
    i = clip_convhull_halfspace(&I, plane, EPS_degenerate, TOL, &IB);
    if (i == -1) goto error;
    s_points_test Btest = test_points_in_halfspace(plane, &B->points, EPS_degenerate, TOL, buff);
    p_newB.N = Btest.Nin + Btest.Nbdy + IB.points.N;
    assert(p_newB.N > 0);
    p_newB.p = malloc(sizeof(s_point) * p_newB.N);
    if (!p_newB.p) goto error;
    jj=0;
    for (int ii=0; ii<B->points.N; ii++)
        if (Btest.indicator[ii] == TEST_IN || Btest.indicator[ii] == TEST_BOUNDARY) p_newB.p[jj++] = B->points.p[ii];  /* Keep points of B OUTSIDE A */
    for (int ii=0; ii<IB.points.N; ii++)
        p_newB.p[jj++] = IB.points.p[ii];
    

    /* Construct new hulls and clean-up */
    s_convh newA, newB;
    if (convhull_from_points(&p_newA, EPS_degenerate, TOL, &newA) != 1) goto error;
    if (convhull_from_points(&p_newB, EPS_degenerate, TOL, &newB) != 1) goto error;

    free_convhull(A);
    free_convhull(B);
    *A = newA;
    *B = newB;

    free(buff);
    free_points(&p_newA);
    free_points(&p_newB);
    free_convhull(&I);
    free_convhull(&IA);
    free_convhull(&IB);

    return 1;

    error:
        if (buff) free(buff);
        if (points_is_valid(&p_newA)) free_points(&p_newA);
        if (points_is_valid(&p_newB)) free_points(&p_newB);
        if (convhull_is_valid(&I)) free_convhull(&I);
        if (convhull_is_valid(&IA)) free_convhull(&IA);
        if (convhull_is_valid(&IB)) free_convhull(&IB);
        return -1;
}


