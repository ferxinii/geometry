#include "points.h"  
#include "gtests.h"
#include "ch_quickhull3D.h"
#include "lists.h"
#include "convh.h" // TODO DEBUG, In the future remove this.

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <errno.h> 
#include <assert.h>
#include <stdbool.h>


/* Helper structs and comparers */
typedef struct indexed_double {
    double val;
    int idx;
} s_indexed_double;

static int cmp_indexed_double_asc(const void *pa, const void *pb) {
	const s_indexed_double *a = (const s_indexed_double*)pa;
	const s_indexed_double *b = (const s_indexed_double*)pb;
	if (a->val > b->val) return 1;
	if (a->val < b->val) return  -1;  
	return 0;
}

typedef struct indexed_edge {
    int e[2];   /* Edge vertex ids, should be in canonical order (e[0] < e[1]) */
    int opp;    /* Opposite vertex id from face where edge comes from */
    int fid;    /* Face id */
} s_indexed_edge;

static int cmp_edge(const void *pa, const void *pb)
{
    const int *a = pa;
    const int *b = pb;
    if (a[0] != b[0]) return a[0] - b[0];
    return a[1] - b[1];
}

static int cmp_indexed_edge(const void *pa, const void *pb)
{
    const s_indexed_edge *a = pa;
    const s_indexed_edge *b = pb;
    return cmp_edge(a->e, b->e);
}


/* Face helpers */
static void vertices_face(const s_points *points, const int *faces, int face_id, s_point out[3])
{
    int i;
    i = faces[face_id*3+0];
    if (i>=points->N || i<0) printf("ii=%d, face_id=%d. face: %d, %d, %d\n", i, face_id, faces[face_id*3+0], faces[face_id*3+1], faces[face_id*3+2]);
    assert(i < points->N && i >= 0);
    out[0] = points->p[i];

    i = faces[face_id*3+1];
    if (i>=points->N || i<0) printf("ii=%d, face_id=%d. face: %d, %d, %d\n", i, face_id, faces[face_id*3+0], faces[face_id*3+1], faces[face_id*3+2]);
    assert(i < points->N && i >= 0);
    out[1] = points->p[i];

    i = faces[face_id*3+2];
    if (i>=points->N || i<0) printf("ii=%d, face_id=%d. face: %d, %d, %d\n", i, face_id, faces[face_id*3+0], faces[face_id*3+1], faces[face_id*3+2]);
    assert(i < points->N && i >= 0);
    out[2] = points->p[i];
}

static void flip_face(int *faces, int face_id)
{
    int tmp = faces[face_id*3 + 1];
    faces[face_id*3 + 1] = faces[face_id*3 + 2];
    faces[face_id*3 + 2] = tmp;
}

static int next_vid_isused_notinface(int Np, bool isused[Np], int face_ids[3], int current_id)
{   /* If current_id -1, return the first vertex used in the convex hull that is not in the face */
    /* Returns -1 if reached the end of the list without finding said point */
    int start = current_id + 1;
    assert(start >= 0);

    for (int ii = start; ii < Np; ++ii)
        if (isused[ii])
            if (ii != face_ids[0] && ii != face_ids[1] && ii != face_ids[2])
                return ii;
    return -1;
}

static int orient_face_if_needed(const s_points *points, bool isused[points->N], int Nfaces, int faces[3*Nfaces], int face_id)
{   /* Face is oriented if its normal points outside the convex hull (Right hand rule) */ 
    /* Returns 1 if reoriented, 0 if not, -1 if error (Could not orient it!) */
    int face_vids[3] = {faces[face_id*3+0], faces[face_id*3+1], faces[face_id*3+2]};
    s_point face_vertices[3];
    vertices_face(points, faces, face_id, face_vertices);
    
    /* Find noncoplanar point which is part of the hull and not in the current new face */
    int p = -1;
    p = next_vid_isused_notinface(points->N, isused, face_vids, p);
    if (p == -1) return -1;

    int o = orientation_robust(face_vertices, points->p[p]);
    while (o == 0) {
        p = next_vid_isused_notinface(points->N, isused, face_vids, p);
        if (p == -1) return -1;  
        o = orientation_robust(face_vertices, points->p[p]);
    }

    /* Orient faces so that each point on the original simplex can't see the opposite face */
    if (o == -1) {
        flip_face(faces, face_id);
        return 1;
    }

    return 0;
}


/* Initialise convex hull and priority of points */
static double tetrahedron_volume(const s_point a, const s_point b, const s_point c, const s_point d)
{
    s_point ab = subtract_points(b, a);
    s_point ac = subtract_points(c, a);
    s_point ad = subtract_points(d, a);
    return dot_prod(ad, cross_prod(ab, ac)) / 6.0;
}

static int find_any_non_coplanar_quad(const s_points *points, const bool *mask_dup, double EPS_degenerate, int out[4])
{   /* brute-force combinations i<j<k<l */
	int N = points->N;
    assert(N >= 4  && "Not enough points");
	for (int i=0; i<N-3; i++) if (!mask_dup[i]) {
    for (int j=i+1; j<N-2; j++) if (!mask_dup[j]) {
    for (int k=j+1; k<N-1; k++) if (!mask_dup[k]) {
        s_point tri[3] = {points->p[i], points->p[j], points->p[k]};
        for (int l=k+1; l<N; l++) if (!mask_dup[l]) {
            if (EPS_degenerate > 0) {
                double vol = tetrahedron_volume(points->p[i], points->p[j], points->p[k], points->p[l]);
                if (fabs(vol) > EPS_degenerate) { 
                    out[0]=i; out[1]=j; out[2]=k; out[3]=l;
                    return 1;
                }
            } else {
                int o = orientation_robust(tri, points->p[l]);
                if (o != 0) {
                    out[0] = i; out[1] = j; out[2] = k; out[3] = l; 
                    return 1;
                }
            }
        }
    }}}
	return 0;  /* all points are degenerate (coplanar/collinear/coincident) */
}

static int initial_tetrahedron(const s_points *points, const bool *mask_dup, bool *isused, double EPS_degenerate, int faces[4*3])
{
    /* Set the indices of the points defining the face  */
    int vertex_ids[4];
    if (!find_any_non_coplanar_quad(points, mask_dup, EPS_degenerate, vertex_ids)) return 0;

    memset(isused, 0, points->N * sizeof(bool));
    isused[vertex_ids[0]] = true;  isused[vertex_ids[1]] = true;
    isused[vertex_ids[2]] = true;  isused[vertex_ids[3]] = true;

    for (int i=0; i<4; i++)
        for (int j=0, k=0; j<4; j++)
            if (vertex_ids[j] != vertex_ids[i]) faces[i*3+k++] = vertex_ids[j];

    /* Ensure that faces are correctly oriented */
    for (int k=0; k<4; k++)
        if (orient_face_if_needed(points, isused, 4, faces, k) == -1) {
            return 0;  /* Face could not be oriented */
        }

    return 1;
}

static int priority_vertices_ignoring_initial_tetra_deduping(const s_points *points, const bool *mark_dup, const bool isused[points->N], double EPS, int out_indices[points->N-4])
{   /* Return 0 if error, Nout>0 if OK */
    /* Coordinates of the centre of the remaining point set */
    int N_aux = 0;
    s_point meanp = {0};
    for (int ii=0; ii<points->N; ii++) {
        if (!mark_dup[ii] && !isused[ii]) {
            meanp = sum_points(meanp, points->p[ii]);
            N_aux++;
        }
    }
    if (N_aux == 0) return 0;
    meanp = scale_point(meanp, 1.0/N_aux);

    /* Relative distance of points from the center */
    s_point span = span_points(points);  /* Used for normalizing, considering ALL input points */
    if (fabs(span.x) < EPS || fabs(span.y) < EPS || fabs(span.z) < EPS) return 0;  /* Points do not span the 3 dimensions. */

    /* This function is only called once, so mallocing in here and freeing is no problem */
    s_indexed_double *reldist = malloc(N_aux * sizeof(s_indexed_double));
    if (!reldist) return 0;
    for (int ii=0, jj=0; ii<points->N; ii++) {
        if (!mark_dup[ii] && !isused[ii]) {
            s_point scaled = {{{ (points->p[ii].x-meanp.x)/span.x, 
                                 (points->p[ii].y-meanp.y)/span.y, 
                                 (points->p[ii].z-meanp.z)/span.z }}};
            reldist[jj].val = norm(scaled);
            reldist[jj].idx = ii;  /* Store original index */
            jj++;
        }
    }
    
    /* Sort by relative distance ASCENDING (we will read from the end) */
    qsort(reldist, N_aux, sizeof(s_indexed_double), &cmp_indexed_double_asc);

    /* Rescue original indexing */
    for (int i=0; i<N_aux; i++) out_indices[i] = reldist[i].idx;

    free(reldist);
    return N_aux;

    // /* Check for duplicate points. IF I DONT DEDUPE, IT IS VERY LIKELY I GET ERRORS! */
    // bool *mask_dup = calloc(N_aux, sizeof(bool));
    // if (!mask_dup) { free(reldist); return 0; };
    // const double TOL2 = TOL_dup * TOL_dup;
    // for (int i=0; i<N_aux; ) {
    //     int j = i + 1;
    //     while (j<N_aux && fabs(reldist[j].val - reldist[j-1].val) <= TOL_dup) j++;
    //
    //     /* Now points [i, j) are in the same distance bucket */
    //     for (int a = i; a < j; a++) if (!mask_dup[a])
    //         for (int b = a + 1; b < j; b++) if (!mask_dup[b])
    //             if (distance_squared(points->p[reldist[a].idx], points->p[reldist[b].idx]) <= TOL2) {
    //                 mask_dup[b] = true;  /* duplicate */
    //             }
    //     i = j;
    // }
    //
    // /* Rescue original indexing of non-duplicate points */
    // int Nout = 0;
    // for (int i=0; i<N_aux; i++) if (!mask_dup[i]) out_indices[Nout++] = reldist[i].idx;
    //
    // for (int ii=0; ii<Nout-1; ii++) {  /* DEBUGGING */
    //     for (int jj=ii+1; jj<Nout; jj++) {
    //         if (distance_squared(points->p[out_indices[ii]], points->p[out_indices[jj]]) <= TOL2) {
    //         printf("distance: %.16f, TOL: %.16f\n", distance(points->p[out_indices[ii]], points->p[out_indices[jj]]), TOL_dup);
    //         printf("ii:%d, reldist:%.16f\n", ii, reldist[ii].val);
    //         printf("ii:%d, reldist:%.16f\n", jj, reldist[jj].val);
    //         exit(1);
    //         }
    //     }
    // }
    //
    // free(reldist);
    // free(mask_dup);
    // return Nout;
}


/* Visible and non visible faces, horizon */
// static int find_plane(const s_points *points, int Nfaces, int faces[Nfaces*3], s_list *cop_fids, double EPS_degenerate, s_point *out_n) 
// {   /* Iterate through coplanar faces until one is not degenerate and extract n, u, v */
//     for (int ii=0; ii<(int)cop_fids->N; ii++) {   
//         int fid; list_get_value(cop_fids, ii, &fid);
//         int f1=faces[fid*3+0], f2=faces[fid*3+1], f3=faces[fid*3+2];
//         s_point plane[3] = { points->p[f1], points->p[f2], points->p[f3] };
//         if (basis_vectors_plane(plane, EPS_degenerate, out_n, NULL, NULL)) return 1;
//     }
//     return 0;
// }



static int find_planes(const s_points *points, int Nfaces, int faces[Nfaces*3], s_list *cop_fids, double EPS_degenerate, double TOL, s_list *out_planes_n, int *out_plane_of_face)     
{   /* Finds normals of the possibly multiple planes (simplified non-general function for this use-case: planes that are coplanar with a single point) */
    out_planes_n->N = 0;

    assert(TOL >= 0);  /* TODO: Necessary? */
    for (int ii = 0; ii < (int)cop_fids->N; ii++) {
        int fid; list_get_value(cop_fids, ii, &fid);
        int f1 = faces[fid*3 + 0], f2 = faces[fid*3 + 1], f3 = faces[fid*3 + 2];
        s_point tri[3] = { points->p[f1], points->p[f2], points->p[f3] };
        s_point n;
        if (!basis_vectors_plane(tri, EPS_degenerate, &n, NULL, NULL)) {
            out_plane_of_face[ii] = -1;
            if (!list_push(out_planes_n, &point_NAN)) return 0;
            continue;
        }
        // printf("n: %g, %g, %g\n", n.x, n.y, n.z);

         /* try to match an existing plane by absolute dot-product */
        int matched = -1;
        for (int p = 0; p < (int)out_planes_n->N; p++) {
            s_point p_n; list_get_value(out_planes_n, p, &p_n);
            double dot = dot_prod(n, p_n);  
            if (fabs(dot) >= 1-TOL) { matched = p; break; }
        }

        if (matched >= 0) out_plane_of_face[ii] = matched;
        else {
            out_plane_of_face[ii] = out_planes_n->N;
            if (!list_push(out_planes_n, &n)) return 0;
        }
    }
    // printf("Nplanes: %d\n", out_planes_n->N);
    // for (int ii=0; ii<(int)out_planes_n->N; ii++) {
    //     s_point test; list_get_value(out_planes_n, ii, &test);
    //     printf("%d: %g, %g, %g\n", ii, test.x, test.y, test.z);
    // }
    return 1;
}


static void drop_to_2D(const s_point p, int coord_to_drop, double out[2])
{
    int i1 = (coord_to_drop + 1) % 3;
    int i2 = (coord_to_drop + 2) % 3;
    out[0] = p.coords[i1];
    out[1] = p.coords[i2];
}

static int coplanar_visible(const s_points *points, s_point p, int Nfaces, int faces[Nfaces*3], s_list *cop_fids, double EPS_degenerate, double TOL, s_list *buff_planes_n, s_list *buff_plane_of_face, s_list *out_indicator)
{   /* Called by mark_visible_faces_from_point, treats coplanar case. */
    // printf("COPLANAR FUNCTION: \n");;
    int N_visible = 0;

    /* Find planes */
    s_list *planes_n = buff_planes_n;
    list_ensure_capacity(buff_plane_of_face, cop_fids->N);
    int *plane_of_face = buff_plane_of_face->items;
    if (!find_planes(points, Nfaces, faces, cop_fids, EPS_degenerate, TOL, planes_n, plane_of_face)) return 0;

    /* Project to 2D */
    for (int pid=0; pid<(int)planes_n->N; pid++) {
        s_point n; list_get_value(planes_n, pid, &n);
        int coord_to_drop = coord_with_largest_component_3D(n);
        double p2D[2]; drop_to_2D(p, coord_to_drop, p2D);
        // printf("N: %g, %g, %g. Coord to drop: %d\n", n.x, n.y, n.z, coord_to_drop);
        
        /* For each plane group, build the edge list only from faces belonging to that group. */
        int faces_in_plane = 0;
        for (int ii = 0; ii < (int)cop_fids->N; ii++)
            if (plane_of_face[ii] == pid) faces_in_plane++;
        if (faces_in_plane == 0) continue;
        // printf("plane %d, Nfaces %d\n", pid, faces_in_plane);

        /* Count how many times each edge appears (canonicalized) */
        s_indexed_edge edges[faces_in_plane * 3];  // TODO allocating on stack? Is this wrong? Use list?
        int Ne = 0;
        for (int ii=0; ii<(int)cop_fids->N; ii++) {  /* Store all edges */
            if (plane_of_face[ii] != pid) continue;

            int fid; list_get_value(cop_fids, ii, &fid);
            int vA = faces[fid*3+0], vB = faces[fid*3+1];
            if (vA > vB) { int t = vA; vA = vB; vB = t; }
            edges[Ne].e[0] = vA; edges[Ne].e[1] = vB; 
            edges[Ne].opp = faces[fid*3+2]; edges[Ne].fid = fid;
            Ne++;

            vA = faces[fid*3+1], vB = faces[fid*3+2];
            if (vA > vB) { int t = vA; vA = vB; vB = t; }
            edges[Ne].e[0] = vA; edges[Ne].e[1] = vB; 
            edges[Ne].opp = faces[fid*3+0]; edges[Ne].fid = fid;
            Ne++;

            vA = faces[fid*3+0], vB = faces[fid*3+2];
            if (vA > vB) { int t = vA; vA = vB; vB = t; }
            edges[Ne].e[0] = vA; edges[Ne].e[1] = vB; 
            edges[Ne].opp = faces[fid*3+1]; edges[Ne].fid = fid;
            Ne++;
        }
        /* sort the edges */
        qsort(edges, Ne, sizeof(s_indexed_edge), cmp_indexed_edge);

        /* Determine visibility of boundary (exterior, count = 1) edges */
        int e = 0;
        while (e < Ne) {
            /* find run of identical edges */
            int run_start = e;
            int run_end = e + 1;
            while (run_end < Ne &&
                edges[run_end].e[0] == edges[run_start].e[0] &&
                edges[run_end].e[1] == edges[run_start].e[1]) {
                run_end++;
            }

            if (/* run_len */(run_end - run_start) == 1) {  /* Edge is unique -> exterior */
                // printf("unique edge %d, %d\n", edges[e].e[0], edges[e].e[1]);
                double pA[2]; drop_to_2D(points->p[edges[e].e[0]], coord_to_drop, pA);
                double pB[2]; drop_to_2D(points->p[edges[e].e[1]], coord_to_drop, pB);
                double pC[2]; drop_to_2D(points->p[edges[e].opp], coord_to_drop, pC);
                if (orient2d(pA, pB, pC) < 0) {  /* Wrong order */
                    double t[2] = {pA[0], pA[1]};
                    pA[0] = pB[0]; pA[1] = pB[1];
                    pB[0] = t[0]; pB[1] = t[1];
                }

                // printf("%g\n", orient2d(pA, pB, p2D));
                if (orient2d(pA, pB, p2D) < 0) {  /* Face is visible! */
                    // printf("%d Visible.\n", edges[e].fid);
                    int t = true;
                    list_change_entry(out_indicator, edges[e].fid, &t);
                    N_visible++;
                }
            }
            e = run_end; 
        }
    }
    return N_visible;
}


static int visible_faces_from_point(const s_points *points, int Nfaces, int faces[Nfaces*3], s_point p, double EPS, double TOL, s_list *coplanar_fids, s_list *buff_planes_n, s_list *buff_plane_of_face, s_list *out_indicator)
{   /* Returns -1 if error, Nvisible >= 0 if OK */
    /* A face is visible if its normal points to the halfspace containing the point,
       or if it is coplanar and the point sees it in the 2D POV */
    if (!list_ensure_capacity(out_indicator, Nfaces))
        { fprintf(stderr, "visible_faces_from_point: error list.\n"); return -1; }
    list_memset0(out_indicator);
    out_indicator->N = Nfaces;
    int N_visible = 0;
    coplanar_fids->N = 0; 

    for (int j = 0; j < Nfaces; ++j) {
        s_point face_pts[3]; vertices_face(points, faces, j, face_pts);
        int o = orientation_robust(face_pts, p);
        if (o < 0) {  /*  Visible, point lies on the side pointed to by face normal  (above the plane) */
            bool t = true;
            list_change_entry(out_indicator, j, &t);
            N_visible++;
        } else if (o == 0) {  /* Point is coplanar. If inside triangle, return, else add to special list. */
            // printf("COPLANAR: p=(%g, %g, %g) with face %d: %d %d %d\n", p.x, p.y, p.z, j, faces[j*3+0], faces[j*3+1], faces[j*3+2]);
            e_geom_test test = test_point_in_triangle_3D(face_pts, p, EPS, 0);
            if (test == TEST_ERROR) printf("ERROR!\n");
            if (test == TEST_BOUNDARY || test == TEST_IN) return 0;  /* If p inside any face, then immediatly exit */
            else if (test == TEST_OUT && !list_push(coplanar_fids, &j)) return 0;
        }
    }

    if (coplanar_fids->N > 0)
        N_visible += coplanar_visible(points, p, Nfaces, faces, coplanar_fids, EPS, TOL, buff_planes_n, buff_plane_of_face,  out_indicator);
    
    return N_visible;
}


static int nonvisible_faces_sharing_edge_with_face(int Nfaces, int faces[Nfaces], bool is_visible[Nfaces], int query_faceid, s_list *out_indicator)
{   /* Returns 0 if error, 1 if OK */
    if (!list_ensure_capacity(out_indicator, Nfaces))
        { fprintf(stderr, "nonvisible_faces_sharing_edge_with_face: error list.\n"); return 0; }
    list_memset0(out_indicator);
    out_indicator->N = Nfaces;

    for (int i=0; i<Nfaces; i++) if (!is_visible[i]) {
        /* Count shared vertices */
        int N_shared_vertices = 0;
        for (int a=0; a<3; a++)
            for (int b=0; b<3; b++)
                if (faces[query_faceid*3+a] == faces[i*3+b])
                    N_shared_vertices++;

        if (N_shared_vertices == 2) {
            bool t = true;
            list_change_entry(out_indicator, i, &t);
        }
        assert(N_shared_vertices < 3);  /* Non-manifold? */
    }
    return 1;
}

static int extract_horizon(int Nfaces, int faces[Nfaces], bool is_visible[Nfaces], s_list *buff_edges, s_list *nvf_sharing_edge, int *out_Nhorizon, s_list *horizon)         
{   /* Returns 0 if eror, 1 if OK.  Size of horizon: 2*out_Nhorizon */
    s_list *edges = buff_edges; 
    edges->N = 0;
    int edge_count = 0;

    for (int i=0; i<Nfaces; i++) {
        if (!is_visible[i]) continue;
        /* Mark faces that are nonvisible and share an edge with current face */
        if (!nonvisible_faces_sharing_edge_with_face(Nfaces, faces, is_visible, i, nvf_sharing_edge)) 
            { fprintf(stderr, "extract_horizon: could not mark faces.\n"); return 0; }
        for (int j=0; j<Nfaces; j++) {
            bool ind; list_get_value(nvf_sharing_edge, j, &ind);
            if (!ind) continue;

            /* But which edge ? */
            int vA = -1, vB = -1;
            for (int a=0; a<3; a++)
                for (int b=0; b<3; b++)
                    if (faces[i*3 + a] == faces[j*3 + b]) {
                        if (vA == -1) vA = faces[i*3 + a];
                        else if (vB == -1) vB = faces[i*3 + a];
                    }
            assert(vA != -1 && vB != -1 && vA != vB);

            if (vA > vB) { int t = vA; vA = vB; vB = t; }  /* canonicalize */
            if (!list_push(edges, &vA)) return 0;
            if (!list_push(edges, &vB)) return 0;
            edge_count++;
        }
    }

    /* sort & store unique the edges */
    qsort(edges->items, edge_count, sizeof(int)*2, cmp_edge);

    horizon->N = 0;
    if (!list_push(horizon, (int*)edges->items)) return 0;
    if (!list_push(horizon, (int*)edges->items + 1)) return 0;
    for (int e = 1; e < edge_count; e++) {
        int *prev = list_get_ptr(horizon, horizon->N-2);
        int *curr = list_get_ptr(edges, 2*e);
        if (curr[0] != prev[0] || curr[1] != prev[1]) {
            if (!list_push(horizon, &curr[0])) return 0;
            if (!list_push(horizon, &curr[1])) return 0;
        }
    }

    *out_Nhorizon = horizon->N/2;
    return 1;
}


/* Add and delete faces */
static void delete_visible_faces(int Nfaces, int faces[Nfaces*3], const bool faces_isvisible[Nfaces], bool *isused)
{
    for (int j=0, l=0; j<Nfaces; j++) 
        if (!faces_isvisible[j]) {
            // i.e. keep those which are non visible
            faces[l*3+0] = faces[j*3+0];
            faces[l*3+1] = faces[j*3+1];
            faces[l*3+2] = faces[j*3+2];
            l++;
        } else {
            isused[faces[j*3+0]] = false;
            isused[faces[j*3+1]] = false;
            isused[faces[j*3+2]] = false;
        }
}

static int add_faces_from_horizon(const s_points *points, bool isused[points->N], int Nfaces, s_list *faces, int Nhorizon, int horizon[Nhorizon*2], int query_pid)
{   /* Returns 0 if error, -1 if OK */
    int start = Nfaces;  /* start is the first row of the new faces */
    int N_realloc_faces = Nfaces + Nhorizon;
    if (N_realloc_faces >= CH_MAX_NUM_FACES) {
        fprintf(stderr, "add_faces_from_horizon: Max faces! (%d + %d)\n", Nfaces, Nhorizon);
        return 0; 
    }

    if (!list_ensure_capacity(faces, N_realloc_faces*3))
        { fprintf(stderr, "add_faces_from_horizon: error list.\n"); return 0; }
    for (int j=0; j<Nhorizon; j++) {
        list_change_entry(faces, Nfaces*3+0, &horizon[j*2+0]);
        list_change_entry(faces, Nfaces*3+1, &horizon[j*2+1]);
        list_change_entry(faces, Nfaces*3+2, &query_pid);

        isused[horizon[j*2+0]] = true;
        isused[horizon[j*2+1]] = true;
        Nfaces++;
    }
    
    /* Orient each new face properly */
    for (int k=start; k<Nfaces; k++)
        if (orient_face_if_needed(points, isused, Nfaces, faces->items, k) == -1) {
            int *facess = faces->items;
            printf("DEBUG: %d, %d, %d\n", facess[k*3+0], facess[k*3+1], facess[k*3+2]);
            printf("DEBUG: %g, %g, %g\n", points->p[facess[k*3+0]].x, points->p[facess[k*3+0]].y, points->p[facess[k*3+0]].z);
            printf("DEBUG: %g, %g, %g\n", points->p[facess[k*3+1]].x, points->p[facess[k*3+1]].y, points->p[facess[k*3+1]].z);
            printf("DEBUG: %g, %g, %g\n", points->p[facess[k*3+2]].x, points->p[facess[k*3+2]].y, points->p[facess[k*3+2]].z);
            printf("Rest of points:\n");
            for (int ii=0; ii<points->N; ii++) if (isused[ii]) {
                printf("%d: %g, %g, %g\n", ii, points->p[ii].x, points->p[ii].y, points->p[ii].z);
            }
            fprintf(stderr, "add_faces_from_horizon: could not orient face?\n");
            return 0;
        }

    return 1;
}


/* DEBUG */
int count_used_vertices(const s_points *p, int Nfaces, int *faces, bool *isused)
{
    (void)Nfaces; (void)faces;
    int count = 0;
    for (int ii=0; ii<p->N; ii++) if (isused[ii]) count++;
    return count;
}

static int edge_pair_cmp(const void *pa, const void *pb) 
{   /* Edge vertices must be ordered!! v[0] < v[1] */
    const int *a = pa;
    const int *b = pb;
    if (a[0] < b[0]) return -1;
    if (a[0] > b[0]) return 1;
    if (a[1] < b[1]) return -1;
    if (a[1] > b[1]) return 1;
    return 0;
}

// static int count_edges_debug(int Nfaces, int *faces, s_list *out_edges)
// {
//     out_edges->N = 0;
//     for (int f=0; f<Nfaces; f++) {
//         int a0 = faces[f*3+0], a1 = faces[f*3+1], a2 = faces[f*3+2];
//         if (!list_push(out_edges, &(int){(a0 < a1) ? a0 : a1})) return 0;
//         if (!list_push(out_edges, &(int){(a0 < a1) ? a1 : a0})) return 0;
//         if (!list_push(out_edges, &(int){(a1 < a2) ? a1 : a2})) return 0;
//         if (!list_push(out_edges, &(int){(a1 < a2) ? a2 : a1})) return 0;
//         if (!list_push(out_edges, &(int){(a2 < a0) ? a2 : a0})) return 0;
//         if (!list_push(out_edges, &(int){(a2 < a0) ? a0 : a2})) return 0;
//     }
//     qsort(out_edges->items, Nfaces*3, sizeof(int)*2, edge_pair_cmp);
//
//     /* unique in place */
//     int unique2 = 0;
//     for (int i=0; i<Nfaces*3; i++) {
//         if (unique2 == 0) {
//             list_change_entry(out_edges, unique2++, list_get_ptr(out_edges, 2*i));
//             list_change_entry(out_edges, unique2++, list_get_ptr(out_edges, 2*i+1));
//         } else {
//             int curr[2]; list_get_value(out_edges, 2*i, &curr[0]); list_get_value(out_edges, 2*i+1, &curr[1]);
//             int prev[2]; list_get_value(out_edges, 2*unique2-2, &prev[0]); list_get_value(out_edges, 2*unique2-1, &prev[1]);
//             if (curr[0] != prev[0] || curr[1] != prev[1]) {
//                 list_change_entry(out_edges, unique2++, list_get_ptr(out_edges, 2*i));
//                 list_change_entry(out_edges, unique2++, list_get_ptr(out_edges, 2*i+1));
//             }
//         }
//     }
//
//     out_edges->N = unique2;
//     return unique2/2;
// }

static int hole_exists(const s_points *p, int Nfaces, int *faces, s_list *out_edges)
{
    (void)p;
    out_edges->N = 0;
    for (int f=0; f<Nfaces; f++) {
        int a0 = faces[f*3+0], a1 = faces[f*3+1], a2 = faces[f*3+2];
        if (!list_push(out_edges, &(int){(a0 < a1) ? a0 : a1})) return 0;
        if (!list_push(out_edges, &(int){(a0 < a1) ? a1 : a0})) return 0;
        if (!list_push(out_edges, &(int){(a1 < a2) ? a1 : a2})) return 0;
        if (!list_push(out_edges, &(int){(a1 < a2) ? a2 : a1})) return 0;
        if (!list_push(out_edges, &(int){(a2 < a0) ? a2 : a0})) return 0;
        if (!list_push(out_edges, &(int){(a2 < a0) ? a0 : a2})) return 0;
    }
    qsort(out_edges->items, Nfaces*3, sizeof(int)*2, edge_pair_cmp);

    /* scan runs of identical pairs, count occurrences */
    // FILE *f = fopen("hole.txt", "w");

    int i = 0;
    int seen_boundary = 0;
    while (i < Nfaces*3) {
        int *start = (int *)list_get_ptr(out_edges, 2*i);
        if (!start) return -1;
        int run_start = i;
        int run_end = i + 1;
        while (run_end < Nfaces*3) {
            int *curr = (int *)list_get_ptr(out_edges, 2*run_end);
            if (!curr) return -1;
            if (curr[0] != start[0] || curr[1] != start[1]) break;
            run_end++;
        }
        int run_len = run_end - run_start; /* number of occurrences -> number of incident faces */
        if (run_len < 2) { seen_boundary = 1; 
            // fprintf(f, "%g, %g, %g, %d\n", p->p[start[0]].x, p->p[start[0]].y, p->p[start[0]].z, start[0]);
            // fprintf(f, "%g, %g, %g, %d\n", p->p[start[1]].x, p->p[start[1]].y, p->p[start[1]].z, start[1]);
        }
        else if (run_len > 2) {printf("Nonmanifold! edge_count > 2.\n");}

        i = run_end;
    }
    // fclose(f);

    if (seen_boundary) return 1;
    return 0;
}
//
//
// static int initialize_normals_convhull(s_convh *convh)
// {   /* 0 if error, 1 if OK. UNNORMALIZED! */
//     s_point ch_CM = point_average(&convh->points);
//     convh->fnormals = malloc(sizeof(s_point) * convh->Nf);
//     if (!convh->fnormals) return 0;
//     for (int ii = 0; ii < convh->Nf; ii++) {
//         s_point v0 = convh->points.p[convh->faces[ii * 3 + 0]];
//         s_point v1 = convh->points.p[convh->faces[ii * 3 + 1]];
//         s_point v2 = convh->points.p[convh->faces[ii * 3 + 2]];
//         s_point n = cross_prod(subtract_points(v1, v0), subtract_points(v2, v0));
//         s_point verts_face[3] = { v0, v1, v2 };
//         int o = orientation_robust(verts_face, ch_CM);
//         if (o == 0) { n = point_NAN; }
//         if (o < 0) {  // Correctly order face vertices (Right hand rule?)
//             n.x = -n.x;  n.y = -n.y;  n.z = -n.z;
//             int tmp = convh->faces[ii * 3 + 1];
//             convh->faces[ii * 3 + 1] = convh->faces[ii * 3 + 2];
//             convh->faces[ii * 3 + 2] = tmp;
//         }
//         convh->fnormals[ii] = n;
//     }
//     return 1;
// }
//
//
// static int interior_face_exists(const s_points *points, int Nfaces, int *faces)
// {
//     s_convh ch;
//     ch.Nf = Nfaces;
//     ch.faces = faces;
//     ch.points = *points;
//     if (!initialize_normals_convhull(&ch)) { puts("ERROR NORMALS!"); exit(1); }
//
//     /* Check if face center of mass is on boundary of convhull? */
//     for (int ii=0; ii<Nfaces; ii++) {
//         s_point p0 = points->p[faces[ii*3+0]];
//         s_point p1 = points->p[faces[ii*3+1]];
//         s_point p2 = points->p[faces[ii*3+2]];
//         s_point CM;
//         CM.x = (p0.x + p1.x + p2.x) / 3.0;
//         CM.y = (p0.y + p1.y + p2.y) / 3.0;
//         CM.z = (p0.z + p1.z + p2.z) / 3.0;
//
//         e_geom_test test = test_point_in_convhull(&ch, CM, 1e-12, 1e-8);
//         if (test == TEST_IN) {
//             s_point closest = closest_point_on_convhull_boundary(&ch, CM, 1e-12);
//             printf("distance: %g\n", distance(closest, CM));
//             free(ch.fnormals);
//             return 1;
//         }
//     }
//     free(ch.fnormals);
//     return 0;
// }
//
//
// static int debug_ensure_face_orientation(const s_points *p, bool *isused, int Nf, int *f)
// {
//     int Nused = 0;
//     s_point CM = {0};
//     for (int ii=0; ii<p->N; ii++) if (isused[ii]) {
//         Nused++;
//         CM = sum_points(CM, p->p[ii]);
//     }
//     CM = scale_point(CM, 1.0/Nused);
//
//     int out = 1;
//     for (int ii=0; ii<Nf; ii++) {
//         s_point face_vertices[3]; vertices_face(p, f, ii, face_vertices);
//         
//         int o = orientation_robust(face_vertices, CM);
//         if (o < 1) {
//             printf("face %d / %d not oriented\n", ii, Nf);
//             out = 0;
//         }
//     }
//
//     FILE *file = fopen("CM.txt", "w");
//     fprintf(file, "%g, %g, %g", CM.x, CM.y, CM.z);
//     fclose(file);
//     return out;
// }

/* Main algorithm */
int quickhull_3d(const s_points *in_vertices, double EPS_degenerate, double TOL_dup, bool buff_isused[in_vertices->N], int **out_faces, int *N_out_faces) 
{   /* Returns:
       -1: computational error (memory, reached max_faces, ...)
       0: degeneracy error
       1: output hull is exact (All faces are big enough)
    */
    if (!in_vertices || !buff_isused || !out_faces || !N_out_faces) return -1;
    *out_faces = NULL;  *N_out_faces = 0;
    s_list faces = {0}, faces_isvisible = {0}, horizon = {0}, AUX_nvf_sharing_edge = {0}, AUX_edges = {0}, AUX_cop_fids = {0}, AUX_planes_n = {0}, AUX_plane_of_face = {0};
    int *pleft = NULL;
    bool *mark_dup = NULL;
    int Nfaces = 0;

    mark_dup = malloc(sizeof(bool) * in_vertices->N);
    if (!mark_dup) { fprintf(stderr, "ch_quickhull3D: error malloc.\n"); goto error; }
    int Ndup = mark_duplicate_points(in_vertices, TOL_dup, mark_dup);
    if (in_vertices->N - Ndup <= 3) { fprintf(stderr, "ch_quickhull3D: N <= 3.\n"); goto error_degenerate; }


    /* The initial convex hull is a tetrahedron with 4 faces (simplex) */
    Nfaces = 4;
    faces = list_initialize(sizeof(int), 12);
    if (!faces.items) { fprintf(stderr, "ch_quickhull3D: error malloc.\n"); goto error; }
    if (!initial_tetrahedron(in_vertices, mark_dup, buff_isused, EPS_degenerate, faces.items)) 
        { fprintf(stderr, "ch_quickhull3D: error initial tetrahedron.\n"); goto error_degenerate; }
    if (in_vertices->N == 4) {
        *out_faces = faces.items;
        *N_out_faces = 4;
        return 1;
    }
    

    /* Initialize the vector of points left. The points with the larger relative
     distance from the center are scanned first, which are in the END. */
    pleft = malloc((in_vertices->N - 4) * sizeof(int));  /* Worst case: no duplicates */
    if (!pleft) { fprintf(stderr, "ch_quickhull3D: error malloc.\n"); goto error; }
    int N_pleft = priority_vertices_ignoring_initial_tetra_deduping(in_vertices, mark_dup, buff_isused, EPS_degenerate, pleft);
    if (N_pleft == 0) { fprintf(stderr, "ch_quickhull3D: error priority vertices.\n"); goto error; }


    // s_list edges = list_initialize(sizeof(int), 0);  /* DEBUG */
    /* The main loop for the quickhull algorithm */
    faces_isvisible = list_initialize(sizeof(bool), 0);
    AUX_nvf_sharing_edge = list_initialize(sizeof(bool), 0); 
    horizon = list_initialize(sizeof(int), 0);
    AUX_edges = list_initialize(sizeof(int), 0);  
    AUX_cop_fids = list_initialize(sizeof(int), 0);  
    AUX_planes_n = list_initialize(sizeof(s_point), 3);
    AUX_plane_of_face = list_initialize(sizeof(int), 0);
    if (!faces_isvisible.items || !AUX_nvf_sharing_edge.items || !horizon.items || 
        !AUX_edges.items || !AUX_cop_fids.items || !AUX_planes_n.items || !AUX_plane_of_face.items) 
        { fprintf(stderr, "ch_quickhull3D: error malloc.\n"); goto error; }
    while (N_pleft > 0) {
        // s_convh ch;
        // ch.Nf = Nfaces;
        // ch.faces = faces.items;
        // ch.points = *in_vertices;
        // write_convhull_to_m(&ch, "prev.m");

        /* Process the LAST element of pleft */
        N_pleft--;
        int current_id = pleft[N_pleft];
        s_point current_p = in_vertices->p[current_id];

        // for (int ii=0; ii<in_vertices->N; ii++) {
        //     if (buff_isused[ii]) {
        //         assert(distance(current_p, in_vertices->p[ii]) >= TOL_dup);
        //     }
        // }


        /* Mark visible faces from this point */
        int N_vf = visible_faces_from_point(in_vertices, Nfaces, faces.items, current_p, EPS_degenerate, TOL_dup, &AUX_cop_fids, &AUX_planes_n, &AUX_plane_of_face, &faces_isvisible);
        if (N_vf == -1) { fprintf(stderr, "ch_quickhull3D: error visible faces.\n"); goto error; }

        // FILE *f = fopen("visible.txt", "w");
        // for (int ii=0; ii<Nfaces; ii++) {
        //     if (((bool*)faces_isvisible.items)[ii]) {
        //         int *facess = faces.items;
        //         // fprintf(f, "%d, %d, %d\n", facess[ii*3+0], facess[ii*3+1], facess[ii*3+2]);
        //         printf("visible: %d, %d, %d\n", facess[ii*3+0], facess[ii*3+1], facess[ii*3+2]);
        //     }
        // }
        // fclose(f);

        /* Proceed if N_visible_faces > 0 */
        if (N_vf == Nfaces) { fprintf(stderr, "quickhull_3d: Point sees all faces.\n"); goto error; }
        if (N_vf == 0) continue;  

        /* Create horizon */
        int Nhorizon;
        if (!extract_horizon(Nfaces, faces.items, faces_isvisible.items, &AUX_edges, &AUX_nvf_sharing_edge, &Nhorizon, &horizon))
            { fprintf(stderr, "quickhull_3d: Error extracting horizon.\n"); goto error; }

        // f = fopen("horizon.txt", "w");
        // for (int ii=0; ii<Nhorizon; ii++) {
        //     int *horizonn = horizon.items;
        //     fprintf(f, "%d, %d\n", horizonn[ii*2+0], horizonn[ii*2+1]);
        // }
        // fclose(f);

        buff_isused[current_id] = true;

        /* Delete visible faces */
        delete_visible_faces(Nfaces, faces.items, faces_isvisible.items, buff_isused);
        Nfaces -= N_vf;
        
        /* Add faces connecting horizon to the new point */
        if(!add_faces_from_horizon(in_vertices, buff_isused, Nfaces, &faces, Nhorizon, horizon.items, current_id))
            { fprintf(stderr, "quickhull_3d: Error adding new faces.\n"); goto error; }

        Nfaces += Nhorizon;

        // for (int ii=0; ii<Nfaces; ii++) {
        //     // if (((bool*)faces_isvisible.items)[ii]) {
        //         int *facess = faces.items;
        //         // fprintf(f, "%d, %d, %d\n", facess[ii*3+0], facess[ii*3+1], facess[ii*3+2]);
        //         printf("face %d: %d, %d, %d\n", ii, facess[ii*3+0], facess[ii*3+1], facess[ii*3+2]);
        //     // }
        // }
        // exit(1);

        /* DEBUG */
    //     int Nv = count_used_vertices(in_vertices, Nfaces, faces.items, buff_isused);
    //     int Nedges = count_edges_debug(Nfaces, faces.items, &edges);
    //     int euler = Nv - Nedges + Nfaces;
    //     if (hole_exists(in_vertices, Nfaces, faces.items, &edges) ||         debug_ensure_face_orientation(in_vertices, buff_isused, Nfaces, faces.items) != 1) {
    //         printf("Just added: %d, HOLE!\n", current_id);
    //         s_convh ch;
    //         ch.Nf = Nfaces;
    //         ch.faces = faces.items;
    //         ch.points = *in_vertices;
    //         write_convhull_to_m(&ch, "hole.m");
    //
    //         FILE *f = fopen("visible.txt", "w");
    //         for (int ii=0; ii<Nfaces; ii++) {
    //             if (((bool*)faces_isvisible.items)[ii]) {
    //                 int *facess = faces.items;
    //                 fprintf(f, "%d, %d, %d\n", facess[ii*3+0], facess[ii*3+1], facess[ii*3+2]);
    //             }
    //         }
    //         fclose(f);
    //
    //         puts("HORIZON:");
    //         for (int ii=0; ii<Nhorizon; ii++) {
    //             int e0; list_get_value(&horizon, 2*ii, &e0);
    //             int e1; list_get_value(&horizon, 2*ii+1, &e1);
    //             printf("(%d, %d)\n", e0, e1);
    //         }
    //
    //         exit(1);
    //     }
    //     if (euler != 2) {
    //         fprintf(stderr, "EULER=%d, N_pleft=%d, hole_exists=%d, interior_face_exists=%d\n", euler, N_pleft, hole_exists(in_vertices, Nfaces, faces.items, &edges), interior_face_exists(in_vertices, Nfaces, faces.items));
    //         fprintf(stderr, "dV = %d, dF = %d, dE = %d\n", 1, Nhorizon - N_vf, Nhorizon);
    //         // assert(euler == 2);
    //
    //         // puts("HORIZON:");
    //         // for (int ii=0; ii<Nhorizon; ii++) {
    //         //     int e0; list_get_value(&horizon, 2*ii, &e0);
    //         //     int e1; list_get_value(&horizon, 2*ii+1, &e1);
    //         //     printf("(%d, %d)\n", e0, e1);
    //         // }
    //         // exit(1);
    //     }
    }
    /* Debugging */
    s_list edges = list_initialize(sizeof(int), 0);
    if (!edges.items) goto error;
    if (hole_exists(in_vertices, Nfaces, faces.items, &edges)) 
        { fprintf(stderr, "ch_quickhull3D: hole exists.\n"); free_list(&edges); goto error; }
    free_list(&edges); 

    // int Nv = count_used_vertices(in_vertices, Nfaces, faces.items, buff_isused);
    // int Nedges = count_edges_debug(Nfaces, faces.items, &edges);
    // int euler = Nv - Nedges + Nfaces;
    // if (euler != 2) 
    //     { fprintf(stderr, "ch_quickhull3D: euler = %d.\n", euler); goto error; }

    

    /* clean-up and exit */
    free(pleft);
    free(mark_dup);
    free_list(&faces_isvisible); 
    free_list(&horizon);
    free_list(&AUX_nvf_sharing_edge);
    free_list(&AUX_planes_n);
    free_list(&AUX_plane_of_face);
    *out_faces = realloc(faces.items, Nfaces * 3 * sizeof(int));
    *N_out_faces = Nfaces;
    return 1;
    
    error_degenerate:
        fprintf(stderr, "Error in 'quickhull_3d'. In vertices are degenerate?\n");
        if (mark_dup) free(mark_dup);
        if (faces.items) free_list(&faces);
        if (pleft) free(pleft);
        if (faces_isvisible.items) free_list(&faces_isvisible); 
        if (horizon.items) free_list(&horizon);
        if (AUX_nvf_sharing_edge.items) free_list(&AUX_nvf_sharing_edge);
        if (AUX_planes_n.items) free_list(&AUX_planes_n);
        if (AUX_plane_of_face.N) free_list(&AUX_plane_of_face);
        *out_faces = NULL;
        *N_out_faces = 0;
        return 0;
        
    error:
        fprintf(stderr, "Error in 'quickhull_3d'. Faces N = %d / %d\n", Nfaces, CH_MAX_NUM_FACES);
        if (mark_dup) free(mark_dup);
        if (faces.items) free_list(&faces);
        if (pleft) free(pleft);
        if (faces_isvisible.items) free_list(&faces_isvisible); 
        if (horizon.items) free_list(&horizon);
        if (AUX_nvf_sharing_edge.items) free_list(&AUX_nvf_sharing_edge);
        if (AUX_planes_n.items) free_list(&AUX_planes_n);
        if (AUX_plane_of_face.N) free_list(&AUX_plane_of_face);
        *out_faces = NULL;
        *N_out_faces = 0;
        return -1;
}

