#include "points.h"  
#include "gtests.h"
#include "ch_quickhull3D.h"
#include "lists.h"

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
    assert(i < points->N && i >= 0);
    out[0] = points->p[i];

    i = faces[face_id*3+1];
    assert(i < points->N && i >= 0);
    out[1] = points->p[i];

    i = faces[face_id*3+2];
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
    if (o < 0) {
        flip_face(faces, face_id);
        return 1;
    }

    return 0;
}


/* Initialise convex hull and priority of points */
static int find_any_non_coplanar_quad(const s_points *points, int out[4])
{   /* brute-force combinations i<j<k<l */
	int N = points->N;
    assert(N >= 4  && "Not enough points");
	for (int i=0; i<N-3; i++)
    for (int j=i+1; j<N-2; j++)
    for (int k=j+1; k<N-1; k++) {
        s_point tri[3] = {points->p[i], points->p[j], points->p[k]};
        for (int l=k+1; l<N; l++) {
            int o = orientation_robust(tri, points->p[l]);
            if (o != 0) {
                out[0] = i;
                out[1] = j;
                out[2] = k;
                out[3] = l;
                return 1;
            }
        }
    }
	return 0;  /* all points are degenerate (coplanar/collinear/coincident) */
}

static int initial_tetrahedron(const s_points *points, bool isused[points->N], int faces[4*3])
{
    /* Set the indices of the points defining the face  */
    int vertex_ids[4];
    if (!find_any_non_coplanar_quad(points, vertex_ids)) return 0;

    memset(isused, 0, points->N * sizeof(bool));
    isused[vertex_ids[0]] = true;  isused[vertex_ids[1]] = true;
    isused[vertex_ids[2]] = true;  isused[vertex_ids[3]] = true;

    for (int i=0; i<4; i++)
        for (int j=0, k=0; j<4; j++)
            if (vertex_ids[j] != vertex_ids[i]) faces[i*3+k++] = vertex_ids[j];

    /* Ensure that faces are correctly oriented */
    for (int k=0; k<4; k++)
        if (orient_face_if_needed(points, isused, 4, faces, k) == -1)
            return 0;  /* Face could not be oriented */

    return 1;
}

static int priority_vertices_ignoring_initial_tetra(const s_points *points, const bool isused[points->N], double EPS, int out_indices[points->N-4])
{   /* Return 0 if error, 1 if OK */
    /* Coordinates of the centre of the remaining point set */
    int N_aux = points->N - 4;
    s_point meanp = {0};
    int aux_count = 0;
    for (int ii=0; ii<points->N; ii++) {
        if (!isused[ii]) {
            meanp = sum_points(meanp, points->p[ii]);
            aux_count++;
        }
    }
    assert(N_aux == aux_count  && "There are more points used than 4.");
    meanp = scale_point(meanp, 1.0/N_aux);

    /* Relative distance of points from the center */
    s_point span = span_points(points);  /* Used for normalizing, considering ALL input points */
    if (fabs(span.x) < EPS || fabs(span.y) < EPS || fabs(span.z) < EPS) return 0;  /* Points do not span the 3 dimensions. */

    /* This function is only called once, so mallocing in here and freeing is no problem */
    s_indexed_double *reldist2 = malloc(N_aux * sizeof(s_indexed_double));
    if (!reldist2) return 0;
    for (int ii=0, jj=0; ii<points->N; ii++) {
        if (!isused[ii]) {
            s_point scaled = {{{ (points->p[ii].x-meanp.x)/span.x, 
                                 (points->p[ii].y-meanp.y)/span.y, 
                                 (points->p[ii].z-meanp.z)/span.z }}};
            reldist2[jj].val = norm_squared(scaled);
            reldist2[jj].idx = ii;  /* Store original index */
            jj++;
        }
    }
    
    /* Sort by relative distance ASCENDING (we will read from the end) */
    qsort(reldist2, N_aux, sizeof(s_indexed_double), &cmp_indexed_double_asc);

    /* Rescue original indexing */
    for (int i=0; i<N_aux; i++) out_indices[i] = reldist2[i].idx;

    free(reldist2);
    return 1;
}


/* Visible and non visible faces, horizon */
static int find_plane(const s_points *points, int Nfaces, int faces[Nfaces*3], s_list *cop_fids, double EPS_degenerate, s_point *out_n) 
{   /* Iterate through coplanar faces until one is not degenerate and extract n, u, v */
    for (int ii=0; ii<(int)cop_fids->N; ii++) {   
        int fid; list_get_value(cop_fids, ii, &fid);
        int f1=faces[fid*3+0], f2=faces[fid*3+1], f3=faces[fid*3+2];
        s_point plane[3] = { points->p[f1], points->p[f2], points->p[f3] };
        if (basis_vectors_plane(plane, EPS_degenerate, out_n, NULL, NULL)) return 1;
    }
    return 0;
}

static void drop_to_2D(const s_point p, int coord_to_drop, double out[2])
{
    int i1 = (coord_to_drop + 1) % 3;
    int i2 = (coord_to_drop + 2) % 3;
    out[0] = p.coords[i1];
    out[1] = p.coords[i2];
}

static int coplanar_visible(const s_points *points, s_point p, int Nfaces, int faces[Nfaces*3], s_list *cop_fids, double EPS_degenerate, s_list *out_indicator)
{   /* Called by mark_visible_faces_from_point, treats coplanar case. */
    int N_visible = 0;

    /* Project to 2D */
    s_point n;
    if (!find_plane(points, Nfaces, faces, cop_fids, EPS_degenerate, &n)) return 0;
    int coord_to_drop = coord_with_largest_component_3D(n);
    double p2D[2]; drop_to_2D(p, coord_to_drop, p2D);
    
    /* Count how many times each edge appears (canonicalized) */
    s_indexed_edge edges[cop_fids->N * 3];  // TEMPORARY TODO
    for (int ii=0; ii<(int)cop_fids->N; ii++) {  /* Store all edges */
        int fid; list_get_value(cop_fids, ii, &fid);
        int vA = faces[fid*3+0], vB = faces[fid*3+1];
        if (vA > vB) { int t = vA; vA = vB; vB = t; }
        edges[ii*3+0].e[0] = vA; edges[ii*3+0].e[1] = vB; 
        edges[ii*3+0].opp = faces[fid*3+2]; edges[ii*3+0].fid = fid;

        vA = faces[fid*3+1], vB = faces[fid*3+2];
        if (vA > vB) { int t = vA; vA = vB; vB = t; }
        edges[ii*3+1].e[0] = vA; edges[ii*3+1].e[1] = vB; 
        edges[ii*3+1].opp = faces[fid*3+0]; edges[ii*3+1].fid = fid;

        vA = faces[fid*3+0], vB = faces[fid*3+2];
        if (vA > vB) { int t = vA; vA = vB; vB = t; }
        edges[ii*3+2].e[0] = vA; edges[ii*3+2].e[1] = vB; 
        edges[ii*3+2].opp = faces[fid*3+1]; edges[ii*3+2].fid = fid;
    }
    /* sort the edges */
    qsort(edges, cop_fids->N*3, sizeof(s_indexed_edge), cmp_indexed_edge);

    /* Determine visibility of boundary (exterior, count = 1) edges */
    int e = 0;
    while (e < (int)cop_fids->N*3) {
        /* find run of identical edges */
        int run_start = e;
        int run_end = e + 1;
        while (run_end < (int)cop_fids->N*3 &&
               edges[run_end].e[0] == edges[run_start].e[0] &&
               edges[run_end].e[1] == edges[run_start].e[1]) {
            run_end++;
        }
        int run_len = run_end - run_start;

        if (run_len == 1) {  /* Edge is unique -> exterior (equivalence class of edges!) */
            double pA[2]; drop_to_2D(points->p[edges[e].e[0]], coord_to_drop, pA);
            double pB[2]; drop_to_2D(points->p[edges[e].e[1]], coord_to_drop, pB);
            double pC[2]; drop_to_2D(points->p[edges[e].opp], coord_to_drop, pC);
            if (orient2d(pA, pB, pC) < 0) {  /* Wrong order */
                double t[2] = {pA[0], pA[1]};
                pA[0] = pB[0]; pA[1] = pB[1];
                pB[0] = t[0]; pB[1] = t[1];
            }

            if (orient2d(pA, pB, p2D) < 0) {  /* Face is visible! */
                int one = 1;
                list_change_entry(out_indicator, edges[e].fid, &one);
                N_visible++;
            }
        }
        e = run_end; 
    }
    return N_visible;
}


static int visible_faces_from_point(const s_points *points, int Nfaces, int faces[Nfaces*3], s_point p, double EPS, s_list *coplanar_fids, s_list *out_indicator)
{   /* Returns -1 if error, Nvisible >= 0 if OK */
    /* A face is visible if its normal points to the halfspace containing the point,
       or if it is coplanar and the point lies strictly outside the triangle */
    if (!list_ensure_capacity(out_indicator, Nfaces))
        { fprintf(stderr, "visible_faces_from_point: error list.\n"); return -1; }
    list_memset0(out_indicator);
    int N_visible = 0;
    coplanar_fids->N = 0; 

    for (int j = 0; j < Nfaces; ++j) {
        s_point face_pts[3];
        vertices_face(points, faces, j, face_pts);
        int o = orientation_robust(face_pts, p);
        if (o < 0) {  /* Point is visible, it lies on the side pointed to by face normal  (above the plane) */
            int one = 1;
            list_change_entry(out_indicator, j, &one);
            N_visible++;
        } else if (o == 0) {  /* Point is coplanar. If inside triangle, return. */
            e_geom_test test = test_point_in_triangle_3D(face_pts, p, EPS, 0);
            if (test == TEST_BOUNDARY || test == TEST_IN) return 0;
            if (!list_push(coplanar_fids, &j)) return 0;
        }
    }

    N_visible += coplanar_visible(points, p, Nfaces, faces, coplanar_fids, EPS, out_indicator);

    return N_visible;
}


static int nonvisible_faces_sharing_edge_with_face(int Nfaces, int faces[Nfaces], int is_visible[Nfaces], int query_faceid, s_list *out_indicator)
{   /* Returns 0 if error, 1 if OK */
    if (!list_ensure_capacity(out_indicator, Nfaces))
        { fprintf(stderr, "nonvisible_faces_sharing_edge_with_face: error list.\n"); return 0; }
    list_memset0(out_indicator);

    for (int i=0; i<Nfaces; i++) {
        if (is_visible[i]) continue;
        
        /* Count shared vertices */
        int N_shared_vertices = 0;
        for (int a=0; a<3; a++)
            for (int b=0; b<3; b++)
                if (faces[query_faceid*3+a] == faces[i*3+b])
                    N_shared_vertices++;

        if (N_shared_vertices == 2) {
            int one = 1;
            list_change_entry(out_indicator, i, &one);
        }
        assert(N_shared_vertices < 3);  /* Non-manifold? */
    }
    return 1;
}

static int extract_horizon(int Nfaces, int faces[Nfaces], int is_visible[Nfaces], s_list *buff_edges, s_list *nvf_sharing_edge, int *out_Nhorizon, s_list *horizon)         
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
            int *ind = list_get_ptr(nvf_sharing_edge, j);
            if (!*ind) continue;

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
        int *prev = list_get_ptr(edges, 2*e-2);
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
static void delete_visible_faces(int Nfaces, int faces[Nfaces*3], int faces_isvisible[Nfaces])
{
    for (int j=0, l=0; j<Nfaces; j++) if (!faces_isvisible[j]) {  
        // i.e. keep those which are non visible
        faces[l*3+0] = faces[j*3+0];
        faces[l*3+1] = faces[j*3+1];
        faces[l*3+2] = faces[j*3+2];
        l++;
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
        Nfaces++;
    }
    
    /* Orient each new face properly */
    for (int k=start; k<Nfaces; k++)
        orient_face_if_needed(points, isused, Nfaces, faces->items, k);
    
    return 1;
}


/* Main algorithm */
int quickhull_3d(const s_points *in_vertices, double EPS_degenerate, int *out_Nused, bool buff_isused[in_vertices->N], int **out_faces, int *N_out_faces) 
{   /* Returns:
       -2 if error initializing tetrahedron. In_vertices degenerate or faces too small?
       -1 if error (memory, reached max_faces, ...)
       1 if output hull is exact (All faces are big enough)
       0 if output hull is non-exact (ignored any face with too small area) 
    */
    if (!in_vertices || !buff_isused || !out_faces || !N_out_faces) return -1;
    *out_Nused = 0; *out_faces = NULL;  *N_out_faces = 0;
    s_list faces = {0}, faces_isvisible = {0}, horizon = {0}, AUX_nvf_sharing_edge = {0}, AUX_edges = {0}, AUX_cop_fids = {0};
    int *pleft = NULL;

    if(in_vertices->N <= 3 ) return -2;


    /* The initial convex hull is a tetrahedron with 4 faces (simplex) */
    int Nfaces = 4;
    faces = list_initialize(sizeof(int), 12);
    if (!faces.items) goto error;
    if (!initial_tetrahedron(in_vertices, buff_isused, faces.items)) goto error_init;
    if (in_vertices->N == 4) {
        *out_Nused = 4;
        *out_faces = faces.items;
        *N_out_faces = 4;
        return 1;
    }
    

    /* Initialize the vector of points left. The points with the larger relative
     distance from the center are scanned first, which are in the END. */
    *out_Nused = 4;
    int N_pleft = in_vertices->N - 4;
    pleft = malloc(N_pleft * sizeof(int));
    if (!pleft) goto error;
    if (!priority_vertices_ignoring_initial_tetra(in_vertices, buff_isused, EPS_degenerate, pleft)) goto error_init;


    /* The main loop for the quickhull algorithm */
    faces_isvisible = list_initialize(sizeof(int), 0);
    horizon = list_initialize(sizeof(int), 0);
    AUX_nvf_sharing_edge = list_initialize(sizeof(int), 0); 
    AUX_edges = list_initialize(sizeof(int), 0);  
    AUX_cop_fids = list_initialize(sizeof(int), 0);  
    if (!faces_isvisible.items || !horizon.items || !AUX_nvf_sharing_edge.items) goto error;
    int out_is_exact = 1;
    while (N_pleft > 0) {
        /* Process the LAST element of pleft */
        N_pleft--;
        int current_id = pleft[N_pleft];
        s_point current_p = in_vertices->p[current_id];

        /* Mark visible faces from this point */
        int N_vf = visible_faces_from_point(in_vertices, Nfaces, faces.items, current_p, EPS_degenerate, &AUX_cop_fids, &faces_isvisible);
        if (N_vf == -1) goto error;

        /* Proceed if N_visible_faces > 0 */
        // assert(N_vf != Nfaces && "Point sees all faces?");
        if (N_vf == Nfaces) { fprintf(stderr, "quickhull_3d: Point sees all faces.\n"); goto error; }
        if (N_vf == 0) continue;  

        /* Create horizon */
        int Nhorizon;
        if (!extract_horizon(Nfaces, faces.items, faces_isvisible.items, &AUX_edges, &AUX_nvf_sharing_edge, &Nhorizon, &horizon))
            goto error;

        buff_isused[current_id] = true;
        (*out_Nused)++;

        /* Delete visible faces */
        delete_visible_faces(Nfaces, faces.items, faces_isvisible.items);
        Nfaces -= N_vf;
        
        /* Add faces connecting horizon to the new point */
        if(!add_faces_from_horizon(in_vertices, buff_isused, Nfaces, &faces, Nhorizon, horizon.items, current_id))
            goto error;
        Nfaces += Nhorizon;
    }
    

    /* clean-up and exit */
    free(pleft);
    free_list(&faces_isvisible); 
    free_list(&horizon);
    free_list(&AUX_nvf_sharing_edge);
    *out_faces = realloc(faces.items, Nfaces * 3 * sizeof(int));
    *N_out_faces = Nfaces;
    return out_is_exact;
    
    error_init:
        fprintf(stderr, "Error in 'quickhull_3d'. Could not setup initial tetrahedron.\n");
        if (faces.items) free_list(&faces);
        if (pleft) free(pleft);
        return -2;

    error:
        fprintf(stderr, "Error in 'quickhull_3d'. Faces N = %d / %d\n", Nfaces, CH_MAX_NUM_FACES);
        if (faces.items) free_list(&faces);
        if (pleft) free(pleft);
        if (faces_isvisible.items) free_list(&faces_isvisible); 
        if (horizon.items) free_list(&horizon);
        if (AUX_nvf_sharing_edge.items) free_list(&AUX_nvf_sharing_edge);
        *out_faces = NULL;
        *N_out_faces = 0;
        return -1;
}

