#include "points.h"  
#include "gtests.h"
#include "ch_quickhull3D.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <errno.h> 
#include <assert.h>


/* Sorting by distance and remembering original index */
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


/* Dealing with memory efficiently */
typedef struct int_list {
    int *list;
    int Nmax;
} s_int_list;

static s_int_list initialize_int_list(int Nmax) {
    if (Nmax <= 0) Nmax = CH_N_INIT_INT_LIST;
    s_int_list out = { .Nmax = Nmax, 
                       .list = malloc(Nmax * sizeof(int)) };
    return out;
}

static int increase_memory_int_list_if_needed(s_int_list *int_list, int N_needed) {
    while (N_needed >= int_list->Nmax) {
        int *tmp = realloc(int_list->list, 2 * int_list->Nmax * sizeof(int));
        if (!tmp) {
            return 0;
        }
        int_list->list = tmp;
        int_list->Nmax *= 2;
    }
    return 1;
}

static void free_int_list(s_int_list *int_list) {
    free(int_list->list);
    memset(int_list, 0, sizeof(s_int_list));
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

static int next_vid_isused_notinface(int Np, int isused[Np], int face_ids[3], int current_id)
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

static int orient_face_if_needed(const s_points *points, int isused[points->N], int Nfaces, int faces[3*Nfaces], int face_id)
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
	for (int i=0; i<N-3; i++) {
		for (int j=i+1; j<N-2; j++) {
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
		}
	}
	return 0;  /* all points are degenerate (coplanar/collinear/coincident) */
}

static int initial_tetrahedron(const s_points *points, int isused[points->N], int faces[4*3])
{
    /* Set the indices of the points defining the face  */
    int vertex_ids[4];
    if (!find_any_non_coplanar_quad(points, vertex_ids)) return 0;

    memset(isused, 0, points->N * sizeof(int));
    isused[vertex_ids[0]] = 1;  isused[vertex_ids[1]] = 1;
    isused[vertex_ids[2]] = 1;  isused[vertex_ids[3]] = 1;

    for (int i=0; i<4; i++) {
        for (int j=0, k=0; j<4; j++){
            if (vertex_ids[j] != vertex_ids[i]) faces[i*3+k++] = vertex_ids[j];
        }
    }

    /* Ensure that faces are correctly oriented */
    for (int k=0; k<4; k++) {
        if (orient_face_if_needed(points, isused, 4, faces, k) == -1) {
            return 0;  /* Face could not be oriented */
        }
    }
    return 1;
}

static int priority_vertices_ignoring_initial_tetra(const s_points *points, const int isused[points->N], double EPS, int out_indices[points->N-4])
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
static int visible_faces_from_point(const s_points *points, int Nfaces, int faces[Nfaces*3], s_point p, double EPS, s_int_list *out_indicator)
{   /* Returns -1 if error, Nvisible > 0 if OK */
    /* A face is visible if its normal points to the halfspace containing the point,
       or if it is coplanar and the point lies strictly outside the triangle */
    if (!increase_memory_int_list_if_needed(out_indicator, Nfaces)) return -1;
    memset(out_indicator->list, 0, Nfaces*sizeof(int));
    int N_visible = 0;
    for (int j = 0; j < Nfaces; ++j) {
        s_point face_pts[3];
        vertices_face(points, faces, j, face_pts);
        int o = orientation_robust(face_pts, p);
        if (o < 0) {  /* Point is visible, it lies on the side pointed to by face normal  (above the plane) */
            out_indicator->list[j] = 1;
            N_visible++;
        } else if (o == 0) {  /* Point is coplanar. Outside or inside triangle? */
            if (test_point_in_triangle_3D(face_pts, p, EPS, 0) == TEST_OUT) {  /* Strictly outside */
                out_indicator->list[j] = 1;
                N_visible++;
            }
        }
    }
    return N_visible;
}

static int nonvisible_faces_sharing_edge_with_face(int Nfaces, int faces[Nfaces], int is_visible[Nfaces], int query_faceid, s_int_list *out_indicator)
{   /* Returns 0 if error, 1 if OK */
    if (!increase_memory_int_list_if_needed(out_indicator, Nfaces)) return 0;
    memset(out_indicator->list, 0, Nfaces*sizeof(int));

    for (int i=0; i<Nfaces; i++) {
        if (is_visible[i]) continue;
        
        /* Count shared vertices */
        int N_shared_vertices = 0;
        for (int a=0; a<3; a++)
            for (int b=0; b<3; b++)
                if (faces[query_faceid*3+a] == faces[i*3+b])
                    N_shared_vertices++;

        if (N_shared_vertices == 2) out_indicator->list[i] = 1;
        assert(N_shared_vertices < 3);  /* Non-manifold? */
    }
    return 1;
}

static int extract_horizon(int Nfaces, int faces[Nfaces], int is_visible[Nfaces], s_int_list *nvf_sharing_edge, int *out_Nhorizon, s_int_list *horizon)         
{   /* Returns 0 if eror, 1 if OK */
    /* size of horizon: 2*out_Nhorizon */
    int Nhorizon = 0;
    
    for (int i=0; i<Nfaces; i++) {
        if (!is_visible[i]) continue;
        /* Mark faces that are nonvisible and share an edge with current face */
        if (!nonvisible_faces_sharing_edge_with_face(Nfaces, faces, is_visible, i, nvf_sharing_edge)) return 0;
        for (int j=0; j<Nfaces; j++) {
            if (!nvf_sharing_edge->list[j]) continue;

            /* This means that non visible face j shares edge with visible face i */
            if (!increase_memory_int_list_if_needed(horizon, (Nhorizon+1)*2)) return 0;

            /* But which edge ? */
            int h = 0;
            for (int a=0; a<3; a++)
                for (int b=0; b<3; b++)
                    if (faces[i*3+a] == faces[j*3+b])
                        horizon->list[Nhorizon*2+h++] = faces[j*3+b];
            assert(h == 2);
            Nhorizon++;
        }
    }

    *out_Nhorizon = Nhorizon;
    return 1;
}


/* Add and delete faces */
static void delete_visible_faces(int Nfaces, int faces[Nfaces*3], int faces_isvisible[Nfaces])
{
    for (int j=0, l=0; j<Nfaces; j++){
        if (!faces_isvisible[j]) {  // i.e. keep those which are non visible
            faces[l*3+0] = faces[j*3+0];
            faces[l*3+1] = faces[j*3+1];
            faces[l*3+2] = faces[j*3+2];
            l++;
        }
    }
}

static int add_faces_from_horizon(const s_points *points, int isused[points->N], int Nfaces, s_int_list *faces, int Nhorizon, int horizon[Nhorizon*2], int query_pid)
{   /* Returns 0 if error, -1 if OK */
    int start = Nfaces;  /* start is the first row of the new faces */
    int N_realloc_faces = Nfaces + Nhorizon;
    if (N_realloc_faces >= CH_MAX_NUM_FACES) return 0;

    if (!increase_memory_int_list_if_needed(faces, N_realloc_faces * 3)) return 0;
    for (int j=0; j<Nhorizon; j++) {
        faces->list[Nfaces*3+0] = horizon[j*2+0];
        faces->list[Nfaces*3+1] = horizon[j*2+1];
        faces->list[Nfaces*3+2] = query_pid;
        Nfaces++;
    }
    
    /* Orient each new face properly */
    for (int k=start; k<Nfaces; k++)
        orient_face_if_needed(points, isused, Nfaces, faces->list, k);
    
    return 1;
}


/* Main algorithm */
int quickhull_3d(const s_points *in_vertices, double EPS_degenerate, int buff_isused[in_vertices->N], int **out_faces, int *N_out_faces) 
{   /* Returns:
       -2 if error initializing tetrahedron. In_vertices degenerate or faces too small?
       -1 if error (memory, reached max_faces, ...)
       1 if output hull is exact (All faces are big enough)
       0 if output hull is non-exact (ignored any face with too small area) 
    */
    if (!in_vertices || !buff_isused || !out_faces || !N_out_faces) return -1;
    *out_faces = NULL;  *N_out_faces = 0;
    memset(buff_isused, 0, in_vertices->N * sizeof(int));
    s_int_list faces = {0}, faces_isvisible = {0}, horizon = {0}, AUX_nvf_sharing_edge = {0};
    int *pleft = NULL;

    if(in_vertices->N <= 3 ) return -2;


    /* The initial convex hull is a tetrahedron with 4 faces (simplex) */
    int Nfaces = 4;
    faces = initialize_int_list(12);
    if (!faces.list) goto error;
    if (!initial_tetrahedron(in_vertices, buff_isused, faces.list)) goto error_init;
    if (in_vertices->N == 4) {
        *out_faces = faces.list;
        *N_out_faces = 4;
        return 1;
    }
    

    /* Initialize the vector of points left. The points with the larger relative
     distance from the center are scanned first, which are in the END. */
    int N_pleft = in_vertices->N - 4;
    pleft = malloc(N_pleft * sizeof(int));
    if (!pleft) goto error;
    if (!priority_vertices_ignoring_initial_tetra(in_vertices, buff_isused, EPS_degenerate, pleft)) goto error_init;


    /* The main loop for the quickhull algorithm */
    faces_isvisible = initialize_int_list(0);
    horizon = initialize_int_list(0);
    AUX_nvf_sharing_edge = initialize_int_list(0);  /* Used internally, but malloced outside for proper memory control */
    if (!faces_isvisible.list || !horizon.list || !AUX_nvf_sharing_edge.list) goto error;
    int out_is_exact = 1;
    while (N_pleft > 0) {
        /* Process the LAST element of pleft */
        N_pleft--;
        int current_id = pleft[N_pleft];
        s_point current_p = in_vertices->p[current_id];

        /* Mark visible faces from this point */
        int N_vf = visible_faces_from_point(in_vertices, Nfaces, faces.list, current_p, EPS_degenerate, &faces_isvisible);
        if (N_vf == -1) goto error;

        /* Proceed if N_visible_faces > 0 */
        assert(N_vf != Nfaces && "Point sees all faces?");
        if (N_vf == 0) continue;  

        /* Create horizon */
        int Nhorizon;
        if (!extract_horizon(Nfaces, faces.list, faces_isvisible.list, &AUX_nvf_sharing_edge, &Nhorizon, &horizon))
            goto error;

        buff_isused[current_id] = 1;

        /* Delete visible faces */
        delete_visible_faces(Nfaces, faces.list, faces_isvisible.list);
        Nfaces -= N_vf;
        
        /* Add faces connecting horizon to the new point */
        if(!add_faces_from_horizon(in_vertices, buff_isused, Nfaces, &faces, Nhorizon, horizon.list, current_id))
            goto error;
        Nfaces += Nhorizon;
    }
    

    /* clean-up and exit */
    free(pleft);
    free_int_list(&faces_isvisible); 
    free_int_list(&horizon);
    free_int_list(&AUX_nvf_sharing_edge);
    *out_faces = realloc(faces.list, Nfaces * 3 * sizeof(int));
    *N_out_faces = Nfaces;
    return out_is_exact;
    
    error_init:
        fprintf(stderr, "Error in 'quickhull_3d'. Could not setup initial tetrahedron.\n");
        if (faces.list) free_int_list(&faces);
        if (pleft) free(pleft);
        return -2;

    error:
        fprintf(stderr, "Error in 'quickhull_3d'. Maybe reached max faces? N = %d / %d\n", Nfaces, CH_MAX_NUM_FACES);
        if (faces.list) free_int_list(&faces);
        if (pleft) free(pleft);
        if (faces_isvisible.list) free_int_list(&faces_isvisible); 
        if (horizon.list) free_int_list(&horizon);
        if (AUX_nvf_sharing_edge.list) free_int_list(&AUX_nvf_sharing_edge);
        *out_faces = NULL;
        *N_out_faces = 0;
        return -1;
}

