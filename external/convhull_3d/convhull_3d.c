// TODO: Change name of dbl_w_idx
// Better mallocing of arrays...
// Deal with degenerate faces... simply exit?
// Sometimes the output results in degenerate vectors: norm(face_normal) < 1e-14

#include "../../include/geometry.h"  // TODO temporal

#include "convhull_3d.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <errno.h> 
#include <assert.h>


#define MIN(a,b) (( (a) < (b) ) ? (a) : (b) )
#define MAX(a,b) (( (a) > (b) ) ? (a) : (b) )

#define CH_MAX_NUM_FACES 50000

/* structs for qsort */
typedef struct dbl_w_idx {
    double val;
    int idx;
} dbl_w_idx;

/* internal functions */
static void* default_memory_resize(void* ptr, size_t size);
static int cmp_dbl_w_idx_desc(const void *pa, const void *pb);
// static void ismember(int*, int*, int*, int, int);
// static void vertices_face(const s_points *points, const int *faces, int face_id, s_point out[3]);

static void* default_memory_resize(void* ptr, size_t size)
{
    if (!ptr) return(malloc(size));
    return realloc(ptr, size);
}

static int cmp_dbl_w_idx_desc(const void *pa, const void *pb) {
	const dbl_w_idx *a = (const dbl_w_idx*)pa;
	const dbl_w_idx *b = (const dbl_w_idx*)pb;
	if (a->val < b->val) return  1;  /* b before a => descending */
	if (a->val > b->val) return -1;
	return 0;
}

// static void ismember
// (
//     int* pLeft,          /* left vector; nLeftElements x 1 */
//     int* pRight,         /* right vector; nRightElements x 1 */
//     int* pOut,           /* 0, unless pRight elements are present in pLeft then 1; nLeftElements x 1 */
//     int nLeftElements,   /* number of elements in pLeft */
//     int nRightElements   /* number of elements in pRight */
// )
// {
//     int i, j;
//     memset(pOut, 0, nLeftElements*sizeof(int));
//     for(i=0; i< nLeftElements; i++)
//         for(j=0; j< nRightElements; j++)
//             if(pLeft[i] == pRight[j] )
//                 pOut[i] = 1;
// }




static void vertices_face(const s_points *points, const int *faces, int face_id, s_point out[3])
{
    int i;
    i = faces[face_id*3+0];  assert(i < points->N && i >= 0);
    out[0] = points->p[i];

    i = faces[face_id*3+1];  assert(i < points->N && i >= 0);
    out[1] = points->p[i];

    i = faces[face_id*3+2];  assert(i < points->N && i >= 0);
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

    for (int ii = start; ii < Np; ++ii) {
        if (isused[ii]) {
            if (ii != face_ids[0] && ii != face_ids[1] && ii != face_ids[2])
                return ii;
        }
    }
    return -1;
}


static int orient_face_if_needed(const s_points *points, int isused[points->N], int Nfaces, int faces[3*Nfaces], int face_id)
{   /* Returns 1 if reoriented, 0 if not, -1 if error (Could not orient it!) */
    int face_vids[3] = {faces[face_id*3+0], faces[face_id*3+1], faces[face_id*3+2]};
    s_point face_vertices[3];
    vertices_face(points, faces, face_id, face_vertices);
    
    /* Find noncoplanar point which is part of the hull and not in the current new face */
    int p = -1;
    p = next_vid_isused_notinface(points->N, isused, face_vids, p);
    if (p == -1) return -1;

    int o = orientation(face_vertices, points->p[p]);
    while (o == 0) {
        p = next_vid_isused_notinface(points->N, isused, face_vids, p);
        if (p == -1) return -1;  
        o = orientation(face_vertices, points->p[p]);
    }

    /* Orient faces so that each point on the original simplex can't see the opposite face */
    if (o < 0) {
        flip_face(faces, face_id);
        return 1;
        /* Check DEBUG (TODO remove assert? I think unnecessary)*/
        vertices_face(points, faces, face_id, face_vertices);
        assert(orientation(face_vertices, points->p[p]) > 0);
    }

    return 0;
}


static int find_any_non_coplanar_quad(const s_points *points, int out[4])
{   /* brute-force combinations i<j<k<l */
	int N = points->N;
    assert(N >= 4  && "Not enough points");
	for (int i=0; i<N-3; i++) {
		for (int j=i+1; j<N-2; j++) {
			for (int k=j+1; k<N-1; k++) {
				s_point tri[3] = {points->p[i], points->p[j], points->p[k]};
				for (int l=k+1; l<N; l++) {
					int o = orientation(tri, points->p[l]);
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


static void priority_vertices_ignoring_initial_tetra(const s_points *points, 
                                                     const int isused[points->N],
                                                     int out_indices[points->N-4])
{
    /* Coordinates of the center of the remaining point set */
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
    double EPS = 1e-12;
    assert(fabs(span.x) > EPS && 
           fabs(span.y) > EPS &&
           fabs(span.z) > EPS && "Points do not span the 3 dimensions. TODO");

    dbl_w_idx *reldist2 = malloc(N_aux * sizeof(dbl_w_idx));
    for (int ii=0, jj=0; ii<points->N; ii++) {
        if (!isused[ii]) {
            s_point scaled = {{{(points->p[ii].x-meanp.x)/span.x, 
                                (points->p[ii].y-meanp.y)/span.y, 
                                (points->p[ii].z-meanp.z)/span.z}}};
            reldist2[jj].val = norm_squared(scaled);
            reldist2[jj].idx = ii;  /* Store original index */
            jj++;
        }
    }
    
    /* Sort by relative distance descending */
    qsort(reldist2, N_aux, sizeof(dbl_w_idx), &cmp_dbl_w_idx_desc);

    /* Rescue original indexing */
    for (int i=0; i<N_aux; i++) out_indices[i] = reldist2[i].idx;

    free(reldist2);
}


static int visible_faces_from_point(const s_points *points, int Nfaces, int faces[Nfaces*3], s_point p, int out_indicator[Nfaces])
{   /* A face is visible if its normal points to the halfspace containing the point,
        or if it is coplanar and the point lies strictly outside the triangle */
    memset(out_indicator, 0, Nfaces*sizeof(int));
    int N_visible = 0;
    for (int j = 0; j < Nfaces; ++j) {
        s_point face_pts[3];
        vertices_face(points, faces, j, face_pts);
        int o = orientation(face_pts, p);
        if (o < 0) {  /* Point is visible, it lies on the side pointed to by face normal  (above the plane) */
            out_indicator[j] = 1;
            N_visible++;
        } else if (o == 0) {  /* Point is coplanar. Outside of inside triangle? */
            if (in_triangle_3d(face_pts, p) == 0) {  /* Strictly outside */
                out_indicator[j] = 1;
                N_visible++;
            }
        }
    }
    return N_visible;
}


static void nonvisible_faces_sharing_edge_with_face(int Nfaces, int faces[Nfaces], int is_visible[Nfaces], int query_faceid, int out_indicator[Nfaces])
{ 
    memset(out_indicator, 0, Nfaces*sizeof(int));

    for (int i=0; i<Nfaces; i++) {
        if (is_visible[i]) continue;
        
        /* Count shared vertices */
        int N_shared_vertices = 0;
        for (int a=0; a<3; a++)
            for (int b=0; b<3; b++)
                if (faces[query_faceid*3+a] == faces[i*3+b])
                    N_shared_vertices++;

        if (N_shared_vertices == 2) out_indicator[i] = 1;
        assert(N_shared_vertices < 3);  /* Non-manifold? */
    }
}


static void extract_horizon(int Nfaces, int faces[Nfaces], int is_visible[Nfaces], int *out_Nhorizon, int **out_horizon)         
{   /* size of horizon: 2*out_Nhorizon */
    free(*out_horizon);

    int Nhorizon = 0;
    int *horizon = NULL;
    int *nvf_sharing_edge = malloc(Nfaces * sizeof(int));  // Indicates if non_visible vertex is shared with the current visible face
    
    for (int i=0; i<Nfaces; i++) {
        if (!is_visible[i]) continue;
        /* Mark faces that are nonvisible and share an edge with current face */
        nonvisible_faces_sharing_edge_with_face(Nfaces, faces, is_visible, i, nvf_sharing_edge);
        for (int j=0; j<Nfaces; j++) {
            if (!nvf_sharing_edge[j]) continue;

            /* This means that non visible face j shares edge with visible face i */
            /* But which edge ? */
            horizon = default_memory_resize(horizon, (Nhorizon+1)*2*sizeof(int));
            int h = 0;
            for (int a=0; a<3; a++)
                for (int b=0; b<3; b++)
                    if (faces[i*3+a] == faces[j*3+b])
                        horizon[Nhorizon*2+h++] = faces[j*3+b];
            assert(h == 2);
            Nhorizon++;
        }
    }

    *out_Nhorizon = Nhorizon;
    *out_horizon = horizon;
    free(nvf_sharing_edge);
}


void quickhull_3d(const s_points *in_vertices, int **out_faces, int *N_out_faces) 
{
    *out_faces = NULL;
    *N_out_faces = 0;

    if(in_vertices == NULL || in_vertices->N <= 3 ) return;

    /* The initial convex hull is a tetrahedron with 4 faces (simplex) */
    int Nfaces = 4;
    int *faces = malloc(4 * 3 * sizeof(int));
    int *isused = malloc(in_vertices->N * sizeof(int));
    if (!initial_tetrahedron(in_vertices, isused, faces)) {
        free(faces);
        free(isused);
        printf("DEBUG: could not construct initial tetrahedron\n");
        return;
    }

    if (in_vertices->N == 4) {
        *out_faces = faces;
        *N_out_faces = 4;
        free(isused);
        return;
    }
    

    /* Initialize the vector of points left. The points with the larger relative
     distance from the center are scanned first. */
    int N_pleft = in_vertices->N - 4;
    int *pleft = malloc(N_pleft * sizeof(int));
    priority_vertices_ignoring_initial_tetra(in_vertices, isused, pleft);


    /* The main loop for the quickhull algorithm */
    /* Use no mallocs inside, only reallocs. (TODO!) */
    int FUCKED = 0;
    int *faces_isvisible = NULL, *vf_ids = NULL, *nvf = NULL, *horizon = NULL;
    while (N_pleft > 0) {
        /* Process the first element of pleft */
        int current_id = pleft[0];
        s_point current_p = in_vertices->p[current_id];

        /* Delete the point selected */
        N_pleft--;
        for(int j=0; j<N_pleft; j++) pleft[j] = pleft[j+1];
        if (N_pleft == 0) free(pleft);
        else pleft = realloc(pleft, N_pleft * sizeof(int));
        

        /* Mark visible faces from this point */
        faces_isvisible = realloc(faces_isvisible, Nfaces * sizeof(int));
        int N_vf = visible_faces_from_point(in_vertices, Nfaces, faces, current_p, faces_isvisible);

        /* Proceed if N_visible_faces > 0 */
        if (N_vf == 0) continue;  
        isused[current_id] = 1;

        /* Create horizon */
        int Nhorizon;
        extract_horizon(Nfaces, faces, faces_isvisible, &Nhorizon, &horizon);


        /* Check if any new face would be degenerate */
        // int N_newfaces = Nhorizon;
        // int degenerate = 0;
        // for (int j=0; j<N_newfaces; j++) {
        //     int v0 = horizon[j*2+0];  int v1 = horizon[j*2+1];  int v2 = current_id;
        //     int degenerate_1 = (orient2d(in_vertices->p[v0].coords, in_vertices->p[v1].coords, in_vertices->p[v2].coords) == 0);
        //     int degenerate_2 = (fabs(norm(cross_prod(subtract_points(in_vertices->p[v1], in_vertices->p[v0]),
        //                                              subtract_points(in_vertices->p[v2], in_vertices->p[v0])))) < 1e-14);
        //     if (degenerate_1 || degenerate_2) {
        //         degenerate = 1;
        //         break;
        //     }
        // }
        // (void)degenerate;
        // if (degenerate) continue;  /* Simply ignore the point */


        /* Delete visible faces */
        if (N_vf == Nfaces) {
            fprintf(stderr, "WARNING: point %d sees ALL faces! N_pleft=%d\n", current_id, N_pleft);
            exit(1);
        }
        for (int j=0, l=0; j<Nfaces; j++){
            if (!faces_isvisible[j]) {  // i.e. keep those which are non visible
                faces[l*3+0] = faces[j*3+0];
                faces[l*3+1] = faces[j*3+1];
                faces[l*3+2] = faces[j*3+2];
                l++;
            }
        }
        Nfaces -= N_vf;
        

        /* Add faces connecting horizon to the new point */
        int N_newfaces = Nhorizon;
        int start = Nfaces;  /* start is the first row of the new faces */
        int N_realloc_faces = Nfaces + N_newfaces;
        if (N_realloc_faces > CH_MAX_NUM_FACES) N_realloc_faces = CH_MAX_NUM_FACES+1;

        faces = realloc(faces, N_realloc_faces * 3 * sizeof(int));
        for (int j=0; j<N_newfaces; j++) {
            faces[Nfaces*3+0] = horizon[j*2+0];
            faces[Nfaces*3+1] = horizon[j*2+1];
            faces[Nfaces*3+2] = current_id;
            Nfaces++;

            if (Nfaces > CH_MAX_NUM_FACES) {
                FUCKED = 1;
                Nfaces = 0;
                break;
            }
        }
        
        /* Orient each new face properly */
        for (int k=start; k<Nfaces; k++) {
            orient_face_if_needed(in_vertices, isused, Nfaces, faces, k);
        }
   
        if(FUCKED){
            break;
        }
    }
    
    /* output */
    if (FUCKED) {
        printf("DEBUG convhull_3d_build FUCKED!\n");
    } else {
        printf("DEBUG convhull_3d_build OK, nfaces=%d\n", Nfaces);
        *out_faces = malloc(Nfaces * 3 * sizeof(int));
        memcpy(*out_faces, faces, Nfaces * 3 * sizeof(int));
        *N_out_faces = Nfaces;
    }
    
    /* clean-up */
    free(faces);
    free(faces_isvisible); 
    free(vf_ids);
    free(nvf);
    free(horizon);
    free(isused);
}





// void convhull_3d_export_obj(
//     ch_vertex* const vertices,
//     const int nVert,
//     int* const faces,
//     const int Nfaces,
//     const int keepOnlyUsedVerticesFLAG,
//     char* const obj_filename)
// {
//     int i, j;
//     char path[256] = "\0";
//     strncpy(path, obj_filename, strlen(obj_filename));
//     FILE* obj_file;
//
//     errno = 0;
//     obj_file = fopen(strcat(path, ".obj"), "wt");
//
//     if (obj_file==NULL) {
//         printf("Error %d \n", errno);
//         printf("It's null");
//     }
//     fprintf(obj_file, "o\n");
//     double scale;
//     s_point v1, v2, normal;
//
//     /* export vertices */
//     if(keepOnlyUsedVerticesFLAG){
//         for (i = 0; i < Nfaces; i++)
//             for(j=0; j<3; j++)
//                 fprintf(obj_file, "v %f %f %f\n", vertices[faces[i*3+j]].x,
//                         vertices[faces[i*3+j]].y, vertices[faces[i*3+j]].z);
//     }
//     else {
//         for (i = 0; i < nVert; i++)
//             fprintf(obj_file, "v %f %f %f\n", vertices[i].x,
//                     vertices[i].y, vertices[i].z);
//     }
//     
//     /* export the face normals */
//     for (i = 0; i < Nfaces; i++){
//         /* calculate cross product between v1-v0 and v2-v0 */
//         v1 = vertices[faces[i*3+1]];
//         v2 = vertices[faces[i*3+2]];
//         v1.x -= vertices[faces[i*3]].x;
//         v1.y -= vertices[faces[i*3]].y;
//         v1.z -= vertices[faces[i*3]].z;
//         v2.x -= vertices[faces[i*3]].x;
//         v2.y -= vertices[faces[i*3]].y;
//         v2.z -= vertices[faces[i*3]].z;
//         normal = cross_prod(&v1, &v2);
//         
//         /* normalise to unit length */
//         scale = ((double)1.0)/(sqrt(pow(normal.x, (double)2.0)+pow(normal.y, (double)2.0)+pow(normal.z, (double)2.0))+(double)2.23e-9);
//         normal.x *= scale;
//         normal.y *= scale;
//         normal.z *= scale;
//         fprintf(obj_file, "vn %f %f %f\n", normal.x, normal.y, normal.z);
//     }
//     
//     /* export the face indices */
//     if(keepOnlyUsedVerticesFLAG){
//         for (i = 0; i < Nfaces; i++){
//             /* vertices are in same order as the faces, and normals are in order */
//             fprintf(obj_file, "f %u//%u %u//%u %u//%u\n",
//                     i*3 + 1, i + 1,
//                     i*3+1 + 1, i + 1,
//                     i*3+2 + 1, i + 1);
//         }
//     }
//     else {
//         /* just normals are in order  */
//         for (i = 0; i < Nfaces; i++){
//             fprintf(obj_file, "f %u//%u %u//%u %u//%u\n",
//                     faces[i*3] + 1, i + 1,
//                     faces[i*3+1] + 1, i + 1,
//                     faces[i*3+2] + 1, i + 1);
//         }
//     }
//     fclose(obj_file);
// }
//
//
// void extract_vertices_from_obj_file (
//     char* const obj_filename,
//     ch_vertex** out_vertices,
//     int* out_nVert)
// {
//     extract_vertices_from_obj_file_alloc(obj_filename, out_vertices, out_nVert, NULL);
// }
//
//
// void extract_vertices_from_obj_file_alloc
// (
//     char* const obj_filename,
//     ch_vertex** out_vertices,
//     int* out_nVert,
//     void* allocator
// )
// {
//     FILE* obj_file;
//     obj_file = fopen(strcat(obj_filename, ".obj"), "r");
//     
//     /* determine number of vertices */
//     unsigned int nVert = 0;
//     char line[256];
//     while (fgets(line, sizeof(line), obj_file)) {
//         char* vexists = strstr(line, "v ");
//         if(vexists!=NULL)
//             nVert++;
//     }
//     (*out_nVert) = nVert;
//     (*out_vertices) = (ch_vertex*)ch_stateful_malloc(allocator, nVert*sizeof(ch_vertex));
//     
//     /* extract the vertices */
//     rewind(obj_file);
//     int i=0;
//     int vertID, prev_char_isDigit, current_char_isDigit;
//     char vert_char[256] = { 0 }; 
//     while (fgets(line, sizeof(line), obj_file)) {
//         char* vexists = strstr(line, "v ");
//         if(vexists!=NULL){
//             prev_char_isDigit = 0;
//             vertID = -1;
//             for(size_t j=0; j<strlen(line)-1; j++){
//                 if(isdigit(line[j])||line[j]=='.'||line[j]=='-'||line[j]=='+'||line[j]=='E'||line[j]=='e'){
//                     vert_char[strlen(vert_char)] = line[j];
//                     current_char_isDigit = 1;
//                 }
//                 else
//                     current_char_isDigit = 0;
//                 if((prev_char_isDigit && !current_char_isDigit) || j ==strlen(line)-2 ){
//                     vertID++;
//                     if(vertID>4){
//                         /* not a valid file */
//                         ch_stateful_free(allocator, (*out_vertices));
//                         (*out_vertices) = NULL;
//                         (*out_nVert) = 0;
//                         return;
//                     }
//                     (*out_vertices)[i].v[vertID] = (double)atof(vert_char);
//                     memset(vert_char, 0, 256 * sizeof(char));
//                 }
//                 prev_char_isDigit = current_char_isDigit;
//             }
//             i++;
//         }
//     }
// }


