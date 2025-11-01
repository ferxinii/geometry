#include "../../include/geometry.h"  // TODO temporal

/*
 Copyright (c) 2017-2021 Leo McCormack
 
 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:
 
 The above copyright notice and this permission notice shall be included in
 all copies or substantial portions of the Software.
 
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 THE SOFTWARE.
*/
/*
 * Filename:
 *     convhull_3d.h
 * Description:
 *     A header only C implementation of the 3-D quickhull algorithm.
 *     The code is largely derived from the "computational-geometry-toolbox"
 *     by George Papazafeiropoulos (c) 2014, originally distributed under
 *     the BSD (2-clause) license.
 *     To include this implementation in a project, simply add this:
 *         #define CONVHULL_3D_ENABLE
 *         #include "convhull_3d.h"
 *     By default, the algorithm uses double floating point precision. To
 *     use single precision (less accurate but quicker), also add this:
 *         #define CONVHULL_3D_USE_SINGLE_PRECISION
 *     If your project has CBLAS linked, then you can also speed things up
 *     a tad by adding this:
 *         #define CONVHULL_3D_USE_CBLAS
 *     The code is C++ compiler safe.
 *     Reference: "The Quickhull Algorithm for Convex Hull, C. Bradford
 *                 Barber, David P. Dobkin and Hannu Huhdanpaa, Geometry
 *                 Center Technical Report GCG53, July 30, 1993"
 * Dependencies:
 *     cblas (optional for speed ups, especially for very large meshes)
 *     (Available in e.g. Apple Accelerate Framework, or Intel MKL)
 * Author, date created:
 *     Leo McCormack, 02.10.2017
 */


#include "convhull_3d.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <float.h>
#include <ctype.h>
#include <string.h>
#include <errno.h> 
#include <assert.h>


#define MIN(a,b) (( (a) < (b) ) ? (a) : (b) )
#define MAX(a,b) (( (a) > (b) ) ? (a) : (b) )

#define CH_MAX_NUM_FACES 50000

/* structs for qsort */
typedef struct float_w_idx {
    double val;
    int idx;
} float_w_idx;

/* internal functions prototypes: */
static void* default_memory_resize(void* ptr, size_t size);
static int cmp_asc_float(const void*, const void*);
static int cmp_desc_float(const void*, const void*);
static void sort_float(double*, double*, int*, int, int);
static void ismember(int*, int*, int*, int, int);
static void vertices_face(const s_points *points, const int *faces, int face_id, s_point out[3]);

/* internal functions definitions: */
static void* default_memory_resize(void* ptr, size_t size)
{
    if (!ptr) return(malloc(size));
    return realloc(ptr, size);
}

static int cmp_asc_float(const void *a,const void *b) {
    struct float_w_idx *a1 = (struct float_w_idx*)a;
    struct float_w_idx *a2 = (struct float_w_idx*)b;
    if((*a1).val<(*a2).val)return -1;
    else if((*a1).val>(*a2).val)return 1;
    else return 0;
}

static int cmp_desc_float(const void *a,const void *b) {
    struct float_w_idx *a1 = (struct float_w_idx*)a;
    struct float_w_idx *a2 = (struct float_w_idx*)b;
    if((*a1).val>(*a2).val)return -1;
    else if((*a1).val<(*a2).val)return 1;
    else return 0;
}

static void sort_float(
    double* in_vec,  /* vector[len] to be sorted */
    double* out_vec, /* if NULL, then in_vec is sorted "in-place" */
    int* new_idices,   /* set to NULL if you don't need them */
    int len,           /* number of elements in vectors, must be consistent with the input data */
    int descendFLAG    /* !1:ascending, 1:descending */
)
{
    int i;
    struct float_w_idx *data;
    
    data = (float_w_idx*)malloc(len*sizeof(float_w_idx));
    for(i=0;i<len;i++) {
        data[i].val=in_vec[i];
        data[i].idx=i;
    }
    if(descendFLAG)
        qsort(data,len,sizeof(data[0]),cmp_desc_float);
    else
        qsort(data,len,sizeof(data[0]),cmp_asc_float);
    for(i=0;i<len;i++){
        if (out_vec!=NULL)
            out_vec[i] = data[i].val;
        else
            in_vec[i] = data[i].val; /* overwrite input vector */
        if(new_idices!=NULL)
            new_idices[i] = data[i].idx;
    }
    free(data);
}


static void ismember
(
    int* pLeft,          /* left vector; nLeftElements x 1 */
    int* pRight,         /* right vector; nRightElements x 1 */
    int* pOut,           /* 0, unless pRight elements are present in pLeft then 1; nLeftElements x 1 */
    int nLeftElements,   /* number of elements in pLeft */
    int nRightElements   /* number of elements in pRight */
)
{
    int i, j;
    memset(pOut, 0, nLeftElements*sizeof(int));
    for(i=0; i< nLeftElements; i++)
        for(j=0; j< nRightElements; j++)
            if(pLeft[i] == pRight[j] )
                pOut[i] = 1;
}

static void vertices_face(const s_points *points, const int *faces, int face_id, s_point out[3])
{
    int i;
    i = faces[face_id*3+0];  assert(i < points->N);
    out[0] = points->p[i];

    i = faces[face_id*3+1];  assert(i < points->N);
    out[1] = points->p[i];

    i = faces[face_id*3+2];  assert(i < points->N);
    out[2] = points->p[i];
}


/* A C version of the 3D quickhull matlab implementation from here:
 * https://www.mathworks.com/matlabcentral/fileexchange/48509-computational-geometry-toolbox?focused=3851550&tab=example
 * (*out_faces) is returned as NULL, if triangulation fails *
 * Original Copyright (c) 2014, George Papazafeiropoulos
 * Distributed under the BSD (2-clause) license
 * Reference: "The Quickhull Algorithm for Convex Hull, C. Bradford Barber, David P. Dobkin
 *             and Hannu Huhdanpaa, Geometry Center Technical Report GCG53, July 30, 1993"
 */
void convhull_3d_build(const s_points *in_vertices, int **out_faces, int *nOut_faces) 
{
    if(in_vertices->N <= 3 || in_vertices == NULL){
        (*out_faces) = NULL;
        (*nOut_faces) = 0;
        return;
    }

    /* The initial convex hull is a tetrahedron with 4 faces (simplex) */
    int nFaces = 4;
    int *faces = calloc(nFaces*3, sizeof(int));
    int *aVec = malloc(nFaces*sizeof(int));
    for (int i=0; i<4; i++) aVec[i] = i;  // {0,1,2,3}
    
    /* Store the planes of the faces */
    s_point *planes_abc = malloc(nFaces*sizeof(s_point));
    double *planes_d = malloc(nFaces*sizeof(double));
    for (int i=0; i<nFaces; i++) {
        /* Set the indices of the points defining the face  */
        for (int j=0, k=0; j<4; j++) {
            if (aVec[j] != i) {
                faces[i*3+k] = aVec[j];
                k++;
            }
        }
        /* Calculate and store the plane coefficients of the face */
        s_point plane[3];
        vertices_face(in_vertices, faces, i, plane);
        plane_equation_from_points(plane, &planes_abc[i], &planes_d[i]);       
    }

    /* Ensure that faces are correctly oriented */
    for (int k=0; k<4; k++) {
        s_point face_k[3];
        vertices_face(in_vertices, faces, k, face_k);

        /* Get the point that is not on the current face (point p) */
        int p = -1;
        for (int v = 0; v < 4; v++) {
            if (v != faces[k*3+0] && v != faces[k*3+1] && v != faces[k*3+2]) {
                p = v;
                break;
            }
        }
        assert(p != -1 && "No opposite vertex found");
        
        int v = orientation(face_k, in_vertices->p[p]);
        assert(v != 0);
        if (v < 0) {
            /* Flip face and plane */
            int tmp = faces[k*3 + 1];
            faces[k*3 + 1] = faces[k*3 + 2];
            faces[k*3 + 2] = tmp;

            planes_abc[k] = scale_point(planes_abc[k], -1);
            planes_d[k] = -planes_d[k];
        }
    }
    

    /* Coordinates of the center of the remaining point set */
    int N_aux_vertices = in_vertices->N - 4;
    s_points NONMALLOCED_aux_vertices = {.N = N_aux_vertices,
                                         .p = in_vertices->p + 4};
    s_point meanp = point_average(&NONMALLOCED_aux_vertices);

    /* Relative distance of points from the center */
    s_point span = span_points(in_vertices);  /* Used for normalizing, considering ALL input points */
    double *reldist2 = malloc(N_aux_vertices * sizeof(double));
    for (int i=0; i<N_aux_vertices; i++) {
        s_point translated = subtract_points(NONMALLOCED_aux_vertices.p[i], meanp);
        s_point scaled = {{{translated.x/span.x, translated.y/span.y, translated.z/span.z}}};
        reldist2[i] = norm_squared(scaled);
    }

    /* Sort from maximum to minimum relative distance */
    double *des_reldist2 = malloc(N_aux_vertices * sizeof(double));
    int *pleft = malloc(N_aux_vertices * sizeof(int));  /* This contains the new indices after sorting */
    sort_float(reldist2, des_reldist2, pleft, N_aux_vertices, 1);
    
    /* Initialize the vector of points left. The points with the larger relative
     distance from the center are scanned first. */
    /* We have to shift indices back by 4 to maintain original indexing. */
    for (int i=0; i<N_aux_vertices; i++) pleft[i] += 4;
    int num_pleft = N_aux_vertices;
    



    /* The main loop for the quickhull algorithm */

    int FUCKED = 0;
    double *points_cf = NULL;  // TODO what is this name?
    int *visible_ind = NULL, *visible_faces = NULL, *nonvisible_faces = NULL, *f0 = NULL;
    int *horizon = NULL;
    int *hVec = NULL, *hVec_mem_face = NULL, *pp = NULL;
    while (num_pleft > 0) {
        /* i is the first point of the points left */
        int i = pleft[0];
        s_point pi = in_vertices->p[i];

        
        /* Delete the point selected */
        num_pleft--;
        for(int j=0; j<num_pleft; j++) pleft[j] = pleft[j+1];
        if(num_pleft == 0) free(pleft);
        else pleft = realloc(pleft, num_pleft*sizeof(int));
        
        /* find visible faces */
        int num_visible_ind = 0;
        visible_ind = default_memory_resize(visible_ind, nFaces * sizeof(int));
        for (int j = 0; j < nFaces; ++j) {
            s_point face_pts[3];
            vertices_face(in_vertices, faces, j, face_pts);
            int o = orientation(face_pts, pi);
            if (o > 0) {  // Point is visible
                visible_ind[j] = 1;
                num_visible_ind++;
            } else if (o == 0) {
                printf("TODO: COPLANAR, VISIBLE OR NOT? OR SOMETHING OTHER?\n");
                exit(1);
            }
        }
        int num_nonvisible_faces = nFaces - num_visible_ind;
        
        /* proceed if there are any visible faces */
        if (num_visible_ind == 0) continue;

        /* Find visible face indices */
        visible_faces = default_memory_resize(visible_faces, num_visible_ind*sizeof(int));
        for (int j=0, k=0; j<nFaces; j++)
            if (visible_ind[j]==1) visible_faces[k++]=j;
        
        /* Find nonvisible faces */
        nonvisible_faces = default_memory_resize(nonvisible_faces, num_nonvisible_faces*3*sizeof(int));
        for (int j=0, k=0; j<nFaces; j++) {
            if (visible_ind[j]==0) {
                nonvisible_faces[k*3+0] = faces[j*3+0];
                nonvisible_faces[k*3+1] = faces[j*3+1];
                nonvisible_faces[k*3+2] = faces[j*3+2];
                k++;
            }
        }

        
        /* Create horizon (Ncount: number of edges of the horizon) */
        int Nhorizon = 0;
        f0 = default_memory_resize(f0, num_nonvisible_faces*3*sizeof(int));  // TODO what is f0?
        for (int j=0; j<num_visible_ind; j++) {
            /* visible face */
            int visible_id = visible_faces[j];
            int visible_vids[3] = {faces[visible_id*3], faces[visible_id*3+1], faces[visible_id*3+2]};
            ismember(nonvisible_faces, visible_vids, f0, num_nonvisible_faces*3, 3);

            /* Find nonvisible faces connected to this particular visible face */
            int *nonv_connected = NULL;
            int nonv_connected_len = 0;
            for (int k=0; k<num_nonvisible_faces; k++) {
                int f0_sum = f0[k*3+0] + f0[k*3+1] + f0[k*3+2];
                if (f0_sum == 2) {
                    nonv_connected_len++;
                    nonv_connected = default_memory_resize(nonv_connected, nonv_connected_len*sizeof(int));
                    nonv_connected[nonv_connected_len-1] = k;
                }
            }

            for (int k=0; k<nonv_connected_len; k++) {
                /* The boundary between the visible face v and the k(th) nonvisible face connected to the face v forms part of the horizon */
                int g[3] = {nonvisible_faces[nonv_connected[k]*3+0],
                            nonvisible_faces[nonv_connected[k]*3+1],
                            nonvisible_faces[nonv_connected[k]*3+2]};
                horizon = default_memory_resize(horizon, (Nhorizon+1)*2*sizeof(int));
                int h = 0;
                if( f0[nonv_connected[k]*3+0] ) horizon[Nhorizon*2+h++] = g[0];
                if( f0[nonv_connected[k]*3+1] ) horizon[Nhorizon*2+h++] = g[1];
                if( f0[nonv_connected[k]*3+2] ) horizon[Nhorizon*2+h++] = g[2];
                Nhorizon++;
            }
        }

        for (int j=0, l=0; j<nFaces; j++){
            if (!visible_ind[j]) {
                /* Delete visible faces */
                faces[l*3+0] = faces[j*3+0];
                faces[l*3+1] = faces[j*3+1];
                faces[l*3+2] = faces[j*3+2];

                /* Delete the corresponding plane coefficients of the faces */
                planes_abc[l*3].x = planes_abc[j*3].x;
                planes_abc[l*3].y = planes_abc[j*3].y;
                planes_abc[l*3].z = planes_abc[j*3].z;
                planes_d[l] = planes_d[j];
                l++;
            }
        }

        /* Update the number of faces */
        nFaces = nFaces - num_visible_ind;
        
        /* start is the first row of the new faces */
        int start = nFaces;
        
        /* Add faces connecting horizon to the new point */
        int N_newfaces = Nhorizon;
        int N_realloc_faces = nFaces + N_newfaces;
        if (N_realloc_faces > CH_MAX_NUM_FACES) N_realloc_faces = CH_MAX_NUM_FACES+1;

        faces = realloc(faces, N_realloc_faces*3*sizeof(int));
        planes_abc = realloc(planes_abc, N_realloc_faces*sizeof(s_point));
        planes_d = realloc(planes_d, N_realloc_faces*sizeof(double));
    
        for (int j=0; j<N_newfaces; j++) {
            faces[nFaces*3+0] = horizon[j*2+0];
            faces[nFaces*3+1] = horizon[j*2+1];
            faces[nFaces*3+2] = i;

            s_point plane[3];
            vertices_face(in_vertices, faces, nFaces, plane);
            plane_equation_from_points(plane, &planes_abc[nFaces], &planes_d[nFaces]);    

            nFaces++;
            if(nFaces > CH_MAX_NUM_FACES){
                FUCKED = 1;
                nFaces = 0;
                break;
            }
        }
        
        /* Orient each new face properly */
        hVec = default_memory_resize(hVec, nFaces*sizeof(int));
        hVec_mem_face = default_memory_resize(hVec_mem_face, nFaces*sizeof(int));
        for (int j=0; j<nFaces; j++) hVec[j] = j;

        for (int k=start; k<nFaces; k++) {
            int face_k[3] = {faces[k*3+0], faces[k*3+1], faces[k*3+2]};
            ismember(hVec, face_k, hVec_mem_face, nFaces, 3);

            int num_p = 0;
            for(int j=0; j<nFaces; j++) 
                if(!hVec_mem_face[j]) num_p++;

            pp = default_memory_resize(pp, num_p*sizeof(int));
            for (int j=0, l=0; j<nFaces; j++) {
                if (!hVec_mem_face[j]) {
                    pp[l++] = hVec[j];
                }
            }

            /* While new point is coplanar, choose another point */
            s_point face_k_vertices[3];
            vertices_face(in_vertices, faces, k, face_k_vertices);
            int index = 0;
            s_point opposite_p = in_vertices->p[pp[0]];
            int o = orientation(face_k_vertices, opposite_p);
            while ( o == 0) {
                assert(index+1 < num_p);
                opposite_p = in_vertices->p[pp[index++]];
                o = orientation(face_k_vertices, opposite_p);
            }

            /* Orient faces so that each point on the original simplex can't see the opposite face */
            if ( o < 0) {
                index--;
                /* Flip face and plane */
                int tmp = faces[k*3 + 1];
                faces[k*3 + 1] = faces[k*3 + 2];
                faces[k*3 + 2] = tmp;
                planes_abc[k] = scale_point(planes_abc[k], -1);
                planes_d[k] = -planes_d[k];

                /* Check */
                vertices_face(in_vertices, faces, k, face_k_vertices);
                assert(orientation(face_k_vertices, opposite_p) > 0);
            }
        }
        
        if(FUCKED){
            break;
        }
    }
    
    /* output */
    if (FUCKED) {
        *out_faces = NULL;
        *nOut_faces = 0;
    } else {
        *out_faces = malloc(nFaces*3*sizeof(int));
        memcpy(*out_faces, faces, nFaces*3*sizeof(int));
        *nOut_faces = nFaces;
    }
    
    /* clean-up */
    free(visible_ind);
    free(pp);
    free(horizon);
    free(f0);
    free(nonvisible_faces);
    free(visible_faces);
    free(hVec);
    free(hVec_mem_face);
    free(visible_ind);
    free(points_cf);
    free(reldist2);
    free(des_reldist2);
    free(pleft);
    free(faces);
    free(aVec);
    free(planes_abc);
    free(planes_d);
}





// void convhull_3d_export_obj(
//     ch_vertex* const vertices,
//     const int nVert,
//     int* const faces,
//     const int nFaces,
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
//         for (i = 0; i < nFaces; i++)
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
//     for (i = 0; i < nFaces; i++){
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
//         for (i = 0; i < nFaces; i++){
//             /* vertices are in same order as the faces, and normals are in order */
//             fprintf(obj_file, "f %u//%u %u//%u %u//%u\n",
//                     i*3 + 1, i + 1,
//                     i*3+1 + 1, i + 1,
//                     i*3+2 + 1, i + 1);
//         }
//     }
//     else {
//         /* just normals are in order  */
//         for (i = 0; i < nFaces; i++){
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


