#include "convh.h"
#include "quickhull3d.h"
#include "geometry.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>


static void initialize_normals_convhull(s_convhull *convh)  // Unnormalized
{
    s_point ch_CM = point_average(&convh->points);

    convh->fnormals = malloc(sizeof(s_point) * convh->Nf);
    assert(convh->fnormals != NULL);

    for (int ii = 0; ii < convh->Nf; ii++) {
        s_point v0 = convh->points.p[convh->faces[ii * 3 + 0]];
        s_point v1 = convh->points.p[convh->faces[ii * 3 + 1]];
        s_point v2 = convh->points.p[convh->faces[ii * 3 + 2]];

        s_point n = cross_prod(subtract_points(v1, v0), subtract_points(v2, v0));

        s_point verts_face[3] = { v0, v1, v2 };
        int o = orientation(verts_face, ch_CM);
        // if (o == 0) {
        //     printf("%f, %f, %f\n", verts_face[0].x, verts_face[0].y, verts_face[0].z);
        //     printf("%f, %f, %f\n", verts_face[1].x, verts_face[1].y, verts_face[1].z);
        //     printf("%f, %f, %f\n", verts_face[2].x, verts_face[2].y, verts_face[2].z);
        //     printf("CM: %f, %f, %f\n", ch_CM.x, ch_CM.y, ch_CM.z);
        // }
        assert(o != 0);
        if (o < 0) {  // Correctly order face vertices (Right hand rule?)
            n.x = -n.x;  n.y = -n.y;  n.z = -n.z;
            int tmp = convh->faces[ii * 3 + 1];
            convh->faces[ii * 3 + 1] = convh->faces[ii * 3 + 2];
            convh->faces[ii * 3 + 2] = tmp;
        }

        convh->fnormals[ii] = n;
    }
}


int convex_hull_winding_valid(const s_convhull *convh)
{
	assert(convh && convh->Nf > 0 && convh->points.N > 0);

	s_point ch_CM = point_average(&convh->points);

	for (int ii=0; ii<convh->Nf; ii++) {
		s_point v0 = convh->points.p[convh->faces[3*ii + 0]];
		s_point v1 = convh->points.p[convh->faces[3*ii + 1]];
		s_point v2 = convh->points.p[convh->faces[3*ii + 2]];

		s_point n = cross_prod(subtract_points(v1, v0), subtract_points(v2, v0));
		s_point fc = {{{
			(v0.x + v1.x + v2.x) / 3.0,
			(v0.y + v1.y + v2.y) / 3.0,
			(v0.z + v1.z + v2.z) / 3.0
		}}};

		s_point to_center = subtract_points(ch_CM, fc);
		double dp = dot_prod(n, to_center);
		if (dp > 0.0) {
			fprintf(stderr, "Face %d likely inward-facing (dot = %g)\n", ii, dp);
            return 0;
		}
	}
    return 1;
}


int is_inside_convhull(const s_convhull *convh, s_point query)
{   
    // 1: inside, 0: outside, -1: in boundary
    assert(convh->Nf > 0 && "is_inside_convhull: Nf <= 0?");

    int prev_sign = 0;
    for (int f = 0; f < convh->Nf; ++f) {
        s_point pf[3] = { convh->points.p[convh->faces[3*f + 0]],
                          convh->points.p[convh->faces[3*f + 1]],
                          convh->points.p[convh->faces[3*f + 2]] };
        
        int sign = orientation(pf, query);

        if (sign == 0) {  // Point is coplanar, so it is either on the face or not
            if (in_triangle_3d(pf, query) != 0) {
                return -1;
            } else return 0;
        }

        // if we've already seen a non-zero sign, it must match
        if (prev_sign == 0) prev_sign = sign;
        else if (sign != prev_sign) return 0;  // outside!
    }

    assert(prev_sign != 0  && "Point is coplanar with all faces? Strange...");
    return 1;
}


static int inarray(const int *arr1, int N, int a)
{
    for (int ii=0; ii<N; ii++) {
        if (arr1[ii] == a) return 1;
    }
    return 0;
}


int is_in_boundary_convhull(const s_convhull *convh, int point_id)
{
    return inarray(convh->faces, convh->Nf * 3, point_id);
}


int mark_boundary_convhull(const s_convhull *convh, int out[convh->points.N])
{
    int count = 0;
    memset(out, 0, sizeof(int) * convh->points.N);
    for (int ii=0; ii<convh->points.N; ii++) {
        if (is_in_boundary_convhull(convh, ii)) { 
            out[ii] = 1;       
            count++;
        }
    }
    return count;
}


s_points boundary_convhull(const s_convhull *convh)
{
    int is_boundary[convh->points.N];
    int count_boundary = mark_boundary_convhull(convh, is_boundary);

    s_points boundary = { .N = count_boundary, 
                          .p = malloc(sizeof(s_point) * count_boundary)};

    int jj = 0;
    for (int ii=0; ii<convh->points.N; ii++) {
        if (is_boundary[ii]) 
            boundary.p[jj++] = convh->points.p[ii];
    }

    return boundary;
}



int points_inside_convhull(const s_convhull *convh, const s_points query, int out_mark[query.N])
{
    memset(out_mark, 0, sizeof(int) * query.N);
    int count = 0;
    for (int ii=0; ii<query.N; ii++) {
        if (is_inside_convhull(convh, query.p[ii])) {
            out_mark[ii] = 1;
            count++;
        }
    }
    return count;
}


s_point random_point_inside_convhull(const s_convhull *convh, s_point min, s_point max)
{
    int MAX_IT = 10000;
    int it = 0;
    s_point out = random_point_uniform_3d(min, max);
    while (is_inside_convhull(convh, out) != 1) {
        out = random_point_uniform_3d(min, max);
        assert(it < MAX_IT && "Reached maximum iters looking for point inside convhull.");
        it++;
    }
    return out;
}


double volume_convhull(const s_convhull *convh)
{
    if (convh->Nf == 0) return 0;

    double vol = 0;
    for (int ii=0; ii<convh->Nf; ii++) {
        double Nx = convh->fnormals[ii].x;
        vol += Nx * (convh->points.p[convh->faces[ii*3 + 0]].x +
                     convh->points.p[convh->faces[ii*3 + 1]].x +
                     convh->points.p[convh->faces[ii*3 + 2]].x);
    }
    return vol / 6;
}


double volume_convhull_from_points(const s_points *points)
{
    s_convhull convh = convhull_from_points(points);   
    if (convh.Nf == 0) return 0;
    double volume = volume_convhull(&convh);
    free_convhull(&convh);
    return volume;
}


s_convhull convhull_from_points(const s_points *points)
{
    s_convhull out = {0};
    if (points->N == 0) return out;

    out.points = remove_duplicate_points(points, 1e-12);

    
    if (quickhull_3d(&out.points, 0, &out.faces, &out.Nf) == -1) {
        fprintf(stderr, "convhull_from_points: Error in quickhull_3d\n");
        out.Nf = 0;
        free_points(&out.points);
        out.fnormals = NULL;
    } else {
        initialize_normals_convhull(&out); 
    }


    return out;
}


void free_convhull(s_convhull *convh)
{
    free((void*)convh->points.p);
    free(convh->faces);
    free(convh->fnormals);
    memset(convh, 0, sizeof(s_convhull));
}


s_convhull copy_convhull(const s_convhull *in)
{
    s_convhull out;
    out.points = copy_points(&in->points);
    out.Nf = in->Nf;
    out.faces = malloc(sizeof(int) * 3 * in->Nf);
    memcpy(out.faces, in->faces, sizeof(int) * 3 * in->Nf);
    out.fnormals = malloc(sizeof(s_point) * in->Nf);
    memcpy(out.fnormals, in->fnormals, sizeof(s_point) * in->Nf);
    return out;
}


void convh_get_face(const s_convhull *convh, int id, s_point out[3])
{
    out[0] = convh->points.p[convh->faces[id*3+0]];
    out[1] = convh->points.p[convh->faces[id*3+1]];
    out[2] = convh->points.p[convh->faces[id*3+2]];
}






int mark_faces_incident_to_vertex(const s_convhull *C, int vid, int out[C->Nf])
{
    int count = 0;
    memset(out, 0, sizeof(int) * C->Nf);
    for (int ii=0; ii<C->Nf; ii++) {
        if (inarray(&C->faces[ii*3], 3, vid)) {
            out[ii] = 1;
            count++;
        }
    }
    return count;
}


typedef struct edge_test {
    int v[2];
} s_edge_test;


static int edge_exists(const s_edge_test *edges, int N, s_edge_test query)
{
    for (int ii=0; ii<N; ii++) {
        if ( (edges[ii].v[0] == query.v[0] && edges[ii].v[1] == query.v[1]) ||
             (edges[ii].v[0] == query.v[1] && edges[ii].v[1] == query.v[0]) )
            return 1;
    }
    return 0;
}


static void extract_edges_face(const s_convhull *C, int fid, s_edge_test out[3])
{
    out[0].v[0] = C->faces[fid*3 + 0];
    out[0].v[1] = C->faces[fid*3 + 1];

    out[1].v[0] = C->faces[fid*3 + 1];
    out[1].v[1] = C->faces[fid*3 + 2];

    out[2].v[0] = C->faces[fid*3 + 2];
    out[2].v[1] = C->faces[fid*3 + 0];
}


static void extract_edges_to_test(const s_convhull *C, int vid, s_edge_test **out, int *Nout)
{
    int marked_faces[C->Nf];
    int Nmarked = mark_faces_incident_to_vertex(C, vid, marked_faces);
    
    s_edge_test *edges = malloc(sizeof(s_edge_test) * Nmarked * 3);
    int Nedges = 0;

    for (int ii=0; ii<C->Nf; ii++) {
        if (!marked_faces[ii]) continue;

        s_edge_test edges_face[3];
        extract_edges_face(C, ii, edges_face);
        if (!edge_exists(edges, Nedges, edges_face[0])) 
            edges[Nedges++] = edges_face[0];
        if (!edge_exists(edges, Nedges, edges_face[1])) 
            edges[Nedges++] = edges_face[1];
        if (!edge_exists(edges, Nedges, edges_face[2])) 
            edges[Nedges++] = edges_face[2];
    }

    *Nout = Nedges;
    *out = realloc(edges, sizeof(s_edge_test) * Nedges);
}


static s_points clip_faces_incident_to_vertex_with_plane(const s_convhull *C, int vid, const s_point plane[3])
{
    int Nedges;
    s_edge_test *edges;
    extract_edges_to_test(C, vid, &edges, &Nedges);
    if (Nedges == 0) {
        s_points empty = {0, NULL};
        return empty;
    }

    s_points out = {.N = 0, 
                    .p = malloc(sizeof(s_point) * Nedges * 2)};

    for (int ii=0; ii<Nedges; ii++) {
        s_point segment[2] = {C->points.p[edges[ii].v[0]], 
                              C->points.p[edges[ii].v[1]]};
        s_point intersections[2];
        int Nintersections = segment_plane_intersection(segment, plane, intersections);
        for (int jj=0; jj<Nintersections; jj++) {
            assert(out.N < Nedges*2 && "Too many cliping vertices");
            out.p[out.N++] = intersections[jj];
        }
    }

    out.p = realloc(out.p, sizeof(s_point) * out.N);
    free(edges);
    return out;
}


static s_convhull clip_convhull_against_halfspace(const s_convhull *C, const s_point plane_ordered[3])
{
    int inside[C->points.N];
    int Nin = points_inside_halfspace(plane_ordered, C->points, inside);
    int Nout = C->points.N - Nin;

    if (Nin == 0) return (s_convhull){0};
    if (Nout == 0) return copy_convhull(C);

    s_points array_points[Nout];
    int jj = 0, N_added_vertices = 0;
    for (int ii=0; ii<C->points.N; ii++) {
        if (inside[ii] == 0) {
            array_points[jj] = clip_faces_incident_to_vertex_with_plane(C, ii, plane_ordered);
            N_added_vertices += array_points[jj].N;
            jj++;
        }
    }

    int Nnew = Nin + N_added_vertices;
    s_points new_points = {.N = Nnew,
                           .p = malloc(sizeof(s_point) * Nnew)};
    int count = 0;
    for (int ii=0; ii<C->points.N; ii++) {
        if (inside[ii] == 1) new_points.p[count++] = C->points.p[ii];
    }
    for (int ii=0; ii<Nout; ii++) {
        for (int jj=0; jj<array_points[ii].N; jj++) {
            new_points.p[count++] = array_points[ii].p[jj];
        }
        free_points(&array_points[ii]);
    }

    s_convhull new_C = convhull_from_points(&new_points);
    free_points(&new_points);

    return new_C;
}


s_convhull intersection_convhulls(const s_convhull *A, const s_convhull *B)
{
    s_convhull I = copy_convhull(A);
    for (int ii=0; ii<B->Nf; ii++) {
        s_point face[3]; 
        convh_get_face(B, ii, face);

        s_convhull tmp = clip_convhull_against_halfspace(&I, face);
        free_convhull(&I);
        I = tmp;
        if (I.points.N == 0) break;
    }
    return I;
}


static void bisector_plane(const s_convhull *A, const s_convhull *B, const s_convhull *I, s_point out[3])
{
    s_point cA = point_average(&A->points);
    s_point cB = point_average(&B->points);
    s_point cI = point_average(&I->points);
    s_point n = subtract_points(cB, cA);
    assert(norm(n) > 1e-9 && "Centroids are coincidental!");
    s_point normal_n = normalize_3d(n);

    int ref_coord = coord_with_smallest_component_3d(n);
    s_point ref = (ref_coord == 0) ?   (s_point){{{1,0,0}}} :
                  ( (ref_coord == 1) ? (s_point){{{0,1,0}}} :
                                       (s_point){{{0,0,1}}} );
    s_point u = normalize_3d(cross_prod(normal_n, ref));
    s_point v = normalize_3d(cross_prod(normal_n, u));
    s_point t1 = sum_points(cI, u);
    s_point t2 = sum_points(cI, v);

    out[0] = cI;  out[1] = t1;  out[2] = t2;
}


int remove_intersection_convhulls(s_convhull *A, s_convhull *B)
{
    s_convhull I = intersection_convhulls(A, B);
    if (volume_convhull(&I) < 1e-9) {
        free_convhull(&I);
        return 0;
    }

    s_point plane[3];
    bisector_plane(A, B, &I, plane);
    
    // Clip I against the plane
    s_convhull IA = clip_convhull_against_halfspace(&I, plane);

    s_point tmp = plane[0];
    plane[0] = plane[1];  plane[1] = tmp;
    s_convhull IB = clip_convhull_against_halfspace(&I, plane);
    assert(IA.points.N > 0 && IB.points.N > 0 &&
           "Invalid clipping of intersection with plane.");


    
    // Remove intersection
    int A_mark_inside_B[A->points.N];
    int N_A_inside_B = points_inside_convhull(B, A->points, A_mark_inside_B);
    int N_newA = A->points.N - N_A_inside_B + IA.points.N;
    s_points p_newA = {.N = N_newA,
                       .p = malloc(sizeof(s_point) * N_newA)};
    int jj=0;
    for (int ii=0; ii<A->points.N; ii++) {
        if (A_mark_inside_B[ii] == 0) p_newA.p[jj++] = A->points.p[ii];
    }
    for (int ii=0; ii<IA.points.N; ii++) {
        p_newA.p[jj++] = IA.points.p[ii];
    }
    s_points p_newA_noduplicate = remove_duplicate_points(&p_newA, 1e-9);
    s_convhull new_A = convhull_from_points(&p_newA_noduplicate);
    free_points(&p_newA);
    free_points(&p_newA_noduplicate);

    int B_mark_inside_A[B->points.N];
    int N_B_inside_A = points_inside_convhull(A, B->points, B_mark_inside_A);
    int N_newB = B->points.N - N_B_inside_A + IB.points.N;
    s_points p_newB = {.N = N_newB,
                       .p = malloc(sizeof(s_point) * N_newB)};
    jj=0;
    for (int ii=0; ii<B->points.N; ii++) {
        if (B_mark_inside_A[ii] == 0) p_newB.p[jj++] = B->points.p[ii];
    }
    for (int ii=0; ii<IB.points.N; ii++) {
        p_newB.p[jj++] = IB.points.p[ii];
    }
    s_points p_newB_noduplicate = remove_duplicate_points(&p_newB, 1e-9);
    s_convhull new_B = convhull_from_points(&p_newB_noduplicate);
    free_points(&p_newB);
    free_points(&p_newB_noduplicate);

    free_convhull(A);
    *A = new_A;
    free_convhull(B);
    *B = new_B;

    free_convhull(&I);
    free_convhull(&IA);
    free_convhull(&IB);

    return 1;
}

