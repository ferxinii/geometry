#include "convh.h"
#include "geometry.h"
#include "convhull_3d.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>


static ch_vertex *malloc_points_to_chvertex(const s_points *points)
{
    if (points->N <= 0) return NULL;
    ch_vertex *out = malloc((size_t)points->N * sizeof(ch_vertex));
    if (!out) { puts("ERROR: malloc_points_to_chvertex"); return NULL; }

    for (int ii=0; ii<points->N; ii++) {
        out[ii].x = (CH_FLOAT) points->p[ii].x;
        out[ii].y = (CH_FLOAT) points->p[ii].y;
        out[ii].z = (CH_FLOAT) points->p[ii].z;
    }
    return out;
}


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

        if (sign == 0) {  // Point is coplanar
            if (point_in_triangle_3d(pf, query)) {
                return -1;
            } else continue;
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


int mark_inside_convhull(const s_convhull *convh, const s_points query, int *out_mark)
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
    ch_vertex *ch_vertices = malloc_points_to_chvertex(points);

    s_convhull out;
    out.points = copy_points(points);

    convhull_3d_build(ch_vertices, points->N, &out.faces, &out.Nf);

    if (!out.faces) {
        fprintf(stderr, "convhull_from_points: Error in convhull_3d_build\n");
        out.Nf = 0;
        free_points(&out.points);
        out.fnormals = NULL;
    } else {
        initialize_normals_convhull(&out); 
    }

    free(ch_vertices);
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


void convh_get_face(const s_convhull *convh, int id, s_point *out)
{
    out[0] = convh->points.p[convh->faces[id*3+0]];
    out[1] = convh->points.p[convh->faces[id*3+1]];
    out[2] = convh->points.p[convh->faces[id*3+2]];
}
