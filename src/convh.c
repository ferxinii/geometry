#include "convh.h"
#include "ch_quickhull3D.h"
#include "points.h"
#include "gtests.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <float.h>


int convhull_is_valid(const s_convh *convh) 
{
    if (!convh) return 0;
    if (convh->Nf == 0) return 0;
    else return 1;
}


static int initialize_normals_convhull(s_convh *convh)
{   /* 0 if error, 1 if OK. UNNORMALIZED! */
    s_point ch_CM = point_average(&convh->points);

    convh->fnormals = malloc(sizeof(s_point) * convh->Nf);
    if (!convh->fnormals) return 0;

    for (int ii = 0; ii < convh->Nf; ii++) {
        s_point v0 = convh->points.p[convh->faces[ii * 3 + 0]];
        s_point v1 = convh->points.p[convh->faces[ii * 3 + 1]];
        s_point v2 = convh->points.p[convh->faces[ii * 3 + 2]];

        s_point n = cross_prod(subtract_points(v1, v0), subtract_points(v2, v0));

        s_point verts_face[3] = { v0, v1, v2 };
        int o = orientation_robust(verts_face, ch_CM);
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
    return 1;
}


int convhull_from_points(const s_points *points, double EPS_degenerate, double TOL, s_convh *out)
{   /* TOL both for deduping points and to establish min_face_area */
    /* Returns  1 if OK,
                0 if could not construct convhull due to the nature of the points,
                -1 if some other error, perhaps memory related
    */
    if (points->N == 0) { *out = convhull_NAN; return 0; }

    s_points points_nodup = copy_points_remove_duplicates(points, TOL);
    int *isused = malloc(points->N * sizeof(int));
    if (!points_is_valid(&points_nodup) || !isused) goto error;

    int i = quickhull_3d(&points_nodup, EPS_degenerate, isused, &out->faces, &out->Nf); 
    if (i == -2) return 0;
    if (i == -1) return -1;
    
    /* Update out.points so that only used points remain */
    int Nnew = 0;
    for (int ii=0; ii<points_nodup.N; ii++) {
        if (isused[ii]) {
            isused[ii] = Nnew;  /* Overwrite in isused the "new index" */
            points_nodup.p[Nnew++] = points_nodup.p[ii];
        } else {
            isused[ii] = -1;
        }
    }
    out->points.N = Nnew;
    out->points.p = realloc(points_nodup.p, Nnew * sizeof(s_point));

    /* Update face indices */
    for (int f=0; f<out->Nf; f++) {
        for (int v=0; v<3; v++) {
            int old = out->faces[f*3+v];
            int newi = isused[old];
            assert(newi >= 0 && "Face refers to unused vertex!");
            out->faces[f*3+v] = newi;
        }
    }

    if (!initialize_normals_convhull(out)) goto error;
    free(isused);
    return 1;

    error:
        fprintf(stderr, "Error in 'convhull_from_points'.\n");
        free(isused);
        free_points(&points_nodup);
        *out = convhull_NAN;
        return -1;
}


int convhull_from_csv(const char *filename, double EPS_degenerate, double TOL, s_convh *out)
{
    s_points p = read_points_from_csv(filename);
    if (!points_is_valid(&p)) { *out = convhull_NAN; return -1; }
    
    int i = convhull_from_points(&p, EPS_degenerate, TOL, out);
    free_points(&p);
    return i;
}


void free_convhull(s_convh *convh)
{
    free((void*)convh->points.p);
    free(convh->faces);
    free(convh->fnormals);
    memset(convh, 0, sizeof(s_convh));
}


s_convh copy_convhull(const s_convh *in)
{
    s_convh out = {0};
    out.points = copy_points(&in->points);
    out.Nf = in->Nf;
    out.faces = malloc(sizeof(int) * 3 * in->Nf);
    if (!out.faces) goto error;
    memcpy(out.faces, in->faces, sizeof(int) * 3 * in->Nf);
    out.fnormals = malloc(sizeof(s_point) * in->Nf);
    if (!out.fnormals) goto error;
    memcpy(out.fnormals, in->fnormals, sizeof(s_point) * in->Nf);
    return out;

    error:
        if (out.faces) free(out.faces);
        if (out.fnormals) free(out.fnormals);
        return convhull_NAN;
}


int convex_hull_winding_valid(const s_convh *convh)
{
    if (!convhull_is_valid(convh)) return 0;

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


static e_geom_test test_point_in_convhull_robust(const s_convh *convh, s_point query, double EPS_degenerate)
{   
    if (!convhull_is_valid(convh)) return TEST_ERROR;

    int prev_sign = 0;
    for (int f = 0; f < convh->Nf; ++f) {
        s_point pf[3] = { convh->points.p[convh->faces[3*f + 0]],
                          convh->points.p[convh->faces[3*f + 1]],
                          convh->points.p[convh->faces[3*f + 2]] };
        
        int sign = orientation_robust(pf, query);

        if (sign == 0) {  // Point is coplanar, but inside face?
            e_geom_test intri = test_point_in_triangle_3D(pf, query, EPS_degenerate, 0);
            if (intri == TEST_IN || intri == TEST_BOUNDARY) return TEST_BOUNDARY;
            else return TEST_OUT;
        }

        /* if we've already seen a non-zero sign, it must match */
        if (prev_sign == 0) prev_sign = sign;
        else if (sign != prev_sign) return TEST_OUT;  // Outside!
    }

    if (prev_sign == 0) return TEST_ERROR;  /* Point is coplanar with all faces? Strange... */
    return TEST_IN;
}

static int sign_double(double x)
{
    if (x>0) return 1;
    else if (x<0) return 0;
    else return 0;
}

e_geom_test test_point_in_convhull(const s_convh *C, s_point query, double EPS_degenerate, double TOL_boundary)
{
    if (!convhull_is_valid(C)) return TEST_ERROR;
    
    if (TOL_boundary == 0) return test_point_in_convhull_robust(C, query, EPS_degenerate);
    int on_boundary = 0, sign_ref = 0;

    for (int f = 0; f < C->Nf; ++f) {
        s_point face[3] = {C->points.p[C->faces[3*f + 0]],
                           C->points.p[C->faces[3*f + 1]],
                           C->points.p[C->faces[3*f + 2]]};

        double s = signed_distance_point_to_plane(query, face, EPS_degenerate);
        if (isnan(s)) continue;
        if (s <= TOL_boundary) {
            e_geom_test inside_face = test_point_in_triangle_3D(face, query, EPS_degenerate, TOL_boundary);
            if (inside_face == TEST_IN || inside_face == TEST_BOUNDARY) on_boundary = 1;
            continue;
        }

        int sign = sign_double(s); 
        if (sign_ref == 0) sign_ref = sign;
        else if (sign != sign_ref) return TEST_OUT;  // Point is outside
    }

    if (on_boundary) return TEST_BOUNDARY;
    if (sign_ref == 0) return TEST_ERROR; // degenerate hull?
    return TEST_IN;
}


s_points_test test_points_in_convhull(const s_convh *convh, const s_points *query, double EPS_degenerate, double TOL_boundary, e_geom_test buff[query->N])
{
    if (buff == NULL) buff = malloc(convh->points.N * sizeof(e_geom_test));

    int Nin = 0, Nout = 0, Nbdy = 0, Nerr = 0;
    for (int ii=0; ii<query->N; ii++) {
        e_geom_test test = test_point_in_convhull(convh, query->p[ii], EPS_degenerate, TOL_boundary);
        if (test == TEST_IN) { buff[ii] = TEST_IN; Nin++; }
        else if (test == TEST_OUT) { buff[ii] = TEST_OUT; Nout++; }
        else if (test == TEST_BOUNDARY) { buff[ii] = TEST_BOUNDARY; Nbdy++; }
        else { buff[ii] = TEST_ERROR; Nerr++; }
    }
    return (s_points_test){ .Nin = Nin, .Nout = Nout, .Nbdy = Nbdy, .Nerr = Nerr, .indicator = buff };
}


static int inarray(const int *arr1, int N, int a)
{
    for (int ii=0; ii<N; ii++)
        if (arr1[ii] == a) return 1;

    return 0;
}


int convhull_id_in_boundary(const s_convh *convh, int point_id)
{
    return inarray(convh->faces, convh->Nf * 3, point_id);
}


int mark_boundary_convhull(const s_convh *convh, int out[convh->points.N])
{
    memset(out, 0, sizeof(int) * convh->points.N);
    for (int ii=0; ii<convh->Nf; ii++) {
        out[convh->faces[ii*3+0]] = 1;       
        out[convh->faces[ii*3+1]] = 1;       
        out[convh->faces[ii*3+2]] = 1;       
    }

    int count = 0;
    for (int ii=0; ii<convh->points.N; ii++) 
        if (out[ii]) count++;
   
    return count;
}


s_points points_boundary_convhull(const s_convh *convh, int buff_isboundary[convh->points.N])
{
    int count_boundary = mark_boundary_convhull(convh, buff_isboundary);

    s_points boundary = { .N = count_boundary, 
                          .p = malloc(sizeof(s_point) * count_boundary)};
    if (!boundary.p) return points_NAN;

    int jj = 0;
    for (int ii=0; ii<convh->points.N; ii++)
        if (buff_isboundary[ii]) boundary.p[jj++] = convh->points.p[ii];

    return boundary;
}


s_point random_point_inside_convhull(const s_convh *convh, double EPS_degenerate, s_point min, s_point max)
{
    if (!point_is_valid(min) || !point_is_valid(max)) 
        bounding_box_points(&convh->points, &min, &max);

    int MAX_IT = 10000;
    int it = 0;
    s_point out = random_point_uniform_3D(min, max);
    while (test_point_in_convhull(convh, out, EPS_degenerate, 0) != TEST_IN) {
        out = random_point_uniform_3D(min, max);
        assert(it < MAX_IT && "Reached maximum iters looking for point inside convhull.");
        it++;
    }
    return out;
}


s_point closest_point_on_convhull_boundary(const s_convh *convh, s_point query, double EPS_degenerate) 
{
    s_point out_point = {{{0, 0, 0}}};
    double best_d2 = DBL_MAX;

    for (int f=0; f<convh->Nf; f++) {
        s_point triangle[3];
        convh_get_face(convh, f, triangle);
        s_point tmp = closest_point_on_triangle(triangle, EPS_degenerate, query);
        double d2 = distance_squared(query, tmp);
        if ( d2 < best_d2) {
            out_point = tmp;
            best_d2 = d2; 
        }
    }

    if (best_d2 != DBL_MAX) return out_point;
    else return point_NAN;
}


double volume_convhull(const s_convh *convh)
{
    if (!convhull_is_valid(convh)) return 0;

    double vol = 0;
    for (int ii=0; ii<convh->Nf; ii++) {
        double Nx = convh->fnormals[ii].x;
        vol += Nx * (convh->points.p[convh->faces[ii*3 + 0]].x +
                     convh->points.p[convh->faces[ii*3 + 1]].x +
                     convh->points.p[convh->faces[ii*3 + 2]].x);
    }
    return vol / 6;
}


s_point convhull_volume_centroid(const s_convh *convh, double EPS_degenerate)
{
    if (!convhull_is_valid(convh)) return point_NAN;

    s_point vertex_centroid = point_average(&convh->points);
    double Vtot = 0;
    s_point pweighted = {0};
    for (int ii=0; ii<convh->Nf; ii++) {
        s_point tetra[4] = {{0}, {0}, {0}, vertex_centroid};
        convh_get_face(convh, ii, tetra);  /* Writes to first 3 entries of tetra */
        s_point c_tet = {{{ (tetra[0].x + tetra[1].x + tetra[2].x + tetra[3].x) / 4.0,
                            (tetra[0].y + tetra[1].y + tetra[2].y + tetra[3].y) / 4.0, 
                            (tetra[0].z + tetra[1].z + tetra[2].z + tetra[3].z) / 4.0 }}}; 
        double vol_tet = signed_volume_tetra(tetra);
        Vtot += vol_tet;
        pweighted = sum_points(pweighted, scale_point(c_tet, vol_tet));
    }
    if (Vtot < EPS_degenerate) return point_NAN;
    return scale_point(pweighted, 1.0 / Vtot);
}


void convh_get_face(const s_convh *convh, int id, s_point out[3])
{
    out[0] = convh->points.p[convh->faces[id*3+0]];
    out[1] = convh->points.p[convh->faces[id*3+1]];
    out[2] = convh->points.p[convh->faces[id*3+2]];
}


int mark_faces_incident_to_vertex(const s_convh *C, int vid, int out[C->Nf])
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


void write_convhull_to_m(const s_convh *convh, const char *filename)
{
    FILE *file = fopen(filename, "wt");
    
    /* save face indices and vertices for verification in matlab: */
    fprintf(file, "vertices = [\n");
    for (int i=0; i<convh->points.N; i++)
        fprintf(file, "%f, %f, %f;\n", convh->points.p[i].x, convh->points.p[i].y, convh->points.p[i].z);

    fprintf(file, "];\n\n\n");

    fprintf(file, "faces = [\n");
    for (int i=0; i < convh->Nf; i++) {
        fprintf(file, " %u, %u, %u;\n", convh->faces[3*i+0]+1, convh->faces[3*i+1]+1, convh->faces[3*i+2]+1);
    }

    fprintf(file, "];\n\n\n");

    fclose(file);
}

