#ifndef GEOMETRY_CH_QUICKHULL3D_H
#define GEOMETRY_CH_QUICKHULL3D_H

#include "points.h"

#define CH_MAX_NUM_FACES 10000
#define CH_N_INIT_INT_LIST 50  /* Memory buffers double memory when needed */

int quickhull_3d(const s_points *in_vertices, double min_face_area, int buff_isused[in_vertices->N], int **out_faces, int *N_out_faces);
/* Returns:
 *     -2 if initalization error. in_vertices cant make a convex hull
 *     -1 if error (memory, ...)
 *     1 if output hull is exact (All faces are big enough)
 *     0 if output hull is non-exact (ignored any face with too small area) 
*/

#endif
