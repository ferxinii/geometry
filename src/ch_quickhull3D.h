#ifndef GEOMETRY_CH_QUICKHULL3D_H
#define GEOMETRY_CH_QUICKHULL3D_H

#include "points.h"
#include <stdbool.h>

#define CH_MAX_NUM_FACES 10000
#define CH_N_INIT_INT_LIST 50  /* Memory buffers double memory when needed */

int quickhull_3d(const s_points *in_vertices, double EPS_degenerate, double TOL_dup, bool buff_isused[in_vertices->N], int **out_faces, int *N_out_faces);
/* Returns:
 *     -1 computational ERROR (memory, ...)
 *     0 degeneracy ERROR
 *     1 OK 
*/

#endif
