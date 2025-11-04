#ifndef CH_QUICKHULL3D_H
#define CH_QUICKHULL3D_H


#include "geometry.h"

#define CH_MAX_NUM_FACES 10000
#define CH_N_INIT_INT_LIST 50

int quickhull_3d(const s_points *in_vertices, double min_face_area, int **out_faces, int *N_out_faces);

/* Returns:
       -1 if error
       1 if output hull is exact (All faces are big enough)
       0 if output hull is non-exact (ignored any face with too small area) 
*/

#endif
