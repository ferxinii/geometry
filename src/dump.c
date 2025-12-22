#include <math.h>
#include "points.h"

static double col_norm(double X[4][3], int p) 
{
    double s = 0.0;
    for (int i = 0; i < 4; ++i) {
        double x = X[i][p];
        s += x * x;
    }
    return sqrt(s);
}

static int one_sided_jacobi_sweeps_4x3(double X[4][3], double V[3][3], int MAX_SWEEPS, double TOL_convergence)
{
    int pairs[3][2] = {{0,1},{0,2},{1,2}};
    int converged = 0;
    for (int sweep=0; sweep<MAX_SWEEPS; sweep++) {
        double max_off = 0;
        for (int pi=0; pi<3; pi++) {
            int p = pairs[pi][0], q = pairs[pi][1];

            /* compute 2x2 Gram entries */
            double app = 0, aqq = 0, apq = 0;
            for (int i=0; i<4; i++) {
                double xp = X[i][p], xq = X[i][q];
                app += xp * xp;
                aqq += xq * xq;
                apq += xp * xq;
            }
            double abs_apq = fabs(apq);
            if (abs_apq > max_off) max_off = abs_apq;
            /* if already nearly orthogonal or one column is zero, skip */
            if (abs_apq <= TOL_convergence) continue;

            /* compute stable Jacobi rotation that diagonalizes [app apq; apq aqq] */
            double tau = (aqq - app) / (2.0 * apq);
            double t;
            if (tau >= 0) t = 1.0 / (tau + sqrt(1 + tau * tau));
            else t = -1.0 / (-tau + sqrt(1 + tau * tau));
            double c = 1.0 / sqrt(1 + t * t);
            double s = t * c;

            /* Apply rotation to columns p and q of X: (for each row) */
            for (int i=0; i<4; i++) {
                double xp = X[i][p], xq = X[i][q];
                X[i][p] =  c * xp - s * xq;
                X[i][q] =  s * xp + c * xq;
            }

            /* Accumulate rotation into V (right singular vectors) */
            for (int i=0; i<3; i++) {
                double vip = V[i][p], viq = V[i][q];
                V[i][p] =  c * vip - s * viq;
                V[i][q] =  s * vip + c * viq;
            }
        }  
        if (max_off <= TOL_convergence) { converged = 1; break; }
    }  
    if (converged) return 1;
    else return 0;
}


int fit_plane_svd_4x3(const s_point pts[4], double TOL_convergence, double EPS_degenerate, s_point *out_normal, s_point *out_centroid)
{   /* 0: DEGENERATE OR ERROR, 1: OK */
    const int MAX_SWEEPS = 100;         /* number of Jacobi sweeps (usually 5-30 needed) */

    *out_centroid = point_average(&(s_points){ .N=4, .p=(s_point*)pts });

    /* 1) build centered X (4x3) */
    double X[4][3], max_abs = 0.0;
    for (int i=0; i<4; i++)
        for (int j=0; j<3; j++) {
            X[i][j] = pts[i].coords[j] - (*out_centroid).coords[j];
            double a = fabs(X[i][j]);
            if (a > max_abs) max_abs = a;
        }
    if (max_abs <= EPS_degenerate) return 0;  /* All points identical -> zero matrix */

    double scale = 1;  /* Scale X to improve numeric behavior (so max_abs = 1) */
    if (max_abs != 0) {
        scale = max_abs;
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 3; ++j)
                X[i][j] /= scale;
    }

    double V[3][3];  /* Initialize V = Id (3x3) */
    for (int i=0; i<3; i++)
        for (int j=0; j<3; j++)
            V[i][j] = (i == j) ? 1.0 : 0.0;
    
    if (!one_sided_jacobi_sweeps_4x3(X, V, MAX_SWEEPS, TOL_convergence)) return 0;

    /* Compute singular values (norms of columns) and find smallest */
    double sigma[3];
    for (int j=0; j<3; j++) sigma[j] = col_norm(X, j) * scale;  /* Scale back ! */

    double sigma_max = sigma[0], sigma_min = sigma[0];
    int idx_min = 0;
    for (int j=1; j<3; j++) {
        if (sigma[j] > sigma_max) sigma_max = sigma[j];
        if (sigma[j] < sigma_min) { sigma_min = sigma[j]; idx_min = j; }
    }
    if (sigma_max < EPS_degenerate) return 0;  /* all points identical after centering */

    /* 6) right singular vector corresponding to smallest singular value
       is the plane normal (column idx_min of V). Ensure it's normalized. */
    double nrm = sqrt(V[0][idx_min]*V[0][idx_min]
                    + V[1][idx_min]*V[1][idx_min]
                    + V[2][idx_min]*V[2][idx_min]);
    if (nrm < EPS_degenerate) return 0; 
    (*out_normal).x = V[0][idx_min] / nrm;
    (*out_normal).y = V[1][idx_min] / nrm;
    (*out_normal).z = V[2][idx_min] / nrm;

    return 1; 
}

int fit_circle_4_points_2D(const double p1[2], const double p2[2], const double p3[2], const double p4[2], double EPS_degenerate, double out_center[2], double *out_radius)
{
	double A[4][3];
	double b[4];
	const double *pts[4] = { p1, p2, p3, p4 };

	for (int ii=0; ii<4; ii++) {  /* Linear system */
		A[ii][0] = pts[ii][0];
		A[ii][1] = pts[ii][1];
		A[ii][2] = 1.0;
		b[ii] = -(pts[ii][0]*pts[ii][0] + pts[ii][1]*pts[ii][1]);
	}

    /* Modified Gram-Schmidt: factor A = Q * R */
	double Q[4][3];
	double R[3][3] = {0};
	/* Copy A into Q initially (we will orthonormalize columns) */
	for (int i = 0; i < 4; ++i)
		for (int j = 0; j < 3; ++j)
			Q[i][j] = A[i][j];

	for (int k = 0; k < 3; ++k) {
		double rkk = 0.0;  /* r_{kk} = norm of column k */
		for (int i = 0; i < 4; ++i)  rkk += Q[i][k] * Q[i][k];
		rkk = sqrt(rkk);

		if (rkk < EPS_degenerate || !isfinite(rkk)) return 0;
		R[k][k] = rkk;
		
		for (int i = 0; i < 4; ++i) Q[i][k] /= rkk;  /* normalize column k */

		for (int j = k + 1; j < 3; ++j) {  /* orthogonalize remaining columns */
			double rkj = 0.0;
			for (int i = 0; i < 4; ++i) rkj += Q[i][k] * Q[i][j];
			R[k][j] = rkj;
			for (int i = 0; i < 4; ++i) Q[i][j] -= rkj * Q[i][k];
		}
	}

	/* Compute y = Q^T * b  (length n) */
	double y[3];
	for (int k = 0; k < 3; ++k) {
		double s = 0.0;
		for (int i = 0; i < 4; ++i) s += Q[i][k] * b[i];
		y[k] = s;
	}

	/* Solve R * x = y  (R is upper-triangular n x n) via back-substitution */
	double xsol[3];
	for (int i = 2; i >= 0; --i) {
		double s = y[i];
		for (int j = i + 1; j < 3; ++j) s -= R[i][j] * xsol[j];
		if (fabs(R[i][i]) < EPS_degenerate) return 0; 
		xsol[i] = s / R[i][i];
	}

	/* Convert to center / radius */
	double xc = -xsol[0] * 0.5;
	double yc = -xsol[1] * 0.5;
	double r = sqrt(xc*xc + yc*yc - xsol[2]);
	if (r < -EPS_degenerate) return 0;
	if (r < 0.0) r = 0.0; /* clamp small negative to zero */

	if (out_center) { out_center[0] = xc; out_center[1] = yc; }
	if (out_radius) *out_radius = r;
	return 1;
}


