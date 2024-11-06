/*
vpoint.cu

V waveform on points
*/

#include "orcwa.h"
#include "orcwa_cuda.h"
#include "finc_cuda.h"

__host__ __device__
static void vpoint_(
	point_t point, double *vpoint,
	real_t *xn, real_t *yn, real_t *zn, real_t *xc, real_t *yc, real_t *zc, real_t t,
	real_t *ex, real_t *ey, real_t *ez,
	param_t *p)
{
	const int i = point.i;
	const int j = point.j;
	const int k = point.k;
	double e = 0;
	double d = 0;
	real_t einc = 0;
	real_t dummy = 0;

	if      ((point.dir == 'X') &&
	         (p->iMin <= i) && (i <  p->iMax) &&
	         (p->jMin <= j) && (j <= p->jMax) &&
	         (p->kMin <= k) && (k <= p->kMax)) {  // MPI
		e = ex[LA(p, i, j, k)];
		d = point.dx;
		if (p->IPlanewave) {
			finc_cuda(xc[i], yn[j], zn[k], t, p->r0, p->ri, p->ei[0], p->ai, p->dt, &einc, &dummy);
		}
	}
	else if ((point.dir == 'Y') &&
	         (p->iMin <= i) && (i <= p->iMax) &&
	         (p->jMin <= j) && (j <  p->jMax) &&
	         (p->kMin <= k) && (k <= p->kMax)) {  // MPI
		e = ey[LA(p, i, j, k)];
		d = point.dy;
		if (p->IPlanewave) {
			finc_cuda(xn[i], yc[j], zn[k], t, p->r0, p->ri, p->ei[1], p->ai, p->dt, &einc, &dummy);
		}
	}
	else if ((point.dir == 'Z') &&
	         (p->iMin <= i) && (i <= p->iMax) &&
	         (p->jMin <= j) && (j <= p->jMax) &&
	         (p->kMin <= k) && (k <  p->kMax)) {  // MPI
		e = ez[LA(p, i, j, k)];
		d = point.dz;
		if (p->IPlanewave) {
			finc_cuda(xn[i], yn[j], zc[k], t, p->r0, p->ri, p->ei[2], p->ai, p->dt, &einc, &dummy);
		}
	}

	*vpoint = (e + einc) * (-d);
}

// gpu
__global__
static void vpoint_gpu(
	int npoint, point_t *point, double *vpoint,
	int itime, int maxiter,
	real_t *xn, real_t *yn, real_t *zn, real_t *xc, real_t *yc, real_t *zc, real_t t,
	real_t *ex, real_t *ey, real_t *ez)
{
	const int tid = threadIdx.x + (blockIdx.x * blockDim.x);
	if (tid < npoint + 2) {
		const int adr = tid * (maxiter + 1) + itime;
		vpoint_(
			point[tid], &vpoint[adr],
			xn, yn, zn, xc, yc, zc, t,
			ex, ey, ez,
			&d_Param);
	}
}

// cpu
static void vpoint_cpu(
	int npoint, point_t *point, double *vpoint,
	int itime, int maxiter,
	real_t *xn, real_t *yn, real_t *zn, real_t *xc, real_t *yc, real_t *zc, real_t t,
	real_t *ex, real_t *ey, real_t *ez)
{
	for (int n = 0; n < npoint + 2; n++) {
		const int adr = n * (maxiter + 1) + itime;
		vpoint_(
			point[n], &vpoint[adr],
			xn, yn, zn, xc, yc, zc, t,
			ex, ey, ez,
			&h_Param);
	}
}

void vpoint(int itime)
{
	if (NPoint <= 0) return;

	const real_t t = (real_t)((itime + 1) * Dt);

	if (GPU) {
		cudaMemcpyToSymbol(d_Param, &h_Param, sizeof(param_t));
		const int block = MIN(256, NPoint + 2);
		const int grid = CEIL(NPoint + 2, block);
		vpoint_gpu<<<grid, block>>>(
			NPoint, d_Point, d_VPoint,
			itime, Solver.maxiter,
			d_Xn, d_Yn, d_Zn, d_Xc, d_Yc, d_Zc, t,
			Ex, Ey, Ez);
		if (UM) cudaDeviceSynchronize();
	}
	else {
		vpoint_cpu(
			NPoint, Point, VPoint,
			itime, Solver.maxiter,
			h_Xn, h_Yn, h_Zn, h_Xc, h_Yc, h_Zc, t,
			Ex, Ey, Ez);
	}
}
