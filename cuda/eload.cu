/*
eload.cu

E on loads (inductor)
*/

#include "orcwa.h"
#include "orcwa_cuda.h"

__host__ __device__
static void eload_(
	inductor_t *ptr, double cdt,
	real_t *ex, real_t *ey, real_t *ez, real_t *hx, real_t *hy, real_t *hz,
	param_t *p)
{
	const int    i  = ptr->i;
	const int    j  = ptr->j;
	const int    k  = ptr->k;
	const double dx = ptr->dx;
	const double dy = ptr->dy;
	const double dz = ptr->dz;

	if      ((ptr->dir == 'X') &&
	         (p->iMin <= i) && (i <  p->iMax) &&
	         (p->jMin <= j) && (j <= p->jMax) &&
	         (p->kMin <= k) && (k <= p->kMax)) {  // MPI
		const double roth = (hz[LA(p, i, j, k)] - hz[LA(p, i,     j - 1, k    )]) / dy
		                  - (hy[LA(p, i, j, k)] - hy[LA(p, i,     j,     k - 1)]) / dz;
		ex[LA(p, i, j, k)] = (real_t)(ptr->e + (cdt * roth) - (ptr->fctr * cdt * cdt * ptr->esum));
		ptr->e = ex[LA(p, i, j, k)];
		ptr->esum += ptr->e;
	}
	else if ((ptr->dir == 'Y') &&
	         (p->iMin <= i) && (i <= p->iMax) &&
	         (p->jMin <= j) && (j <  p->jMax) &&
	         (p->kMin <= k) && (k <= p->kMax)) {  // MPI
		const double roth = (hx[LA(p, i, j, k)] - hx[LA(p, i,     j,     k - 1)]) / dz
		                  - (hz[LA(p, i, j, k)] - hz[LA(p, i - 1, j,     k    )]) / dx;
		ey[LA(p, i, j, k)] = (real_t)(ptr->e + (cdt * roth) - (ptr->fctr * cdt * cdt * ptr->esum));
		ptr->e = ey[LA(p, i, j, k)];
		ptr->esum += ptr->e;
	}
	else if ((ptr->dir == 'Z') &&
	         (p->iMin <= i) && (i <= p->iMax) &&
	         (p->jMin <= j) && (j <= p->jMax) &&
	         (p->kMin <= k) && (k <  p->kMax)) {  // MPI
		const double roth = (hy[LA(p, i, j, k)] - hy[LA(p, i - 1, j,     k    )]) / dx
		                  - (hx[LA(p, i, j, k)] - hx[LA(p, i,     j - 1, k    )]) / dy;
		ez[LA(p, i, j, k)] = (real_t)(ptr->e + (cdt * roth) - (ptr->fctr * cdt * cdt * ptr->esum));
		ptr->e = ez[LA(p, i, j, k)];
		ptr->esum += ptr->e;
	}
}

// gpu
__global__
static void eload_gpu(
	int ninductor, inductor_t *inductor, double cdt,
	real_t *ex, real_t *ey, real_t *ez, real_t *hx, real_t *hy, real_t *hz)
{
	const int n = threadIdx.x + (blockIdx.x * blockDim.x);
	if (n < ninductor) {
		eload_(
			&inductor[n], cdt,
			ex, ey, ez, hx, hy, hz,
			&d_Param);
	}
}

// cpu
static void eload_cpu(
	int ninductor, inductor_t *inductor, double cdt,
	real_t *ex, real_t *ey, real_t *ez, real_t *hx, real_t *hy, real_t *hz)
{
	for (int n = 0; n < ninductor; n++) {
		eload_(
			&inductor[n], cdt,
			ex, ey, ez, hx, hy, hz,
			&h_Param);
	}
}

void eload(void)
{
	if (NInductor <= 0) return;

	if (GPU) {
		cudaMemcpyToSymbol(d_Param, &h_Param, sizeof(param_t));
		const int block = MIN(256, NInductor);
		const int grid = CEIL(NInductor, block);
		eload_gpu<<<grid, block>>>(
			NInductor, d_Inductor, C * Dt,
			Ex, Ey, Ez, Hx, Hy, Hz);
		if (UM) cudaDeviceSynchronize();
	}
	else {
		eload_cpu(
			NInductor, Inductor, C * Dt,
			Ex, Ey, Ez, Hx, Hy, Hz);
	}
}
