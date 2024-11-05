/*
updateHx.cu

update Hx
*/

#include "orcwa.h"
#include "orcwa_cuda.h"
#include "finc_cuda.h"


__host__ __device__
static void updateHx_f_vector(
	int i, int j, int k,
	real_t hx[], const real_t ey[], const real_t ez[],
	const real_t k1[], const real_t k2[],
	real_t ryc, real_t rzc, param_t *p)
{
	const int64_t n = LA(p, i, j, k);

	hx[n] = k1[n] * hx[n]
	      - k2[n] * (ryc * (ez[n + p->Nj] - ez[n])
	               - rzc * (ey[n + p->Nk] - ey[n]));
}


__host__ __device__
static void updateHx_f_no_vector(
	int i, int j, int k,
	real_t hx[], const real_t ey[], const real_t ez[], const id_t ihx[],
	const real_t d1[], const real_t d2[],
	real_t ryc, real_t rzc, param_t *p)
{
	const int64_t n = LA(p, i, j, k);
	const int64_t m = ihx[n];

	hx[n] = d1[m] * hx[n]
	      - d2[m] * (ryc * (ez[n + p->Nj] - ez[n])
	               - rzc * (ey[n + p->Nk] - ey[n]));
}


__host__ __device__
static void updateHx_p_vector(
	int i, int j, int k,
	real_t hx[], const real_t ey[], const real_t ez[], const id_t ihx[],
	const real_t k1[], const real_t k2[],
	real_t ryc, real_t rzc, param_t *p,
	real_t x, real_t y, real_t z, real_t t)
{
	const int64_t n = LA(p, i, j, k);
	const int64_t m = ihx[n];

	if (m == 0) {
		hx[n] -= ryc * (ez[n + p->Nj] - ez[n])
		       - rzc * (ey[n + p->Nk] - ey[n]);
	}
	else {
		real_t fi, dfi;
		finc_cuda(x, y, z, t, p->r0, p->ri, p->hi[0], p->ai, p->dt, &fi, &dfi);
		if (m == PEC) {
			hx[n] = -fi;
		}
		else {
			hx[n] = k1[n] * hx[n]
			      - k2[n] * (ryc * (ez[n + p->Nj] - ez[n])
			               - rzc * (ey[n + p->Nk] - ey[n]))
			      - (k1[n] - k2[n]) * dfi
			      - (1 - k1[n]) * fi;
		}
	}
}


__host__ __device__
static void updateHx_p_no_vector(
	int i, int j, int k,
	real_t hx[], const real_t ey[], const real_t ez[], const id_t ihx[],
	const real_t d1[], const real_t d2[],
	real_t ryc, real_t rzc, param_t *p,
	real_t x, real_t y, real_t z, real_t t)
{
	const int64_t n = LA(p, i, j, k);
	const int64_t m = ihx[n];

	if (m == 0) {
		hx[n] -= ryc * (ez[n + p->Nj] - ez[n])
		       - rzc * (ey[n + p->Nk] - ey[n]);
	}
	else {
		real_t fi, dfi;
		finc_cuda(x, y, z, t, p->r0, p->ri, p->hi[0], p->ai, p->dt, &fi, &dfi);
		if (m == PEC) {
			hx[n] = -fi;
		}
		else {
			hx[n] = d1[m] * hx[n]
			      - d2[m] * (ryc * (ez[n + p->Nj] - ez[n])
			               - rzc * (ey[n + p->Nk] - ey[n]))
			      - (d1[m] - d2[m]) * dfi
			      - (1 - d1[m]) * fi;
		}
	}
}


__global__
static void updateHx_gpu(int vector,
	real_t hx[], const real_t ey[], const real_t ez[], const id_t ihx[],
	const real_t d1[], const real_t d2[], const real_t k1[], const real_t k2[],
	const real_t ryc[], const real_t rzc[], const real_t xn[], const real_t yc[], const real_t zc[], real_t t)
{
	const int i = d_Param.iMin + threadIdx.z + (blockIdx.z * blockDim.z);
	const int j = d_Param.jMin + threadIdx.y + (blockIdx.y * blockDim.y);
	const int k = d_Param.kMin + threadIdx.x + (blockIdx.x * blockDim.x);
	if ((i <= d_Param.iMax) &&
	    (j <  d_Param.jMax) &&
	    (k <  d_Param.kMax)) {
		if (d_Param.NFeed) {
			if (vector) {
				updateHx_f_vector(
					i, j, k,
					hx, ey, ez,
					k1, k2,
					ryc[j], rzc[k], &d_Param);
			}
			else {
				updateHx_f_no_vector(
					i, j, k,
					hx, ey, ez, ihx,
					d1, d2,
					ryc[j], rzc[k], &d_Param);
			}
		}
		else if (d_Param.IPlanewave) {
			if (vector) {
				updateHx_p_vector(
					i, j, k,
					hx, ey, ez, ihx,
					k1, k2,
					ryc[j], rzc[k], &d_Param,
					xn[i], yc[j], zc[k], t);
			}
			else {
				updateHx_p_no_vector(
					i, j, k,
					hx, ey, ez, ihx,
					d1, d2,
					ryc[j], rzc[k], &d_Param,
					xn[i], yc[j], zc[k], t);
			}
		}
	}
}


static void updateHx_cpu(
	real_t hx[], const real_t ey[], const real_t ez[], const id_t ihx[],
	const real_t d1[], const real_t d2[],
	const real_t ryc[], const real_t rzc[], const real_t xn[], const real_t yc[], const real_t zc[], real_t t)
{
	for (int i = h_Param.iMin; i <= h_Param.iMax; i++) {
	for (int j = h_Param.jMin; j <  h_Param.jMax; j++) {
	for (int k = h_Param.kMin; k <  h_Param.kMax; k++) {
		if (h_Param.NFeed) {
			updateHx_f_no_vector(
				i, j, k,
				hx, ey, ez, ihx,
				d1, d2,
				ryc[j], rzc[k], &h_Param);
		}
		else if (h_Param.IPlanewave) {
			updateHx_p_no_vector(
				i, j, k,
				hx, ey, ez, ihx,
				d1, d2,
				ryc[j], rzc[k], &h_Param,
				xn[i], yc[j], zc[k], t);
		}
	}
	}
	}
}


void updateHx(double t)
{
	if (GPU) {
		cudaMemcpyToSymbol(d_Param, &h_Param, sizeof(param_t));
		dim3 grid(
			CEIL(kMax - kMin + 0, updateBlock.x),
			CEIL(jMax - jMin + 0, updateBlock.y),
			CEIL(iMax - iMin + 1, updateBlock.z));
		updateHx_gpu<<<grid, updateBlock>>>(VECTOR,
			Hx, Ey, Ez, d_iHx,
			d_D1, d_D2, d_K1Hx, d_K2Hx,
			d_RYc, d_RZc, d_Xn, d_Yc, d_Zc, (real_t)t);
		if (UM) cudaDeviceSynchronize();
	}
	else {
		// CPU : no-vector only
		updateHx_cpu(
			Hx, Ey, Ez, iHx,
			D1, D2,
			RYc, RZc, h_Xn, h_Yc, h_Zc, (real_t)t);
	}
}
