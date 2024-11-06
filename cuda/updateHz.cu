/*
updateHz.cu

update Hz
*/

#include "orcwa.h"
#include "orcwa_cuda.h"
#include "finc_cuda.h"


__host__ __device__
static void updateHz_f_vector(
	int i, int j, int k,
	real_t hz[], const real_t ex[], const real_t ey[],
	const real_t k1[], const real_t k2[],
	real_t rxc, real_t ryc, param_t *p)
{
	const int64_t n = LA(p, i, j, k);

	hz[n] = k1[n] * hz[n]
	      - k2[n] * (rxc * (ey[n + p->Ni] - ey[n])
	               - ryc * (ex[n + p->Nj] - ex[n]));
}


__host__ __device__
static void updateHz_f_no_vector(
	int i, int j, int k,
	real_t hz[], const real_t ex[], const real_t ey[], const id_t ihz[],
	const real_t d1[], const real_t d2[],
	real_t rxc, real_t ryc, param_t *p)
{
	const int64_t n = LA(p, i, j, k);
	const int64_t m = ihz[n];

	hz[n] = d1[m] * hz[n]
	      - d2[m] * (rxc * (ey[n + p->Ni] - ey[n])
	               - ryc * (ex[n + p->Nj] - ex[n]));
}


__host__ __device__
static void updateHz_p_vector(
	int i, int j, int k,
	real_t hz[], const real_t ex[], const real_t ey[], const id_t ihz[],
	const real_t k1[], const real_t k2[],
	real_t rxc, real_t ryc, param_t *p,
	real_t x, real_t y, real_t z, real_t t)
{
	const int64_t n = LA(p, i, j, k);
	const int64_t m = ihz[n];

	if (m == 0) {
		hz[n] -= rxc * (ey[n + p->Ni] - ey[n])
		       - ryc * (ex[n + p->Nj] - ex[n]);
	}
	else {
		real_t fi, dfi;
		finc_cuda(x, y, z, t, p->r0, p->ri, p->hi[2], p->ai, p->dt, &fi, &dfi);
		if (m == PEC) {
			hz[n] = -fi;
		}
		else {
			hz[n] = k1[n] * hz[n]
			      - k2[n] * (rxc * (ey[n + p->Ni] - ey[n])
			               - ryc * (ex[n + p->Nj] - ex[n]))
			      - (k1[n] - k2[n]) * dfi
			      - (1 - k1[n]) * fi;
		}
	}
}


__host__ __device__
static void updateHz_p_no_vector(
	int i, int j, int k,
	real_t hz[], const real_t ex[], const real_t ey[], const id_t ihz[],
	const real_t d1[], const real_t d2[],
	real_t rxc, real_t ryc, param_t *p,
	real_t x, real_t y, real_t z, real_t t)
{
	const int64_t n = LA(p, i, j, k);
	const int64_t m = ihz[n];

	if (m == 0) {
		hz[n] -= rxc * (ey[n + p->Ni] - ey[n])
		       - ryc * (ex[n + p->Nj] - ex[n]);
	}
	else {
		real_t fi, dfi;
		finc_cuda(x, y, z, t, p->r0, p->ri, p->hi[2], p->ai, p->dt, &fi, &dfi);
		if (m == PEC) {
			hz[n] = -fi;
		}
		else {
			hz[n] = d1[m] * hz[n]
			      - d2[m] * (rxc * (ey[n + p->Ni] - ey[n])
			               - ryc * (ex[n + p->Nj] - ex[n]))
			      - (d1[m] - d2[m]) * dfi
			      - (1 - d1[m]) * fi;
		}
	}
}


__global__
static void updateHz_gpu(int vector,
	real_t hz[], const real_t ex[], const real_t ey[], const id_t ihz[],
	const real_t d1[], const real_t d2[], const real_t k1[], const real_t k2[],
	const real_t rxc[], const real_t ryc[], const real_t xc[], const real_t yc[], const real_t zn[], real_t t)
{
	const int i = d_Param.iMin + threadIdx.z + (blockIdx.z * blockDim.z);
	const int j = d_Param.jMin + threadIdx.y + (blockIdx.y * blockDim.y);
	const int k = d_Param.kMin + threadIdx.x + (blockIdx.x * blockDim.x);
	if ((i <  d_Param.iMax) &&
	    (j <  d_Param.jMax) &&
	    (k <= d_Param.kMax)) {
		if (d_Param.NFeed) {
			if (vector) {
				updateHz_f_vector(
					i, j, k,
					hz, ex, ey,
					k1, k2,
					rxc[i], ryc[j], &d_Param);
			}
			else {
				updateHz_f_no_vector(
					i, j, k,
					hz, ex, ey, ihz,
					d1, d2,
					rxc[i], ryc[j], &d_Param);
			}
		}
		else if (d_Param.IPlanewave) {
			if (vector) {
				updateHz_p_vector(
					i, j, k,
					hz, ex, ey, ihz,
					k1, k2,
					rxc[i], ryc[j], &d_Param,
					xc[i], yc[j], zn[k], t);
			}
			else {
				updateHz_p_no_vector(
					i, j, k,
					hz, ex, ey, ihz,
					d1, d2,
					rxc[i], ryc[j], &d_Param,
					xc[i], yc[j], zn[k], t);
			}
		}
	}
}


static void updateHz_cpu(
	real_t hz[], const real_t ex[], const real_t ey[], const id_t ihz[],
	const real_t d1[], const real_t d2[],
	const real_t rxc[], const real_t ryc[], const real_t xc[], const real_t yc[], const real_t zn[], real_t t)
{
	for (int i = h_Param.iMin; i <  h_Param.iMax; i++) {
	for (int j = h_Param.jMin; j <  h_Param.jMax; j++) {
	for (int k = h_Param.kMin; k <= h_Param.kMax; k++) {
		if (h_Param.NFeed) {
			updateHz_f_no_vector(
				i, j, k,
				hz, ex, ey, ihz,
				d1, d2,
				rxc[i], ryc[j], &h_Param);
		}
		else if (h_Param.IPlanewave) {
			updateHz_p_no_vector(
				i, j, k,
				hz, ex, ey, ihz,
				d1, d2,
				rxc[i], ryc[j], &h_Param,
				xc[i], yc[j], zn[k], t);
		}
	}
	}
	}
}


void updateHz(double t)
{
	if (GPU) {
		cudaMemcpyToSymbol(d_Param, &h_Param, sizeof(param_t));
		dim3 grid(
			CEIL(kMax - kMin + 1, updateBlock.x),
			CEIL(jMax - jMin + 0, updateBlock.y),
			CEIL(iMax - iMin + 0, updateBlock.z));
		updateHz_gpu<<<grid, updateBlock>>>(VECTOR,
			Hz, Ex, Ey, d_iHz,
			d_D1, d_D2, d_K1Hz, d_K2Hz,
			d_RXc, d_RYc, d_Xc, d_Yc, d_Zn, (real_t)t);
		if (UM) cudaDeviceSynchronize();
	}
	else {
		// CPU : no-vector only
		updateHz_cpu(
			Hz, Ex, Ey, iHz,
			D1, D2,
			RXc, RYc, h_Xc, h_Yc, h_Zn, (real_t)t);
	}
}
