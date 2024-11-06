/*
updateHy.cu

update Hy
*/

#include "orcwa.h"
#include "orcwa_cuda.h"
#include "finc_cuda.h"


__host__ __device__
static void updateHy_f_vector(
	int i, int j, int k,
	real_t hy[], const real_t ez[], const real_t ex[],
	const real_t k1[], const real_t k2[],
	real_t rzc, real_t rxc, param_t *p)
{
	const int64_t n = LA(p, i, j, k);

	hy[n] = k1[n] * hy[n]
	      - k2[n] * (rzc * (ex[n + p->Nk] - ex[n])
	               - rxc * (ez[n + p->Ni] - ez[n]));
}


__host__ __device__
static void updateHy_f_no_vector(
	int i, int j, int k,
	real_t hy[], const real_t ez[], const real_t ex[], const id_t ihy[],
	const real_t d1[], const real_t d2[],
	real_t rzc, real_t rxc, param_t *p)
{
	const int64_t n = LA(p, i, j, k);
	const int64_t m = ihy[n];

	hy[n] = d1[m] * hy[n]
	      - d2[m] * (rzc * (ex[n + p->Nk] - ex[n])
	               - rxc * (ez[n + p->Ni] - ez[n]));
}


__host__ __device__
static void updateHy_p_vector(
	int i, int j, int k,
	real_t hy[], const real_t ez[], const real_t ex[], const id_t ihy[],
	const real_t k1[], const real_t k2[],
	real_t rzc, real_t rxc, param_t *p,
	real_t x, real_t y, real_t z, real_t t)
{
	const int64_t n = LA(p, i, j, k);
	const int64_t m = ihy[n];

	if (m == 0) {
		hy[n] -= rzc * (ex[n + p->Nk] - ex[n])
		       - rxc * (ez[n + p->Ni] - ez[n]);
	}
	else {
		real_t fi, dfi;
		finc_cuda(x, y, z, t, p->r0, p->ri, p->hi[1], p->ai, p->dt, &fi, &dfi);
		if (m == PEC) {
			hy[n] = -fi;
		}
		else {
			hy[n] = k1[n] * hy[n]
			      - k2[n] * (rzc * (ex[n + p->Nk] - ex[n])
			               - rxc * (ez[n + p->Ni] - ez[n]))
			      - (k1[n] - k2[n]) * dfi
			      - (1 - k1[n]) * fi;
		}
	}
}


__host__ __device__
static void updateHy_p_no_vector(
	int i, int j, int k,
	real_t hy[], const real_t ez[], const real_t ex[], const id_t ihy[],
	const real_t d1[], const real_t d2[],
	real_t rzc, real_t rxc, param_t *p,
	real_t x, real_t y, real_t z, real_t t)
{
	const int64_t n = LA(p, i, j, k);
	const int64_t m = ihy[n];

	if (m == 0) {
		hy[n] -= rzc * (ex[n + p->Nk] - ex[n])
		       - rxc * (ez[n + p->Ni] - ez[n]);
	}
	else {
		real_t fi, dfi;
		finc_cuda(x, y, z, t, p->r0, p->ri, p->hi[1], p->ai, p->dt, &fi, &dfi);
		if (m == PEC) {
			hy[n] = -fi;
		}
		else {
			hy[n] = d1[m] * hy[n]
			      - d2[m] * (rzc * (ex[n + p->Nk] - ex[n])
			               - rxc * (ez[n + p->Ni] - ez[n]))
			      - (d1[m] - d2[m]) * dfi
			      - (1 - d1[m]) * fi;
		}
	}
}


__global__
static void updateHy_gpu(int vector,
	real_t hy[], const real_t ez[], const real_t ex[], const id_t ihy[],
	const real_t d1[], const real_t d2[], const real_t k1[], const real_t k2[],
	const real_t rzc[], const real_t rxc[], const real_t xc[], const real_t yn[], const real_t zc[], real_t t)
{
	const int i = d_Param.iMin + threadIdx.z + (blockIdx.z * blockDim.z);
	const int j = d_Param.jMin + threadIdx.y + (blockIdx.y * blockDim.y);
	const int k = d_Param.kMin + threadIdx.x + (blockIdx.x * blockDim.x);
	if ((i <  d_Param.iMax) &&
	    (j <= d_Param.jMax) &&
	    (k <  d_Param.kMax)) {
		if (d_Param.NFeed) {
			if (vector) {
				updateHy_f_vector(
					i, j, k,
					hy, ez, ex,
					k1, k2,
					rzc[k], rxc[i], &d_Param);
			}
			else {
				updateHy_f_no_vector(
					i, j, k,
					hy, ez, ex, ihy,
					d1, d2,
					rzc[k], rxc[i], &d_Param);
			}
		}
		else if (d_Param.IPlanewave) {
			if (vector) {
				updateHy_p_vector(
					i, j, k,
					hy, ez, ex, ihy,
					k1, k2,
					rzc[k], rxc[i], &d_Param,
					xc[i], yn[j], zc[k], t);
			}
			else {
				updateHy_p_no_vector(
					i, j, k,
					hy, ez, ex, ihy,
					d1, d2,
					rzc[k], rxc[i], &d_Param,
					xc[i], yn[j], zc[k], t);
			}
		}
	}
}


static void updateHy_cpu(
	real_t hy[], const real_t ez[], const real_t ex[], const id_t ihy[],
	const real_t d1[], const real_t d2[],
	const real_t rzc[], const real_t rxc[], const real_t xc[], const real_t yn[], const real_t zc[], real_t t)
{
	for (int i = h_Param.iMin; i <  h_Param.iMax; i++) {
	for (int j = h_Param.jMin; j <= h_Param.jMax; j++) {
	for (int k = h_Param.kMin; k <  h_Param.kMax; k++) {
		if (h_Param.NFeed) {
			updateHy_f_no_vector(
				i, j, k,
				hy, ez, ex, ihy,
				d1, d2,
				rzc[k], rxc[i], &h_Param);
		}
		else if (h_Param.IPlanewave) {
			updateHy_p_no_vector(
				i, j, k,
				hy, ez, ex, ihy,
				d1, d2,
				rzc[k], rxc[i], &h_Param,
				xc[i], yn[j], zc[k], t);
		}
	}
	}
	}
}


void updateHy(double t)
{
	if (GPU) {
		cudaMemcpyToSymbol(d_Param, &h_Param, sizeof(param_t));
		dim3 grid(
			CEIL(kMax - kMin + 0, updateBlock.x),
			CEIL(jMax - jMin + 1, updateBlock.y),
			CEIL(iMax - iMin + 0, updateBlock.z));
		updateHy_gpu<<<grid, updateBlock>>>(VECTOR,
			Hy, Ez, Ex, d_iHy,
			d_D1, d_D2, d_K1Hy, d_K2Hy,
			d_RZc, d_RXc, d_Xc, d_Yn, d_Zc, (real_t)t);
		if (UM) cudaDeviceSynchronize();
	}
	else {
		// CPU : no-vector only
		updateHy_cpu(
			Hy, Ez, Ex, iHy,
			D1, D2,
			RZc, RXc, h_Xc, h_Yn, h_Zc, (real_t)t);
	}
}
