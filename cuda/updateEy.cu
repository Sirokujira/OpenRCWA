/*
updateEy.cu

update Ey
*/

#include "orcwa.h"
#include "orcwa_cuda.h"
#include "finc_cuda.h"


__host__ __device__
static void updateEy_f_vector(
	int i, int j, int k,
	real_t ey[], const real_t hz[], const real_t hx[],
	const real_t k1[], const real_t k2[],
	real_t rzn, real_t rxn, param_t *p)
{
	const int64_t n = LA(p, i, j, k);

	ey[n] = k1[n] * ey[n]
	      + k2[n] * (rzn * (hx[n] - hx[n - p->Nk])
	               - rxn * (hz[n] - hz[n - p->Ni]));
}


__host__ __device__
static void updateEy_f_no_vector(
	int i, int j, int k,
	real_t ey[], const real_t hz[], const real_t hx[], const id_t iey[],
	const real_t c1[], const real_t c2[],
	real_t rzn, real_t rxn, param_t *p)
{
	const int64_t n = LA(p, i, j, k);
	const int64_t m = iey[n];

	ey[n] = c1[m] * ey[n]
	      + c2[m] * (rzn * (hx[n] - hx[n - p->Nk])
	               - rxn * (hz[n] - hz[n - p->Ni]));
}


__host__ __device__
static void updateEy_p_vector(
	int i, int j, int k,
	real_t ey[], const real_t hz[], const real_t hx[], const id_t iey[],
	const real_t k1[], const real_t k2[],
	real_t rzn, real_t rxn, param_t *p,
	real_t x, real_t y, real_t z, real_t t)
{
	const int64_t n = LA(p, i, j, k);
	const int64_t m = iey[n];

	if (m == 0) {
		ey[n] += rzn * (hx[n] - hx[n - p->Nk])
		       - rxn * (hz[n] - hz[n - p->Ni]);
	}
	else {
		real_t fi, dfi;
		finc_cuda(x, y, z, t, p->r0, p->ri, p->ei[1], p->ai, p->dt, &fi, &dfi);
		if (m == PEC) {
			ey[n] = -fi;
		}
		else {
			ey[n] = k1[n] * ey[n]
			      + k2[n] * (rzn * (hx[n] - hx[n - p->Nk])
			               - rxn * (hz[n] - hz[n - p->Ni]))
			      - (k1[n] - k2[n]) * dfi
			      - (1 - k1[n]) * fi;
		}
	}
}


__host__ __device__
static void updateEy_p_no_vector(
	int i, int j, int k,
	real_t ey[], const real_t hz[], const real_t hx[], const id_t iey[],
	const real_t c1[], const real_t c2[],
	real_t rzn, real_t rxn, param_t *p,
	real_t x, real_t y, real_t z, real_t t)
{
	const int64_t n = LA(p, i, j, k);
	const int64_t m = iey[n];

	if (m == 0) {
		ey[n] += rzn * (hx[n] - hx[n - p->Nk])
		       - rxn * (hz[n] - hz[n - p->Ni]);
	}
	else {
		real_t fi, dfi;
		finc_cuda(x, y, z, t, p->r0, p->ri, p->ei[1], p->ai, p->dt, &fi, &dfi);
		if (m == PEC) {
			ey[n] = -fi;
		}
		else {
			ey[n] = c1[m] * ey[n]
			      + c2[m] * (rzn * (hx[n] - hx[n - p->Nk])
			               - rxn * (hz[n] - hz[n - p->Ni]))
			      - (c1[m] - c2[m]) * dfi
			      - (1 - c1[m]) * fi;
		}
	}
}


__global__
static void updateEy_gpu(int vector,
	real_t ey[], const real_t hz[], const real_t hx[], const id_t iey[],
	const real_t c1[], const real_t c2[], const real_t k1[], const real_t k2[],
	const real_t rzn[], const real_t rxn[], const real_t xn[], const real_t yc[], const real_t zn[], real_t t)
{
	const int i = d_Param.iMin + threadIdx.z + (blockIdx.z * blockDim.z);
	const int j = d_Param.jMin + threadIdx.y + (blockIdx.y * blockDim.y);
	const int k = d_Param.kMin + threadIdx.x + (blockIdx.x * blockDim.x);
	if ((i <= d_Param.iMax) &&
	    (j <  d_Param.jMax) &&
	    (k <= d_Param.kMax)) {
		if (d_Param.NFeed) {
			if (vector) {
				updateEy_f_vector(
					i, j, k,
					ey, hz, hx,
					k1, k2,
					rzn[k], rxn[i], &d_Param);
			}
			else {
				updateEy_f_no_vector(
					i, j, k,
					ey, hz, hx, iey,
					c1, c2,
					rzn[k], rxn[i], &d_Param);
			}
		}
		else if (d_Param.IPlanewave) {
			if (vector) {
				updateEy_p_vector(
					i, j, k,
					ey, hz, hx, iey,
					k1, k2,
					rzn[k], rxn[i], &d_Param,
					xn[i], yc[j], zn[k], t);
			}
			else {
				updateEy_p_no_vector(
					i, j, k,
					ey, hz, hx, iey,
					c1, c2,
					rzn[k], rxn[i], &d_Param,
					xn[i], yc[j], zn[k], t);
			}
		}
	}
}


static void updateEy_cpu(
	real_t ey[], const real_t hz[], const real_t hx[], const id_t iey[],
	const real_t c1[], const real_t c2[],
	const real_t rzn[], const real_t rxn[], const real_t xn[], const real_t yc[], const real_t zn[], real_t t)
{
	for (int i = h_Param.iMin; i <= h_Param.iMax; i++) {
	for (int j = h_Param.jMin; j <  h_Param.jMax; j++) {
	for (int k = h_Param.kMin; k <= h_Param.kMax; k++) {
		if (h_Param.NFeed) {
			updateEy_f_no_vector(
				i, j, k,
				ey, hz, hx, iey,
				c1, c2,
				rzn[k], rxn[i], &h_Param);
		}
		else if (h_Param.IPlanewave) {
			updateEy_p_no_vector(
				i, j, k,
				ey, hz, hx, iey,
				c1, c2,
				rzn[k], rxn[i], &h_Param,
				xn[i], yc[j], zn[k], t);
		}
	}
	}
	}
}


void updateEy(double t)
{
	if (GPU) {
		cudaMemcpyToSymbol(d_Param, &h_Param, sizeof(param_t));
		dim3 grid(
			CEIL(kMax - kMin + 1, updateBlock.x),
			CEIL(jMax - jMin + 0, updateBlock.y),
			CEIL(iMax - iMin + 1, updateBlock.z));
		updateEy_gpu<<<grid, updateBlock>>>(VECTOR,
			Ey, Hz, Hx, d_iEy,
			d_C1, d_C2, d_K1Ey, d_K2Ey,
			d_RZn, d_RXn, d_Xn, d_Yc, d_Zn, (real_t)t);
		if (UM) cudaDeviceSynchronize();
	}
	else {
		// CPU : no-vector only
		updateEy_cpu(
			Ey, Hz, Hx, iEy,
			C1, C2,
			RZn, RXn, h_Xn, h_Yc, h_Zn, (real_t)t);
	}
}
