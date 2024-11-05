/*
updateEz.cu

update Ez
*/

#include "orcwa.h"
#include "orcwa_cuda.h"
#include "finc_cuda.h"


__host__ __device__
static void updateEz_f_vector(
	int i, int j, int k,
	real_t ez[], const real_t hx[], const real_t hy[],
	const real_t k1[], const real_t k2[],
	real_t rxn, real_t ryn, param_t *p)
{
	const int64_t n = LA(p, i, j, k);

	ez[n] = k1[n] * ez[n]
	      + k2[n] * (rxn * (hy[n] - hy[n - p->Ni])
	               - ryn * (hx[n] - hx[n - p->Nj]));
}


__host__ __device__
static void updateEz_f_no_vector(
	int i, int j, int k,
	real_t ez[], const real_t hx[], const real_t hy[], const id_t iez[],
	const real_t c1[], const real_t c2[],
	real_t rxn, real_t ryn, param_t *p)
{
	const int64_t n = LA(p, i, j, k);
	const int64_t m = iez[n];

	ez[n] = c1[m] * ez[n]
	      + c2[m] * (rxn * (hy[n] - hy[n - p->Ni])
	               - ryn * (hx[n] - hx[n - p->Nj]));
}


__host__ __device__
static void updateEz_p_vector(
	int i, int j, int k,
	real_t ez[], const real_t hx[], const real_t hy[], const id_t iez[],
	const real_t k1[], const real_t k2[],
	real_t rxn, real_t ryn, param_t *p,
	real_t x, real_t y, real_t z, real_t t)
{
	const int64_t n = LA(p, i, j, k);
	const int64_t m = iez[n];

	if (m == 0) {
		ez[n] += rxn * (hy[n] - hy[n - p->Ni])
		       - ryn * (hx[n] - hx[n - p->Nj]);
	}
	else {
		real_t fi, dfi;
		finc_cuda(x, y, z, t, p->r0, p->ri, p->ei[2], p->ai, p->dt, &fi, &dfi);
		if (m == PEC) {
			ez[n] = -fi;
		}
		else {
			ez[n] = k1[n] * ez[n]
			      + k2[n] * (rxn * (hy[n] - hy[n - p->Ni])
			               - ryn * (hx[n] - hx[n - p->Nj]))
			      - (k1[n] - k2[n]) * dfi
			      - (1 - k1[n]) * fi;
		}
	}
}


__host__ __device__
static void updateEz_p_no_vector(
	int i, int j, int k,
	real_t ez[], const real_t hx[], const real_t hy[], const id_t iez[],
	const real_t c1[], const real_t c2[],
	real_t rxn, real_t ryn, param_t *p,
	real_t x, real_t y, real_t z, real_t t)
{
	const int64_t n = LA(p, i, j, k);
	const int64_t m = iez[n];

	if (m == 0) {
		ez[n] += rxn * (hy[n] - hy[n - p->Ni])
		       - ryn * (hx[n] - hx[n - p->Nj]);
	}
	else {
		real_t fi, dfi;
		finc_cuda(x, y, z, t, p->r0, p->ri, p->ei[2], p->ai, p->dt, &fi, &dfi);
		if (m == PEC) {
			ez[n] = -fi;
		}
		else {
			ez[n] = c1[m] * ez[n]
			      + c2[m] * (rxn * (hy[n] - hy[n - p->Ni])
			               - ryn * (hx[n] - hx[n - p->Nj]))
			      - (c1[m] - c2[m]) * dfi
			      - (1 - c1[m]) * fi;
		}
	}
}


__global__
static void updateEz_gpu(int vector,
	real_t ez[], const real_t hx[], const real_t hy[], const id_t iez[],
	const real_t c1[], const real_t c2[], const real_t k1[], const real_t k2[],
	const real_t rxn[], const real_t ryn[], const real_t xn[], const real_t yn[], const real_t zc[], real_t t)
{
	const int i = d_Param.iMin + threadIdx.z + (blockIdx.z * blockDim.z);
	const int j = d_Param.jMin + threadIdx.y + (blockIdx.y * blockDim.y);
	const int k = d_Param.kMin + threadIdx.x + (blockIdx.x * blockDim.x);
	if ((i <= d_Param.iMax) &&
	    (j <= d_Param.jMax) &&
	    (k <  d_Param.kMax)) {
		if (d_Param.NFeed) {
			if (vector) {
				updateEz_f_vector(
					i, j, k,
					ez, hx, hy,
					k1, k2,
					rxn[i], ryn[j], &d_Param);
			}
			else {
				updateEz_f_no_vector(
					i, j, k,
					ez, hx, hy, iez,
					c1, c2,
					rxn[i], ryn[j], &d_Param);
			}
		}
		else if (d_Param.IPlanewave) {
			if (vector) {
				updateEz_p_vector(
					i, j, k,
					ez, hx, hy, iez,
					k1, k2,
					rxn[i], ryn[j], &d_Param,
					xn[i], yn[j], zc[k], t);
			}
			else {
				updateEz_p_no_vector(
					i, j, k,
					ez, hx, hy, iez,
					c1, c2,
					rxn[i], ryn[j], &d_Param,
					xn[i], yn[j], zc[k], t);
			}
		}
	}
}


static void updateEz_cpu(
	real_t ez[], const real_t hx[], const real_t hy[], const id_t iez[],
	const real_t c1[], const real_t c2[],
	const real_t rxn[], const real_t ryn[], const real_t xn[], const real_t yn[], const real_t zc[], real_t t)
{
	for (int i = h_Param.iMin; i <= h_Param.iMax; i++) {
	for (int j = h_Param.jMin; j <= h_Param.jMax; j++) {
	for (int k = h_Param.kMin; k <  h_Param.kMax; k++) {
		if (h_Param.NFeed) {
			updateEz_f_no_vector(
				i, j, k,
				ez, hx, hy, iez,
				c1, c2,
				rxn[i], ryn[j], &h_Param);
		}
		else if (h_Param.IPlanewave) {
			updateEz_p_no_vector(
				i, j, k,
				ez, hx, hy, iez,
				c1, c2,
				rxn[i], ryn[j], &h_Param,
				xn[i], yn[j], zc[k], t);
		}
	}
	}
	}
}


void updateEz(double t)
{
	if (GPU) {
		cudaMemcpyToSymbol(d_Param, &h_Param, sizeof(param_t));
		dim3 grid(
			CEIL(kMax - kMin + 0, updateBlock.x),
			CEIL(jMax - jMin + 1, updateBlock.y),
			CEIL(iMax - iMin + 1, updateBlock.z));
		updateEz_gpu<<<grid, updateBlock>>>(VECTOR,
			Ez, Hx, Hy, d_iEz,
			d_C1, d_C2, d_K1Ez, d_K2Ez,
			d_RXn, d_RYn, d_Xn, d_Yn, d_Zc, (real_t)t);
		if (UM) cudaDeviceSynchronize();
	}
	else {
		// CPU : no-vector only
		updateEz_cpu(
			Ez, Hx, Hy, iEz,
			C1, C2,
			RXn, RYn, h_Xn, h_Yn, h_Zc, (real_t)t);
	}
}
