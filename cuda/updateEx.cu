/*
updateEx.cu

update Ex
*/

#include "orcwa.h"
#include "orcwa_cuda.h"
#include "finc_cuda.h"


__host__ __device__
static void updateEx_f_vector(
	int i, int j, int k,
	real_t ex[], const real_t hy[], const real_t hz[],
	const real_t k1[], const real_t k2[],
	real_t ryn, real_t rzn, param_t *p)
{
	const int64_t n = LA(p, i, j, k);

	ex[n] = k1[n] * ex[n]
	      + k2[n] * (ryn * (hz[n] - hz[n - p->Nj])
	               - rzn * (hy[n] - hy[n - p->Nk]));
}


__host__ __device__
static void updateEx_f_no_vector(
	int i, int j, int k,
	real_t ex[], const real_t hy[], const real_t hz[], const id_t iex[],
	const real_t c1[], const real_t c2[],
	real_t ryn, real_t rzn, param_t *p)
{
	const int64_t n = LA(p, i, j, k);
	const int64_t m = iex[n];

	ex[n] = c1[m] * ex[n]
	      + c2[m] * (ryn * (hz[n] - hz[n - p->Nj])
	               - rzn * (hy[n] - hy[n - p->Nk]));
}


__host__ __device__
static void updateEx_p_vector(
	int i, int j, int k,
	real_t ex[], const real_t hy[], const real_t hz[], const id_t iex[],
	const real_t k1[], const real_t k2[],
	real_t ryn, real_t rzn, param_t *p,
	real_t x, real_t y, real_t z, real_t t)
{
	const int64_t n = LA(p, i, j, k);
	const int64_t m = iex[n];

	if (m == 0) {
		ex[n] += ryn * (hz[n] - hz[n - p->Nj])
		       - rzn * (hy[n] - hy[n - p->Nk]);
	}
	else {
		real_t fi, dfi;
		finc_cuda(x, y, z, t, p->r0, p->ri, p->ei[0], p->ai, p->dt, &fi, &dfi);
		if (m == PEC) {
			ex[n] = -fi;
		}
		else {
			ex[n] = k1[n] * ex[n]
			      + k2[n] * (ryn * (hz[n] - hz[n - p->Nj])
			               - rzn * (hy[n] - hy[n - p->Nk]))
			      - (k1[n] - k2[n]) * dfi
			      - (1 - k1[n]) * fi;
		}
	}
}


__host__ __device__
static void updateEx_p_no_vector(
	int i, int j, int k,
	real_t ex[], const real_t hy[], const real_t hz[], const id_t iex[],
	const real_t c1[], const real_t c2[],
	real_t ryn, real_t rzn, param_t *p,
	real_t x, real_t y, real_t z, real_t t)
{
	const int64_t n = LA(p, i, j, k);
	const int64_t m = iex[n];

	if (m == 0) {
		ex[n] += ryn * (hz[n] - hz[n - p->Nj])
		       - rzn * (hy[n] - hy[n - p->Nk]);
	}
	else {
		real_t fi, dfi;
		finc_cuda(x, y, z, t, p->r0, p->ri, p->ei[0], p->ai, p->dt, &fi, &dfi);
		if (m == PEC) {
			ex[n] = -fi;
		}
		else {
			ex[n] = c1[m] * ex[n]
			      + c2[m] * (ryn * (hz[n] - hz[n - p->Nj])
			               - rzn * (hy[n] - hy[n - p->Nk]))
			      - (c1[m] - c2[m]) * dfi
			      - (1 - c1[m]) * fi;
		}
	}
}


__global__
static void updateEx_gpu(int vector,
	real_t ex[], const real_t hy[], const real_t hz[], const id_t iex[],
	const real_t c1[], const real_t c2[], const real_t k1[], const real_t k2[],
	const real_t ryn[], const real_t rzn[], const real_t xc[], const real_t yn[], const real_t zn[], real_t t)
{
	const int i = d_Param.iMin + threadIdx.z + (blockIdx.z * blockDim.z);
	const int j = d_Param.jMin + threadIdx.y + (blockIdx.y * blockDim.y);
	const int k = d_Param.kMin + threadIdx.x + (blockIdx.x * blockDim.x);
	if ((i <  d_Param.iMax) &&
	    (j <= d_Param.jMax) &&
	    (k <= d_Param.kMax)) {
		if (d_Param.NFeed) {
			if (vector) {
				updateEx_f_vector(
					i, j, k,
					ex, hy, hz,
					k1, k2,
					ryn[j], rzn[k], &d_Param);
			}
			else {
				updateEx_f_no_vector(
					i, j, k,
					ex, hy, hz, iex,
					c1, c2,
					ryn[j], rzn[k], &d_Param);
			}
		}
		else if (d_Param.IPlanewave) {
			if (vector) {
				updateEx_p_vector(
					i, j, k,
					ex, hy, hz, iex,
					k1, k2,
					ryn[j], rzn[k], &d_Param,
					xc[i], yn[j], zn[k], t);
			}
			else {
				updateEx_p_no_vector(
					i, j, k,
					ex, hy, hz, iex,
					c1, c2,
					ryn[j], rzn[k], &d_Param,
					xc[i], yn[j], zn[k], t);
			}
		}
	}
}


static void updateEx_cpu(
	real_t ex[], const real_t hy[], const real_t hz[], const id_t iex[],
	const real_t c1[], const real_t c2[],
	const real_t ryn[], const real_t rzn[], const real_t xc[], const real_t yn[], const real_t zn[], real_t t)
{
	for (int i = h_Param.iMin; i <  h_Param.iMax; i++) {
	for (int j = h_Param.jMin; j <= h_Param.jMax; j++) {
	for (int k = h_Param.kMin; k <= h_Param.kMax; k++) {
		if (h_Param.NFeed) {
			updateEx_f_no_vector(
				i, j, k,
				ex, hy, hz, iex,
				c1, c2,
				ryn[j], rzn[k], &h_Param);
		}
		else if (h_Param.IPlanewave) {
			updateEx_p_no_vector(
				i, j, k,
				ex, hy, hz, iex,
				c1, c2,
				ryn[j], rzn[k], &h_Param,
				xc[i], yn[j], zn[k], t);
		}
	}
	}
	}
}


void updateEx(double t)
{
	if (GPU) {
		cudaMemcpyToSymbol(d_Param, &h_Param, sizeof(param_t));
		dim3 grid(
			CEIL(kMax - kMin + 1, updateBlock.x),
			CEIL(jMax - jMin + 1, updateBlock.y),
			CEIL(iMax - iMin + 0, updateBlock.z));
		updateEx_gpu<<<grid, updateBlock>>>(VECTOR,
			Ex, Hy, Hz, d_iEx,
			d_C1, d_C2, d_K1Ex, d_K2Ex,
			d_RYn, d_RZn, d_Xc, d_Yn, d_Zn, (real_t)t);
		if (UM) cudaDeviceSynchronize();
	}
	else {
		// CPU : no-vector only
		updateEx_cpu(
			Ex, Hy, Hz, iEx,
			C1, C2,
			RYn, RZn, h_Xc, h_Yn, h_Zn, (real_t)t);
	}
}
