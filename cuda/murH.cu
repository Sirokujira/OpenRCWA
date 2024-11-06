/*
murH_gpu.cu
*/

#include "orcwa.h"
#include "orcwa_cuda.h"


__host__ __device__
static void murh(real_t *h, mur_t *q, int64_t ni, int64_t nj, int64_t nk, int64_t n0)
{
	const int64_t m0 = (ni * q->i)  + (nj * q->j)  + (nk * q->k)  + n0;
	const int64_t m1 = (ni * q->i1) + (nj * q->j1) + (nk * q->k1) + n0;

	h[m0] = q->f + q->g * (h[m1] - h[m0]);
	q->f = h[m1];
}


__global__
static void murh_gpu(
	int64_t num, real_t *h, mur_t *fmur,
	int64_t ni, int64_t nj, int64_t nk, int64_t n0)
{
	const int64_t n = threadIdx.x + (blockIdx.x * blockDim.x);
	if (n < num) {
		murh(h, &fmur[n], ni, nj, nk, n0);
	}
}


static void murh_cpu(
	int64_t num, real_t *h, mur_t *fmur,
	int64_t ni, int64_t nj, int64_t nk, int64_t n0)
{
	for (int64_t n = 0; n < num; n++) {
		murh(h, &fmur[n], ni, nj, nk, n0);
	}
}


void murH(int64_t num, mur_t *fmur, real_t *h)
{
	if (num <= 0) return;

	if (GPU) {
		murh_gpu<<<(int)CEIL(num, murBlock), murBlock>>>(
			num, h, fmur,
			Ni, Nj, Nk, N0);
		if (UM) cudaDeviceSynchronize();
	}
	else {
		murh_cpu(
			num, h, fmur,
			Ni, Nj, Nk, N0);
	}
}
