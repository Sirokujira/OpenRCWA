/*
pmlHy.cu

PML ABC for Hy
*/

#include "orcwa.h"
#include "orcwa_cuda.h"

__host__ __device__
static void pmlhy(
	int64_t n, int64_t nk, int64_t ni,
	real_t *hy, real_t *ez, real_t *ex, real_t *hyz, real_t *hyx,
	real_t rz, real_t rx, real_t gpmlzc, real_t gpmlxc, real_t rm)
{
	*hyz = (*hyz - (rz * (ex[n + nk] - ex[n]))) / (1 + (gpmlzc * rm));
	*hyx = (*hyx + (rx * (ez[n + ni] - ez[n]))) / (1 + (gpmlxc * rm));
	hy[n] = *hyz + *hyx;
}

__global__
static void pmlhy_gpu(
	int nz, int nx,
	int64_t ni, int64_t nj, int64_t nk, int64_t n0,
	real_t *hy, real_t *ez, real_t *ex, real_t *hyz, real_t *hyx,
	int l, int64_t numpmlhy,
	pml_t *fpmlhy, real_t *rpmlh, real_t *rpml, real_t *rzc, real_t *rxc, real_t *gpmlzc, real_t *gpmlxc)
{
	int64_t n = threadIdx.x + (blockIdx.x * blockDim.x);
	if (n < numpmlhy) {
		const int  i = fpmlhy[n].i;
		const int  j = fpmlhy[n].j;
		const int  k = fpmlhy[n].k;
		const id_t m = fpmlhy[n].m;
		const real_t rz = rzc[MIN(MAX(k, 0), nz - 1)] * rpmlh[m];
		const real_t rx = rxc[MIN(MAX(i, 0), nx - 1)] * rpmlh[m];
		const int64_t nc = (ni * i) + (nj * j) + (nk * k) + n0;
		pmlhy(
			nc, nk, ni,
			hy, ez, ex, &hyz[n], &hyx[n],
			rz, rx, gpmlzc[k + l], gpmlxc[i + l], rpml[m]);
	}
}

static void pmlhy_cpu(
	int nz, int nx,
	int64_t ni, int64_t nj, int64_t nk, int64_t n0,
	real_t *hy, real_t *ez, real_t *ex, real_t *hyz, real_t *hyx,
	int l, int64_t numpmlhy,
	pml_t *fpmlhy, real_t *rpmlh, real_t *rpml, real_t *rzc, real_t *rxc, real_t *gpmlzc, real_t *gpmlxc)
{
	for (int64_t n = 0; n < numpmlhy; n++) {
		const int  i = fpmlhy[n].i;
		const int  j = fpmlhy[n].j;
		const int  k = fpmlhy[n].k;
		const id_t m = fpmlhy[n].m;
		const real_t rz = rzc[MIN(MAX(k, 0), nz - 1)] * rpmlh[m];
		const real_t rx = rxc[MIN(MAX(i, 0), nx - 1)] * rpmlh[m];
		const int64_t nc = (ni * i) + (nj * j) + (nk * k) + n0;
		pmlhy(
			nc, nk, ni,
			hy, ez, ex, &hyz[n], &hyx[n],
			rz, rx, gpmlzc[k + l], gpmlxc[i + l], rpml[m]);
	}
}

void pmlHy()
{
	if (GPU) {
		pmlhy_gpu<<<(int)CEIL(numPmlHy, pmlBlock), pmlBlock>>>(
			Nz, Nx,
			Ni, Nj, Nk, N0,
			Hy, Ez, Ex, Hyz, Hyx,
			cPML.l, numPmlHy,
			d_fPmlHy, d_rPmlH, d_rPml, d_RZc, d_RXc, d_gPmlZc, d_gPmlXc);
		if (UM) cudaDeviceSynchronize();
	}
	else {
		pmlhy_cpu(
			Nz, Nx,
			Ni, Nj, Nk, N0,
			Hy, Ez, Ex, Hyz, Hyx,
			cPML.l, numPmlHy,
			fPmlHy, rPmlH, rPml, RZc, RXc, gPmlZc, gPmlXc);
	}
}
