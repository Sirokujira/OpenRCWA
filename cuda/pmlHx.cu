/*
pmlHx.cu

PML ABC for Hx
*/

#include "orcwa.h"
#include "orcwa_cuda.h"

__host__ __device__
static void pmlhx(
	int64_t n, int64_t nj, int64_t nk,
	real_t *hx, real_t *ey, real_t *ez, real_t *hxy, real_t *hxz,
	real_t ry, real_t rz, real_t gpmlyc, real_t gpmlzc, real_t rm)
{
	*hxy = (*hxy - (ry * (ez[n + nj] - ez[n]))) / (1 + (gpmlyc * rm));
	*hxz = (*hxz + (rz * (ey[n + nk] - ey[n]))) / (1 + (gpmlzc * rm));
	hx[n] = *hxy + *hxz;
}

__global__
static void pmlhx_gpu(
	int ny, int nz,
	int64_t ni, int64_t nj, int64_t nk, int64_t n0,
	real_t *hx, real_t *ey, real_t *ez, real_t *hxy, real_t *hxz,
	int l, int64_t numpmlhx,
	pml_t *fpmlhx, real_t *rpmlh, real_t *rpml, real_t *ryc, real_t *rzc, real_t *gpmlyc, real_t *gpmlzc)
{
	int64_t n = threadIdx.x + (blockIdx.x * blockDim.x);
	if (n < numpmlhx) {
		const int  i = fpmlhx[n].i;
		const int  j = fpmlhx[n].j;
		const int  k = fpmlhx[n].k;
		const id_t m = fpmlhx[n].m;
		const real_t ry = ryc[MIN(MAX(j, 0), ny - 1)] * rpmlh[m];
		const real_t rz = rzc[MIN(MAX(k, 0), nz - 1)] * rpmlh[m];
		const int64_t nc = (ni * i) + (nj * j) + (nk * k) + n0;
		pmlhx(
			nc, nj, nk,
			hx, ey, ez, &hxy[n], &hxz[n],
			ry, rz, gpmlyc[j + l], gpmlzc[k + l], rpml[m]);
	}
}

static void pmlhx_cpu(
	int ny, int nz,
	int64_t ni, int64_t nj, int64_t nk, int64_t n0,
	real_t *hx, real_t *ey, real_t *ez, real_t *hxy, real_t *hxz,
	int l, int64_t numpmlhx,
	pml_t *fpmlhx, real_t *rpmlh, real_t *rpml, real_t *ryc, real_t *rzc, real_t *gpmlyc, real_t *gpmlzc)
{
	for (int64_t n = 0; n < numpmlhx; n++) {
		const int  i = fpmlhx[n].i;
		const int  j = fpmlhx[n].j;
		const int  k = fpmlhx[n].k;
		const id_t m = fpmlhx[n].m;
		const real_t ry = ryc[MIN(MAX(j, 0), ny - 1)] * rpmlh[m];
		const real_t rz = rzc[MIN(MAX(k, 0), nz - 1)] * rpmlh[m];
		const int64_t nc = (ni * i) + (nj * j) + (nk * k) + n0;
		pmlhx(
			nc, nj, nk,
			hx, ey, ez, &hxy[n], &hxz[n],
			ry, rz, gpmlyc[j + l], gpmlzc[k + l], rpml[m]);
	}
}

void pmlHx()
{
	if (GPU) {
		pmlhx_gpu<<<(int)CEIL(numPmlHx, pmlBlock), pmlBlock>>>(
			Ny, Nz,
			Ni, Nj, Nk, N0,
			Hx, Ey, Ez, Hxy, Hxz,
			cPML.l, numPmlHx,
			d_fPmlHx, d_rPmlH, d_rPml, d_RYc, d_RZc, d_gPmlYc, d_gPmlZc);
		if (UM) cudaDeviceSynchronize();
	}
	else {
		pmlhx_cpu(
			Ny, Nz,
			Ni, Nj, Nk, N0,
			Hx, Ey, Ez, Hxy, Hxz,
			cPML.l, numPmlHx,
			fPmlHx, rPmlH, rPml, RYc, RZc, gPmlYc, gPmlZc);
	}
}
