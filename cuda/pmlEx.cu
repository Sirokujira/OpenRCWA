/*
pmlEx.cu

PML ABC for Ex
*/

#include "orcwa.h"
#include "orcwa_cuda.h"

__host__ __device__
static void pmlex(
	int64_t n, int64_t nj, int64_t nk,
	real_t *ex, real_t *hy, real_t *hz, real_t *exy, real_t *exz,
	real_t ry, real_t rz, real_t gpmlyn, real_t gpmlzn, real_t rm)
{
	*exy = (*exy + (ry * (hz[n] - hz[n - nj]))) / (1 + (gpmlyn * rm));
	*exz = (*exz - (rz * (hy[n] - hy[n - nk]))) / (1 + (gpmlzn * rm));
	ex[n] = *exy + *exz;
}

__global__
static void pmlex_gpu(
	int ny, int nz,
	int64_t ni, int64_t nj, int64_t nk, int64_t n0,
	real_t *ex, real_t *hy, real_t *hz, real_t *exy, real_t *exz,
	int l, int64_t numpmlex,
	pml_t *fpmlex, real_t *rpmle, real_t *rpml, real_t *ryn, real_t *rzn, real_t *gpmlyn, real_t *gpmlzn)
{
	int64_t n = threadIdx.x + (blockIdx.x * blockDim.x);
	if (n < numpmlex) {
		const int  i = fpmlex[n].i;
		const int  j = fpmlex[n].j;
		const int  k = fpmlex[n].k;
		const id_t m = fpmlex[n].m;
		const real_t ry = ryn[MIN(MAX(j, 0), ny    )] * rpmle[m];
		const real_t rz = rzn[MIN(MAX(k, 0), nz    )] * rpmle[m];
		const int64_t nc = (ni * i) + (nj * j) + (nk * k) + n0;
		pmlex(
			nc, nj, nk,
			ex, hy, hz, &exy[n], &exz[n],
			ry, rz, gpmlyn[j + l], gpmlzn[k + l], rpml[m]);
	}
}

static void pmlex_cpu(
	int ny, int nz,
	int64_t ni, int64_t nj, int64_t nk, int64_t n0,
	real_t *ex, real_t *hy, real_t *hz, real_t *exy, real_t *exz,
	int l, int64_t numpmlex,
	pml_t *fpmlex, real_t *rpmle, real_t *rpml, real_t *ryn, real_t *rzn, real_t *gpmlyn, real_t *gpmlzn)
{
	for (int64_t n = 0; n < numpmlex; n++) {
		const int  i = fpmlex[n].i;
		const int  j = fpmlex[n].j;
		const int  k = fpmlex[n].k;
		const id_t m = fpmlex[n].m;
		const real_t ry = ryn[MIN(MAX(j, 0), ny    )] * rpmle[m];
		const real_t rz = rzn[MIN(MAX(k, 0), nz    )] * rpmle[m];
		const int64_t nc = (ni * i) + (nj * j) + (nk * k) + n0;
		pmlex(
			nc, nj, nk,
			ex, hy, hz, &exy[n], &exz[n],
			ry, rz, gpmlyn[j + l], gpmlzn[k + l], rpml[m]);
	}
}

void pmlEx()
{
	if (GPU) {
		pmlex_gpu<<<(int)CEIL(numPmlEx, pmlBlock), pmlBlock>>>(
			Ny, Nz,
			Ni, Nj, Nk, N0,
			Ex, Hy, Hz, Exy, Exz,
			cPML.l, numPmlEx,
			d_fPmlEx, d_rPmlE, d_rPml, d_RYn, d_RZn, d_gPmlYn, d_gPmlZn);
		if (UM) cudaDeviceSynchronize();
	}
	else {
		pmlex_cpu(
			Ny, Nz,
			Ni, Nj, Nk, N0,
			Ex, Hy, Hz, Exy, Exz,
			cPML.l, numPmlEx,
			fPmlEx, rPmlE, rPml, RYn, RZn, gPmlYn, gPmlZn);
	}
}
