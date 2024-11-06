/*
pmlEz.cu

PML ABC for Ez
*/

#include "orcwa.h"
#include "orcwa_cuda.h"

__host__ __device__
static void pmlez(
	int64_t n, int64_t ni, int64_t nj,
	real_t *ez, real_t *hx, real_t *hy, real_t *ezx, real_t *ezy,
	real_t rx, real_t ry, real_t gpmlxn, real_t gpmlyn, real_t rm)
{
	*ezx = (*ezx + (rx * (hy[n] - hy[n - ni]))) / (1 + (gpmlxn * rm));
	*ezy = (*ezy - (ry * (hx[n] - hx[n - nj]))) / (1 + (gpmlyn * rm));
	ez[n] = *ezx + *ezy;
}

__global__
static void pmlez_gpu(
	int nx, int ny,
	int64_t ni, int64_t nj, int64_t nk, int64_t n0,
	real_t *ez, real_t *hx, real_t *hy, real_t *ezx, real_t *ezy,
	int l, int64_t numpmlez,
	pml_t *fpmlez, real_t *rpmle, real_t *rpml, real_t *rxn, real_t *ryn, real_t *gpmlxn, real_t *gpmlyn)
{
	int64_t n = threadIdx.x + (blockIdx.x * blockDim.x);
	if (n < numpmlez) {
		const int  i = fpmlez[n].i;
		const int  j = fpmlez[n].j;
		const int  k = fpmlez[n].k;
		const id_t m = fpmlez[n].m;
		const real_t rx = rxn[MIN(MAX(i, 0), nx    )] * rpmle[m];
		const real_t ry = ryn[MIN(MAX(j, 0), ny    )] * rpmle[m];
		const int64_t nc = (ni * i) + (nj * j) + (nk * k) + n0;
		pmlez(
			nc, ni, nj,
			ez, hx, hy, &ezx[n], &ezy[n],
			rx, ry, gpmlxn[i + l], gpmlyn[j + l], rpml[m]);
	}
}

static void pmlez_cpu(
	int nx, int ny,
	int64_t ni, int64_t nj, int64_t nk, int64_t n0,
	real_t *ez, real_t *hx, real_t *hy, real_t *ezx, real_t *ezy,
	int l, int64_t numpmlez,
	pml_t *fpmlez, real_t *rpmle, real_t *rpml, real_t *rxn, real_t *ryn, real_t *gpmlxn, real_t *gpmlyn)
{
	for (int64_t n = 0; n < numpmlez; n++) {
		const int  i = fpmlez[n].i;
		const int  j = fpmlez[n].j;
		const int  k = fpmlez[n].k;
		const id_t m = fpmlez[n].m;
		const real_t rx = rxn[MIN(MAX(i, 0), nx    )] * rpmle[m];
		const real_t ry = ryn[MIN(MAX(j, 0), ny    )] * rpmle[m];
		const int64_t nc = (ni * i) + (nj * j) + (nk * k) + n0;
		pmlez(
			nc, ni, nj,
			ez, hx, hy, &ezx[n], &ezy[n],
			rx, ry, gpmlxn[i + l], gpmlyn[j + l], rpml[m]);
	}
}

void pmlEz()
{
	if (GPU) {
		pmlez_gpu<<<(int)CEIL(numPmlEz, pmlBlock), pmlBlock>>>(
			Nx, Ny,
			Ni, Nj, Nk, N0,
			Ez, Hx, Hy, Ezx, Ezy,
			cPML.l, numPmlEz,
			d_fPmlEz, d_rPmlE, d_rPml, d_RXn, d_RYn, d_gPmlXn, d_gPmlYn);
		if (UM) cudaDeviceSynchronize();
	}
	else {
		pmlez_cpu(
			Nx, Ny,
			Ni, Nj, Nk, N0,
			Ez, Hx, Hy, Ezx, Ezy,
			cPML.l, numPmlEz,
			fPmlEz, rPmlE, rPml, RXn, RYn, gPmlXn, gPmlYn);
	}
}
