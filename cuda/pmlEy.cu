/*
pmlEy.cu

PML ABC for Ey
*/

#include "orcwa.h"
#include "orcwa_cuda.h"

__host__ __device__
static void pmley(
	int64_t n, int64_t nk, int64_t ni,
	real_t *ey, real_t *hz, real_t *hx, real_t *eyz, real_t *eyx,
	real_t rz, real_t rx, real_t gpmlzn, real_t gpmlxn, real_t rm)
{
	*eyz = (*eyz + (rz * (hx[n] - hx[n - nk]))) / (1 + (gpmlzn * rm));
	*eyx = (*eyx - (rx * (hz[n] - hz[n - ni]))) / (1 + (gpmlxn * rm));
	ey[n] = *eyz + *eyx;
}

__global__
static void pmley_gpu(
	int nz, int nx,
	int64_t ni, int64_t nj, int64_t nk, int64_t n0,
	real_t *ey, real_t *hz, real_t *hx, real_t *eyz, real_t *eyx,
	int l, int64_t numpmley,
	pml_t *fpmley, real_t *rpmle, real_t *rpml, real_t *rzn, real_t *rxn, real_t *gpmlzn, real_t *gpmlxn)
{
	int64_t n = threadIdx.x + (blockIdx.x * blockDim.x);
	if (n < numpmley) {
		const int  i = fpmley[n].i;
		const int  j = fpmley[n].j;
		const int  k = fpmley[n].k;
		const id_t m = fpmley[n].m;
		const real_t rz = rzn[MIN(MAX(k, 0), nz    )] * rpmle[m];
		const real_t rx = rxn[MIN(MAX(i, 0), nx    )] * rpmle[m];
		const int64_t nc = (ni * i) + (nj * j) + (nk * k) + n0;
		pmley(
			nc, nk, ni,
			ey, hz, hx, &eyz[n], &eyx[n],
			rz, rx, gpmlzn[k + l], gpmlxn[i + l], rpml[m]);
	}
}

static void pmley_cpu(
	int nz, int nx,
	int64_t ni, int64_t nj, int64_t nk, int64_t n0,
	real_t *ey, real_t *hz, real_t *hx, real_t *eyz, real_t *eyx,
	int l, int64_t numpmley,
	pml_t *fpmley, real_t *rpmle, real_t *rpml, real_t *rzn, real_t *rxn, real_t *gpmlzn, real_t *gpmlxn)
{
	for (int64_t n = 0; n < numpmley; n++) {
		const int  i = fpmley[n].i;
		const int  j = fpmley[n].j;
		const int  k = fpmley[n].k;
		const id_t m = fpmley[n].m;
		const real_t rz = rzn[MIN(MAX(k, 0), nz    )] * rpmle[m];
		const real_t rx = rxn[MIN(MAX(i, 0), nx    )] * rpmle[m];
		const int64_t nc = (ni * i) + (nj * j) + (nk * k) + n0;
		pmley(
			nc, nk, ni,
			ey, hz, hx, &eyz[n], &eyx[n],
			rz, rx, gpmlzn[k + l], gpmlxn[i + l], rpml[m]);
	}
}

void pmlEy()
{
	if (GPU) {
		pmley_gpu<<<(int)CEIL(numPmlEy, pmlBlock), pmlBlock>>>(
			Nz, Nx,
			Ni, Nj, Nk, N0,
			Ey, Hz, Hx, Eyz, Eyx,
			cPML.l, numPmlEy,
			d_fPmlEy, d_rPmlE, d_rPml, d_RZn, d_RXn, d_gPmlZn, d_gPmlXn);
		if (UM) cudaDeviceSynchronize();
	}
	else {
		pmley_cpu(
			Nz, Nx,
			Ni, Nj, Nk, N0,
			Ey, Hz, Hx, Eyz, Eyx,
			cPML.l, numPmlEy,
			fPmlEy, rPmlE, rPml, RZn, RXn, gPmlZn, gPmlXn);
	}
}
