/*
pmlHz.c (OpenMP)

PML for Hz
*/

#include "orcwa.h"

void pmlHz(void)
{
	const int lx = cPML.l;
	const int ly = cPML.l;

	int64_t n;
#ifdef _OPENMP
#pragma omp parallel for
#endif
#ifdef __CLANG_FUJITSU
#pragma clang loop vectorize(assume_safety)
#endif
	for (n = 0; n < numPmlHz; n++) {
		const int  i = fPmlHz[n].i;
		const int  j = fPmlHz[n].j;
		const int  k = fPmlHz[n].k;
		const id_t m = fPmlHz[n].m;

		const real_t dey = EY(i + 1, j, k) - EY(i, j, k);
		const real_t rx = RXc[MIN(MAX(i, 0), Nx - 1)] * rPmlH[m];
		Hzx[n] = (Hzx[n] - (rx * dey)) / (1 + (gPmlXc[i + lx] * rPml[m]));

		const real_t dex = EX(i, j + 1, k) - EX(i, j, k);
		const real_t ry = RYc[MIN(MAX(j, 0), Ny - 1)] * rPmlH[m];
		Hzy[n] = (Hzy[n] + (ry * dex)) / (1 + (gPmlYc[j + ly] * rPml[m]));

		HZ(i, j, k) = Hzx[n] + Hzy[n];
	}
}
