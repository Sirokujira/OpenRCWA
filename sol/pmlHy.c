/*
pmlHy.c (OpenMP)

PML for Hy
*/

#include "orcwa.h"

void pmlHy(void)
{
	const int lz = cPML.l;
	const int lx = cPML.l;

	int64_t n;
#ifdef _OPENMP
#pragma omp parallel for
#endif
#ifdef __CLANG_FUJITSU
#pragma clang loop vectorize(assume_safety)
#endif
	for (n = 0; n < numPmlHy; n++) {
		const int  i = fPmlHy[n].i;
		const int  j = fPmlHy[n].j;
		const int  k = fPmlHy[n].k;
		const id_t m = fPmlHy[n].m;

		const real_t dex = EX(i, j, k + 1) - EX(i, j, k);
		const real_t rz = RZc[MIN(MAX(k, 0), Nz - 1)] * rPmlH[m];
		Hyz[n] = (Hyz[n] - (rz * dex)) / (1 + (gPmlZc[k + lz] * rPml[m]));

		const real_t dez = EZ(i + 1, j, k) - EZ(i, j, k);
		const real_t rx = RXc[MIN(MAX(i, 0), Nx - 1)] * rPmlH[m];
		Hyx[n] = (Hyx[n] + (rx * dez)) / (1 + (gPmlXc[i + lx] * rPml[m]));

		HY(i, j, k) = Hyz[n] + Hyx[n];
	}
}
