/*
murH.c (OpenMP)

Mur for Hx/Hy/Hz
*/

#include "orcwa.h"

void murH(int64_t num, mur_t *fmur, real_t *h)
{
	int64_t n;
#ifdef _OPENMP
#pragma omp parallel for
#endif
#ifdef __NEC__
#pragma _NEC ivdep
#endif
#ifdef __CLANG_FUJITSU
#pragma clang loop vectorize(assume_safety)
#endif
	for (n = 0; n < num; n++) {
		const int i  = fmur[n].i;
		const int j  = fmur[n].j;
		const int k  = fmur[n].k;
		const int i1 = fmur[n].i1;
		const int j1 = fmur[n].j1;
		const int k1 = fmur[n].k1;
		h[NA(i, j, k)] = fmur[n].f
		               + fmur[n].g * (h[NA(i1, j1, k1)] - h[NA(i, j, k)]);
		fmur[n].f = h[NA(i1, j1, k1)];
	}
}
