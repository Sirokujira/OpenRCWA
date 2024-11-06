/*
pbcx.c

PBC for X boundary
*/

#include "orcwa.h"

void pbcx(void)
{
	const int id1 = -1;
	const int id2 = 0;
	const int id3 = Nx - 1;
	const int id4 = Nx;
	int j;

	// Hy
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (    j = jMin - 0; j <= jMax; j++) {
	for (int k = kMin - 1; k <= kMax; k++) {
		HY(id1, j, k) = HY(id3, j, k);
		HY(id4, j, k) = HY(id2, j, k);
	}
	}

	// Hz
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (    j = jMin - 1; j <= jMax; j++) {
	for (int k = kMin - 0; k <= kMax; k++) {
		HZ(id1, j, k) = HZ(id3, j, k);
		HZ(id4, j, k) = HZ(id2, j, k);
	}
	}
}
