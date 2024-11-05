/*
pbcy.c

PBC for Y boundary
*/

#include "orcwa.h"

void pbcy(void)
{
	const int id1 = -1;
	const int id2 = 0;
	const int id3 = Ny - 1;
	const int id4 = Ny;
	int k;

	// Hz
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (    k = kMin - 0; k <= kMax; k++) {
	for (int i = iMin - 1; i <= iMax; i++) {
		HZ(i, id1, k) = HZ(i, id3, k);
		HZ(i, id4, k) = HZ(i, id2, k);
	}
	}

	// Hx
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (    k = kMin - 1; k <= kMax; k++) {
	for (int i = iMin - 0; i <= iMax; i++) {
		HX(i, id1, k) = HX(i, id3, k);
		HX(i, id4, k) = HX(i, id2, k);
	}
	}
}
