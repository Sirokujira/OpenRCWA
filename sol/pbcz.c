/*
pbcz.c

PBC for Z boundary
*/

#include "orcwa.h"

void pbcz(void)
{
	const int id1 = -1;
	const int id2 = 0;
	const int id3 = Nz - 1;
	const int id4 = Nz;
	int i;

	// Hx
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (    i = iMin - 0; i <= iMax; i++) {
	for (int j = jMin - 1; j <= jMax; j++) {
		HX(i, j, id1) = HX(i, j, id3);
		HX(i, j, id4) = HX(i, j, id2);
	}
	}

	// Hy
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (    i = iMin - 1; i <= iMax; i++) {
	for (int j = jMin - 0; j <= jMax; j++) {
		HY(i, j, id1) = HY(i, j, id3);
		HY(i, j, id4) = HY(i, j, id2);
	}
	}
}
