/*
average.c (OpenMP)

E/H average
*/

#include "orcwa.h"

void average(double fsum[])
{
	double se = 0;
	double sh = 0;
	int    i;
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : se, sh)
#endif
	for (    i = iMin; i < iMax; i++) {
	for (int j = jMin; j < jMax; j++) {
	for (int k = kMin; k < kMax; k++) {
		se +=
			+ fabs(
				+ EX(i    , j    , k    )
				+ EX(i    , j + 1, k    )
				+ EX(i    , j    , k + 1)
				+ EX(i    , j + 1, k + 1))
			+ fabs(
				+ EY(i    , j    , k    )
				+ EY(i    , j    , k + 1)
				+ EY(i + 1, j    , k    )
				+ EY(i + 1, j    , k + 1))
			+ fabs(
				+ EZ(i    , j    , k    )
				+ EZ(i + 1, j    , k    )
				+ EZ(i    , j + 1, k    )
				+ EZ(i + 1, j + 1, k    ));
		sh +=
			+ fabs(
				+ HX(i    , j    , k    )
				+ HX(i + 1, j    , k    ))
			+ fabs(
				+ HY(i    , j    , k    )
				+ HY(i    , j + 1, k    ))
			+ fabs(
				+ HZ(i    , j    , k    )
				+ HZ(i    , j    , k + 1));
	}
	}
	}

	fsum[0] = se / (4.0 * Nx * Ny * Nz);
	fsum[1] = sh / (2.0 * Nx * Ny * Nz);
}
