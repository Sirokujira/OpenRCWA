/*
dispersionEy.c

update Ey (dispersion)
*/

#include "orcwa.h"
#include "finc.h"


void dispersionEy(double t)
{
	int64_t n;
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (n = 0; n < numDispersionEy; n++) {
		const int     i = mDispersionEy[n].i;
		const int     j = mDispersionEy[n].j;
		const int     k = mDispersionEy[n].k;
		const real_t f1 = mDispersionEy[n].f1;
		const real_t f2 = mDispersionEy[n].f2;
		const real_t f3 = mDispersionEy[n].f3;

		real_t fi = 0;
		if (IPlanewave) {
			real_t dfi;
			finc(Xn[i], Yc[j], Zn[k], t, Planewave.r0, Planewave.ri, Planewave.ei[1], Planewave.ai, Dt, &fi, &dfi);
		}

		EY(i, j, k) += f1 * DispersionEy[n];

		DispersionEy[n] = f2 * (EY(i, j, k) + fi)
		                + f3 * DispersionEy[n];
	}
}
