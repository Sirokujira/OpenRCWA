/*
dispersionEz.c

update Ez (dispersion)
*/

#include "orcwa.h"
#include "finc.h"


void dispersionEz(double t)
{
	int64_t n;
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (n = 0; n < numDispersionEz; n++) {
		const int     i = mDispersionEz[n].i;
		const int     j = mDispersionEz[n].j;
		const int     k = mDispersionEz[n].k;
		const real_t f1 = mDispersionEz[n].f1;
		const real_t f2 = mDispersionEz[n].f2;
		const real_t f3 = mDispersionEz[n].f3;

		real_t fi = 0;
		if (IPlanewave) {
			real_t dfi;
			finc(Xn[i], Yn[j], Zc[k], t, Planewave.r0, Planewave.ri, Planewave.ei[2], Planewave.ai, Dt, &fi, &dfi);
		}

		EZ(i, j, k) += f1 * DispersionEz[n];

		DispersionEz[n] = f2 * (EZ(i, j, k) + fi)
		                + f3 * DispersionEz[n];
	}
}
