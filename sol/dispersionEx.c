/*
dispersionEx.c

update Ex (dispersion)
*/

#include "orcwa.h"
#include "finc.h"


void dispersionEx(double t)
{
	int64_t n;
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (n = 0; n < numDispersionEx; n++) {
		const int     i = mDispersionEx[n].i;
		const int     j = mDispersionEx[n].j;
		const int     k = mDispersionEx[n].k;
		const real_t f1 = mDispersionEx[n].f1;
		const real_t f2 = mDispersionEx[n].f2;
		const real_t f3 = mDispersionEx[n].f3;

		real_t fi = 0;
		if (IPlanewave) {
			real_t dfi;
			finc(Xc[i], Yn[j], Zn[k], t, Planewave.r0, Planewave.ri, Planewave.ei[0], Planewave.ai, Dt, &fi, &dfi);
		}

		EX(i, j, k) += f1 * DispersionEx[n];

		DispersionEx[n] = f2 * (EX(i, j, k) + fi)
		                + f3 * DispersionEx[n];
	}
}
