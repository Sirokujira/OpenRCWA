/*
updateEy.c

update Ey
*/

#include "orcwa.h"
#include "finc.h"

static void updateEy_f_vector(void)
{
	assert(Nk == 1);

	int i;
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (    i = iMin; i <= iMax; i++) {
	for (int j = jMin; j <  jMax; j++) {
		int64_t n = NA(i, j, kMin);
		for (int k = kMin; k <= kMax; k++) {
			Ey[n] = K1Ey[n] * Ey[n]
			      + K2Ey[n] * (RZn[k] * (Hx[n] - Hx[n - Nk])
			                 - RXn[i] * (Hz[n] - Hz[n - Ni]));
			n++;
		}
	}
	}
}


static void updateEy_f_no_vector(void)
{
	assert(Nk == 1);

	int i;
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (    i = iMin; i <= iMax; i++) {
	for (int j = jMin; j <  jMax; j++) {
		int64_t n = NA(i, j, kMin);
		for (int k = kMin; k <= kMax; k++) {
			const int64_t m = iEy[n];
			Ey[n] = C1[m] * Ey[n]
			      + C2[m] * (RZn[k] * (Hx[n] - Hx[n - Nk])
			               - RXn[i] * (Hz[n] - Hz[n - Ni]));
			n++;
		}
	}
	}
}


static void updateEy_p_vector(double t)
{
	assert(Nk == 1);

	int i;
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (    i = iMin; i <= iMax; i++) {
	for (int j = jMin; j <  jMax; j++) {
		int64_t n = NA(i, j, kMin);
		for (int k = kMin; k <= kMax; k++) {
			real_t fi, dfi;
			finc(Xn[i], Yc[j], Zn[k], t, Planewave.r0, Planewave.ri, Planewave.ei[1], Planewave.ai, Dt, &fi, &dfi);
			Ey[n] = K1Ey[n] * Ey[n]
			      + K2Ey[n] * (RZn[k] * (Hx[n] - Hx[n - Nk])
			                 - RXn[i] * (Hz[n] - Hz[n - Ni]))
			      - (K1Ey[n] - K2Ey[n]) * dfi
			      - (1 - K1Ey[n]) * fi;
			n++;
		}
	}
	}
}


static void updateEy_p_no_vector(double t)
{
	assert(Nk == 1);

	int i;
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (    i = iMin; i <= iMax; i++) {
	for (int j = jMin; j <  jMax; j++) {
		int64_t n = NA(i, j, kMin);
		for (int k = kMin; k <= kMax; k++) {
			const int64_t m = iEy[n];
			if (m == 0) {
				Ey[n] += RZn[k] * (Hx[n] - Hx[n - Nk])
				       - RXn[i] * (Hz[n] - Hz[n - Ni]);
			}
			else {
				real_t fi, dfi;
				finc(Xn[i], Yc[j], Zn[k], t, Planewave.r0, Planewave.ri, Planewave.ei[1], Planewave.ai, Dt, &fi, &dfi);
				if (m == PEC) {
					Ey[n] = -fi;
				}
				else {
					Ey[n] = C1[m] * Ey[n]
					      + C2[m] * (RZn[k] * (Hx[n] - Hx[n - Nk])
					               - RXn[i] * (Hz[n] - Hz[n - Ni]))
					      - (C1[m] - C2[m]) * dfi
					      - (1 - C1[m]) * fi;
				}
			}
			n++;
		}
	}
	}
}


void updateEy(double t)
{
	if (NFeed) {
		if (VECTOR) {
			updateEy_f_vector();
		}
		else {
			updateEy_f_no_vector();
		}
	}
	else if (IPlanewave) {
		if (VECTOR) {
			updateEy_p_vector(t);
		}
		else {
			updateEy_p_no_vector(t);
		}
	}
}
