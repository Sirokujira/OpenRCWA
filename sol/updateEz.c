/*
updateEz.c

update Ez
*/

#include "orcwa.h"
#include "finc.h"

static void updateEz_f_vector(void)
{
	assert(Nk == 1);

	int i;
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (    i = iMin; i <= iMax; i++) {
	for (int j = jMin; j <= jMax; j++) {
		int64_t n = NA(i, j, kMin);
		for (int k = kMin; k <  kMax; k++) {
			Ez[n] = K1Ez[n] * Ez[n]
			      + K2Ez[n] * (RXn[i] * (Hy[n] - Hy[n - Ni])
			                 - RYn[j] * (Hx[n] - Hx[n - Nj]));
			n++;
		}
	}
	}
}


static void updateEz_f_no_vector(void)
{
	assert(Nk == 1);

	int i;
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (    i = iMin; i <= iMax; i++) {
	for (int j = jMin; j <= jMax; j++) {
		int64_t n = NA(i, j, kMin);
		for (int k = kMin; k <  kMax; k++) {
			const int64_t m = iEz[n];
			Ez[n] = C1[m] * Ez[n]
			      + C2[m] * (RXn[i] * (Hy[n] - Hy[n - Ni])
			               - RYn[j] * (Hx[n] - Hx[n - Nj]));
			n++;
		}
	}
	}
}


static void updateEz_p_vector(double t)
{
	assert(Nk == 1);

	int i;
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (    i = iMin; i <= iMax; i++) {
	for (int j = jMin; j <= jMax; j++) {
		int64_t n = NA(i, j, kMin);
		for (int k = kMin; k <  kMax; k++) {
			real_t fi, dfi;
			finc(Xn[i], Yn[j], Zc[k], t, Planewave.r0, Planewave.ri, Planewave.ei[2], Planewave.ai, Dt, &fi, &dfi);
			Ez[n] = K1Ez[n] * Ez[n]
			      + K2Ez[n] * (RXn[i] * (Hy[n] - Hy[n - Ni])
			                 - RYn[j] * (Hx[n] - Hx[n - Nj]))
			      - (K1Ez[n] - K2Ez[n]) * dfi
			      - (1 - K1Ez[n]) * fi;
			n++;
		}
	}
	}
}


static void updateEz_p_no_vector(double t)
{
	assert(Nk == 1);

	int i;
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (    i = iMin; i <= iMax; i++) {
	for (int j = jMin; j <= jMax; j++) {
		int64_t n = NA(i, j, kMin);
		for (int k = kMin; k <  kMax; k++) {
			const int m = iEz[n];
			if (m == 0) {
				Ez[n] += RXn[i] * (Hy[n] - Hy[n - Ni])
				       - RYn[j] * (Hx[n] - Hx[n - Nj]);
			}
			else {
				real_t fi, dfi;
				finc(Xn[i], Yn[j], Zc[k], t, Planewave.r0, Planewave.ri, Planewave.ei[2], Planewave.ai, Dt, &fi, &dfi);
				if (m == PEC) {
					Ez[n] = -fi;
				}
				else {
					Ez[n] = C1[m] * Ez[n]
					      + C2[m] * (RXn[i] * (Hy[n] - Hy[n - Ni])
					               - RYn[j] * (Hx[n] - Hx[n - Nj]))
					      - (C1[m] - C2[m]) * dfi
					      - (1 - C1[m]) * fi;
				}
			}
			n++;
		}
	}
	}
}


void updateEz(double t)
{
	if (NFeed) {
		if (VECTOR) {
			updateEz_f_vector();
		}
		else {
			updateEz_f_no_vector();
		}
	}
	else if (IPlanewave) {
		if (VECTOR) {
			updateEz_p_vector(t);
		}
		else {
			updateEz_p_no_vector(t);
		}
	}
}
