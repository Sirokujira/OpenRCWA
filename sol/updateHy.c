/*
updateHy.c

update Hy
*/

#include "orcwa.h"
#include "finc.h"

static void updateHy_f_vector(void)
{
	assert(Nk == 1);

	int i;
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (    i = iMin; i <  iMax; i++) {
	for (int j = jMin; j <= jMax; j++) {
		int64_t n = NA(i, j, kMin);
		for (int k = kMin; k <  kMax; k++) {
			Hy[n] = K1Hy[n] * Hy[n]
			      - K2Hy[n] * (RZc[k] * (Ex[n + Nk] - Ex[n])
			                 - RXc[i] * (Ez[n + Ni] - Ez[n]));
			n++;
		}
	}
	}
}


static void updateHy_f_no_vector(void)
{
	assert(Nk == 1);

	int i;
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (    i = iMin; i <  iMax; i++) {
	for (int j = jMin; j <= jMax; j++) {
		int64_t n = NA(i, j, kMin);
		for (int k = kMin; k <  kMax; k++) {
			const int64_t m = iHy[n];
			Hy[n] = D1[m] * Hy[n]
			      - D2[m] * (RZc[k] * (Ex[n + Nk] - Ex[n])
			               - RXc[i] * (Ez[n + Ni] - Ez[n]));
			n++;
		}
	}
	}
}


static void updateHy_p_vector(double t)
{
	assert(Nk == 1);

	int i;
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (    i = iMin; i <  iMax; i++) {
	for (int j = jMin; j <= jMax; j++) {
		int64_t n = NA(i, j, kMin);
		for (int k = kMin; k <  kMax; k++) {
			real_t fi, dfi;
			finc(Xc[i], Yn[j], Zc[k], t, Planewave.r0, Planewave.ri, Planewave.hi[1], Planewave.ai, Dt, &fi, &dfi);
			Hy[n] = K1Hy[n] * Hy[n]
			      - K2Hy[n] * (RZc[k] * (Ex[n + Nk] - Ex[n])
			                 - RXc[i] * (Ez[n + Ni] - Ez[n]))
			      - (K1Hy[n] - K2Hy[n]) * dfi
			      - (1 - K1Hy[n]) * fi;
			n++;
		}
	}
	}
}


static void updateHy_p_no_vector(double t)
{
	assert(Nk == 1);

	int i;
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (    i = iMin; i <  iMax; i++) {
	for (int j = jMin; j <= jMax; j++) {
		int64_t n = NA(i, j, kMin);
		for (int k = kMin; k <  kMax; k++) {
			const int64_t m = iHy[n];
			if (m == 0) {
				Hy[n] -= RZc[k] * (Ex[n + Nk] - Ex[n])
				       - RXc[i] * (Ez[n + Ni] - Ez[n]);
			}
			else {
				real_t fi, dfi;
				finc(Xc[i], Yn[j], Zc[k], t, Planewave.r0, Planewave.ri, Planewave.hi[1], Planewave.ai, Dt, &fi, &dfi);
				if (m == PEC) {
					Hy[n] = -fi;
				}
				else {
					Hy[n] = D1[m] * Hy[n]
					      - D2[m] * (RZc[k] * (Ex[n + Nk] - Ex[n])
					               - RXc[i] * (Ez[n + Ni] - Ez[n]))
					      - (D1[m] - D2[m]) * dfi
					      - (1 - D1[m]) * fi;
				}
			}
			n++;
		}
	}
	}
}


void updateHy(double t)
{
	if (NFeed) {
		if (VECTOR) {
			updateHy_f_vector();
		}
		else {
			updateHy_f_no_vector();
		}
	}
	else if (IPlanewave) {
		if (VECTOR) {
			updateHy_p_vector(t);
		}
		else {
			updateHy_p_no_vector(t);
		}
	}
}
