/*
updateHx.c

update Hx
*/

#include "orcwa.h"
#include "finc.h"

static void updateHx_f_vector(void)
{
	assert(Nk == 1);

	int i;
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (    i = iMin; i <= iMax; i++) {
	for (int j = jMin; j <  jMax; j++) {
		int64_t n = NA(i, j, kMin);
		for (int k = kMin; k <  kMax; k++) {
			Hx[n] = K1Hx[n] * Hx[n]
			      - K2Hx[n] * (RYc[j] * (Ez[n + Nj] - Ez[n])
			                 - RZc[k] * (Ey[n + Nk] - Ey[n]));
			n++;
		}
	}
	}
}


static void updateHx_f_no_vector(void)
{
	assert(Nk == 1);

	int i;
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (    i = iMin; i <= iMax; i++) {
	for (int j = jMin; j <  jMax; j++) {
		int64_t n = NA(i, j, kMin);
		for (int k = kMin; k <  kMax; k++) {
			const int64_t m = iHx[n];
			Hx[n] = D1[m] * Hx[n]
			      - D2[m] * (RYc[j] * (Ez[n + Nj] - Ez[n])
			               - RZc[k] * (Ey[n + Nk] - Ey[n]));
			n++;
		}
	}
	}
}


static void updateHx_p_vector(double t)
{
	assert(Nk == 1);

	int i;
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (    i = iMin; i <= iMax; i++) {
	for (int j = jMin; j <  jMax; j++) {
		int64_t n = NA(i, j, kMin);
		for (int k = kMin; k <  kMax; k++) {
			real_t fi, dfi;
			finc(Xn[i], Yc[j], Zc[k], t, Planewave.r0, Planewave.ri, Planewave.hi[0], Planewave.ai, Dt, &fi, &dfi);
			Hx[n] = K1Hx[n] * Hx[n]
			      - K2Hx[n] * (RYc[j] * (Ez[n + Nj] - Ez[n])
			                 - RZc[k] * (Ey[n + Nk] - Ey[n]))
			      - (K1Hx[n] - K2Hx[n]) * dfi
			      - (1 - K1Hx[n]) * fi;
			n++;
		}
	}
	}
}


static void updateHx_p_no_vector(double t)
{
	assert(Nk == 1);

	int i;
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (    i = iMin; i <= iMax; i++) {
	for (int j = jMin; j <  jMax; j++) {
		int64_t n = NA(i, j, kMin);
		for (int k = kMin; k <  kMax; k++) {
			const int64_t m = iHx[n];
			if (m == 0) {
				Hx[n] -= RYc[j] * (Ez[n + Nj] - Ez[n])
				       - RZc[k] * (Ey[n + Nk] - Ey[n]);
			}
			else {
				real_t fi, dfi;
				finc(Xn[i], Yc[j], Zc[k], t, Planewave.r0, Planewave.ri, Planewave.hi[0], Planewave.ai, Dt, &fi, &dfi);
				if (m == PEC) {
					Hx[n] = -fi;
				}
				else {
					Hx[n] = D1[m] * Hx[n]
					      - D2[m] * (RYc[j] * (Ez[n + Nj] - Ez[n])
					               - RZc[k] * (Ey[n + Nk] - Ey[n]))
					      - (D1[m] - D2[m]) * dfi
					      - (1 - D1[m]) * fi;
				}
			}
			n++;
		}
	}
	}
}


void updateHx(double t)
{
	if (NFeed) {
		if (VECTOR) {
			updateHx_f_vector();
		}
		else {
			updateHx_f_no_vector();
		}
	}
	else if (IPlanewave) {
		if (VECTOR) {
			updateHx_p_vector(t);
		}
		else {
			updateHx_p_no_vector(t);
		}
	}
}
