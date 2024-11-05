/*
updateHz.c

update Hz
*/

#include "orcwa.h"
#include "finc.h"

static void updateHz_f_vector(void)
{
	assert(Nk == 1);

	int i;
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (    i = iMin; i <  iMax; i++) {
	for (int j = jMin; j <  jMax; j++) {
		int64_t n = NA(i, j, kMin);
		for (int k = kMin; k <= kMax; k++) {
			Hz[n] = K1Hz[n] * Hz[n]
			      - K2Hz[n] * (RXc[i] * (Ey[n + Ni] - Ey[n])
			                 - RYc[j] * (Ex[n + Nj] - Ex[n]));
			n++;
		}
	}
	}
}


static void updateHz_f_no_vector(void)
{
	assert(Nk == 1);

	int i;
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (    i = iMin; i <  iMax; i++) {
	for (int j = jMin; j <  jMax; j++) {
		int64_t n = NA(i, j, kMin);
		for (int k = kMin; k <= kMax; k++) {
			const int64_t m = iHz[n];
			Hz[n] = D1[m] * Hz[n]
			      - D2[m] * (RXc[i] * (Ey[n + Ni] - Ey[n])
			               - RYc[j] * (Ex[n + Nj] - Ex[n]));
			n++;
		}
	}
	}
}


static void updateHz_p_vector(double t)
{
	assert(Nk == 1);

	int i;
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (    i = iMin; i <  iMax; i++) {
	for (int j = jMin; j <  jMax; j++) {
		int64_t n = NA(i, j, kMin);
		for (int k = kMin; k <= kMax; k++) {
			real_t fi, dfi;
			finc(Xc[i], Yc[j], Zn[k], t, Planewave.r0, Planewave.ri, Planewave.hi[2], Planewave.ai, Dt, &fi, &dfi);
			Hz[n] = K1Hz[n] * Hz[n]
			      - K2Hz[n] * (RXc[i] * (Ey[n + Ni] - Ey[n])
			                 - RYc[j] * (Ex[n + Nj] - Ex[n]))
			      - (K1Hz[n] - K2Hz[n]) * dfi
			      - (1 - K1Hz[n]) * fi;
			n++;
		}
	}
	}
}


static void updateHz_p_no_vector(double t)
{
	assert(Nk == 1);

	int i;
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (    i = iMin; i <  iMax; i++) {
	for (int j = jMin; j <  jMax; j++) {
		int64_t n = NA(i, j, kMin);
		for (int k = kMin; k <= kMax; k++) {
			const int64_t m = iHz[n];
			if (m == 0) {
				Hz[n] -= RXc[i] * (Ey[n + Ni] - Ey[n])
				       - RYc[j] * (Ex[n + Nj] - Ex[n]);
			}
			else {
				real_t fi, dfi;
				finc(Xc[i], Yc[j], Zn[k], t, Planewave.r0, Planewave.ri, Planewave.hi[2], Planewave.ai, Dt, &fi, &dfi);
				if (m == PEC) {
					Hz[n] = -fi;
				}
				else {
					Hz[n] = D1[m] * Hz[n]
					      - D2[m] * (RXc[i] * (Ey[n + Ni] - Ey[n])
					               - RYc[j] * (Ex[n + Nj] - Ex[n]))
					      - (D1[m] - D2[m]) * dfi
					      - (1 - D1[m]) * fi;
				}
			}
			n++;
		}
	}
	}
}


void updateHz(double t)
{
	if (NFeed) {
		if (VECTOR) {
			updateHz_f_vector();
		}
		else {
			updateHz_f_no_vector();
		}
	}
	else if (IPlanewave) {
		if (VECTOR) {
			updateHz_p_vector(t);
		}
		else {
			updateHz_p_no_vector(t);
		}
	}
}
