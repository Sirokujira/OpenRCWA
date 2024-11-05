/*
updateEx.c

update Ex
*/

#include "orcwa.h"
#include "finc.h"

static void updateEx_f_vector(void)
{
	assert(Nk == 1);

	int i;
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (    i = iMin; i <  iMax; i++) {
	for (int j = jMin; j <= jMax; j++) {
		int64_t n = NA(i, j, kMin);
		for (int k = kMin; k <= kMax; k++) {
			Ex[n] = K1Ex[n] * Ex[n]
			      + K2Ex[n] * (RYn[j] * (Hz[n] - Hz[n - Nj])
			                 - RZn[k] * (Hy[n] - Hy[n - Nk]));
			n++;
		}
	}
	}
}


static void updateEx_f_no_vector(void)
{
	assert(Nk == 1);

	int i;
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (    i = iMin; i <  iMax; i++) {
	for (int j = jMin; j <= jMax; j++) {
		int64_t n = NA(i, j, kMin);
		for (int k = kMin; k <= kMax; k++) {
			const int64_t m = iEx[n];
			Ex[n] = C1[m] * Ex[n]
			      + C2[m] * (RYn[j] * (Hz[n] - Hz[n - Nj])
			               - RZn[k] * (Hy[n] - Hy[n - Nk]));
			n++;
		}
	}
	}
}


static void updateEx_p_vector(double t)
{
	assert(Nk == 1);

	int i;
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (    i = iMin; i <  iMax; i++) {
	for (int j = jMin; j <= jMax; j++) {
		int64_t n = NA(i, j, kMin);
		for (int k = kMin; k <= kMax; k++) {
			real_t fi, dfi;
			finc(Xc[i], Yn[j], Zn[k], t, Planewave.r0, Planewave.ri, Planewave.ei[0], Planewave.ai, Dt, &fi, &dfi);
			Ex[n] = K1Ex[n] * Ex[n]
			      + K2Ex[n] * (RYn[j] * (Hz[n] - Hz[n - Nj])
			                 - RZn[k] * (Hy[n] - Hy[n - Nk]))
			      - (K1Ex[n] - K2Ex[n]) * dfi
			      - (1 - K1Ex[n]) * fi;
			n++;
		}
	}
	}
}


static void updateEx_p_no_vector(double t)
{
	assert(Nk == 1);

	int i;
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (    i = iMin; i <  iMax; i++) {
	for (int j = jMin; j <= jMax; j++) {
		int64_t n = NA(i, j, kMin);
		for (int k = kMin; k <= kMax; k++) {
			const int64_t m = iEx[n];
			if (m == 0) {
				Ex[n] += RYn[j] * (Hz[n] - Hz[n - Nj])
				       - RZn[k] * (Hy[n] - Hy[n - Nk]);
			}
			else {
				real_t fi, dfi;
				finc(Xc[i], Yn[j], Zn[k], t, Planewave.r0, Planewave.ri, Planewave.ei[0], Planewave.ai, Dt, &fi, &dfi);
				if (m == PEC) {
					Ex[n] = -fi;
				}
				else {
					Ex[n] = C1[m] * Ex[n]
					      + C2[m] * (RYn[j] * (Hz[n] - Hz[n - Nj])
					               - RZn[k] * (Hy[n] - Hy[n - Nk]))
					      - (C1[m] - C2[m]) * dfi
					      - (1 - C1[m]) * fi;
				}
			}
			n++;
		}
	}
	}
}


void updateEx(double t)
{
	if (NFeed) {
		if (VECTOR) {
			updateEx_f_vector();
		}
		else {
			updateEx_f_no_vector();
		}
	}
	else if (IPlanewave) {
		if (VECTOR) {
			updateEx_p_vector(t);
		}
		else {
			updateEx_p_no_vector(t);
		}
	}
}
