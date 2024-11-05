/*
dftNear3d.c (OpenMP)

DFT of near field
*/

#include "orcwa.h"

void dftNear3d(int itime)
{
	for (int ifreq = 0; ifreq < NFreq2; ifreq++) {
		const int64_t n0 = ifreq * NN;
		const int id = (itime * NFreq2) + ifreq;

		const float ef_r = (float)cEdft[id].r;
		const float ef_i = (float)cEdft[id].i;
		const float hf_r = (float)cHdft[id].r;
		const float hf_i = (float)cHdft[id].i;

		int i;

		// Ex
#ifdef _OPENMP
#pragma omp parallel for
#endif
		for (    i = iMin; i < iMax + 0; i++) {
		for (int j = jMin; j < jMax + 1; j++) {
		for (int k = kMin; k < kMax + 1; k++) {
			const int64_t n = NA(i, j, k);
			cEx_r[n0 + n] += (float)Ex[n] * ef_r;
			cEx_i[n0 + n] += (float)Ex[n] * ef_i;
		}
		}
		}

		// Ey
#ifdef _OPENMP
#pragma omp parallel for
#endif
		for (    i = iMin; i < iMax + 1; i++) {
		for (int j = jMin; j < jMax + 0; j++) {
		for (int k = kMin; k < kMax + 1; k++) {
			const int64_t n = NA(i, j, k);
			cEy_r[n0 + n] += (float)Ey[n] * ef_r;
			cEy_i[n0 + n] += (float)Ey[n] * ef_i;
		}
		}
		}

		// Ez
#ifdef _OPENMP
#pragma omp parallel for
#endif
		for (    i = iMin; i < iMax + 1; i++) {
		for (int j = jMin; j < jMax + 1; j++) {
		for (int k = kMin; k < kMax + 0; k++) {
			const int64_t n = NA(i, j, k);
			cEz_r[n0 + n] += (float)Ez[n] * ef_r;
			cEz_i[n0 + n] += (float)Ez[n] * ef_i;
		}
		}
		}

		// Hx
#ifdef _OPENMP
#pragma omp parallel for
#endif
		for (    i = iMin - 1; i < iMax + 2; i++) {
		for (int j = jMin - 1; j < jMax + 1; j++) {
		for (int k = kMin - 1; k < kMax + 1; k++) {
			const int64_t n = NA(i, j, k);
			cHx_r[n0 + n] += (float)Hx[n] * hf_r;
			cHx_i[n0 + n] += (float)Hx[n] * hf_i;
		}
		}
		}

		// Hy
#ifdef _OPENMP
#pragma omp parallel for
#endif
		for (    i = iMin - 1; i < iMax + 1; i++) {
		for (int j = jMin - 1; j < jMax + 2; j++) {
		for (int k = kMin - 1; k < kMax + 1; k++) {
			const int64_t n = NA(i, j, k);
			cHy_r[n0 + n] += (float)Hy[n] * hf_r;
			cHy_i[n0 + n] += (float)Hy[n] * hf_i;
		}
		}
		}

		// Hz
#ifdef _OPENMP
#pragma omp parallel for
#endif
		for (    i = iMin - 1; i < iMax + 1; i++) {
		for (int j = jMin - 1; j < jMax + 1; j++) {
		for (int k = kMin - 1; k < kMax + 2; k++) {
			const int64_t n = NA(i, j, k);
			cHz_r[n0 + n] += (float)Hz[n] * hf_r;
			cHz_i[n0 + n] += (float)Hz[n] * hf_i;
		}
		}
		}
	}
}
