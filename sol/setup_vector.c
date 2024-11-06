/*
setup_vector.c

setup (VECTOR): direct access
*/

#include "orcwa.h"

/*
// higher material ?
// t0, t1 = type (1 or 2)
// e0, e1 = epsr + amur
static inline id_t higher_vector(id_t id0, id_t id1, real_t t0, real_t t1, real_t e0, real_t e1)
{
	id_t ret = id0;

	if      ((id0 == id1) ||
	         (id0 == PEC) ||
	         (fabs(t0 - 2) < EPS)) {
		;
	}
	else if ((id1 == PEC) ||
	         (fabs(t1 - 2) < EPS) ||
	         (e0 < e1)) {
		ret = id1;
	}

	return ret;
}
*/


// highest material
// t1, t2, t3, t4 = type (1 or 2)
static inline id_t highest_vector(
	id_t id0, id_t id1, id_t id2, id_t id3, id_t id4,
	real_t t1, real_t t2, real_t t3, real_t t4)
{
	id_t ret = id0;

	if      ((id1 == id0) && (id2 == id0) && (id3 == id0) && (id4 == id0)) {
		;
	}
	else if ((id0 == PEC) || (id1 == PEC) || (id2 == PEC) || (id3 == PEC) || (id4 == PEC)) {
		ret = PEC;
	}
	else if (t1 > 1.5) {
		ret = id1;
	}
	else if (t2 > 1.5) {
		ret = id2;
	}
	else if (t3 > 1.5) {
		ret = id3;
	}
	else if (t4 > 1.5) {
		ret = id4;
	}

	return ret;
}


// setup material parameters (vector)
void setupId_vector(void)
{
	// clear
	size_t size = NN * sizeof(real_t);
	memset(K1Ex, 0, size);
	memset(K1Ey, 0, size);
	memset(K1Ez, 0, size);
	memset(K1Hx, 0, size);
	memset(K1Hy, 0, size);
	memset(K1Hz, 0, size);

	int i;

	// E
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (    i = iMin; i <  iMax; i++) {
	for (int j = jMin; j <= jMax; j++) {
	for (int k = kMin; k <= kMax; k++) {
		const int64_t n = NA(i, j, k);
		const int64_t m = IEX(i, j, k);
		K1Ex[n] = (real_t)Material[m].type;
	}
	}
	}
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (    i = iMin; i <= iMax; i++) {
	for (int j = jMin; j <  jMax; j++) {
	for (int k = kMin; k <= kMax; k++) {
		const int64_t n = NA(i, j, k);
		const int64_t m = IEY(i, j, k);
		K1Ey[n] = (real_t)Material[m].type;
	}
	}
	}
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (    i = iMin; i <= iMax; i++) {
	for (int j = jMin; j <= jMax; j++) {
	for (int k = kMin; k <  kMax; k++) {
		const int64_t n = NA(i, j, k);
		const int64_t m = IEZ(i, j, k);
		K1Ez[n] = (real_t)Material[m].type;
	}
	}
	}

	// H
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (    i = iMin; i <= iMax; i++) {
	for (int j = jMin; j <  jMax; j++) {
	for (int k = kMin; k <  kMax; k++) {
		const int64_t n = NA(i, j, k);
		const int64_t m = IHX(i, j, k);
		K1Hx[n] = (real_t)Material[m].type;
	}
	}
	}
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (    i = iMin; i <  iMax; i++) {
	for (int j = jMin; j <= jMax; j++) {
	for (int k = kMin; k <  kMax; k++) {
		const int64_t n = NA(i, j, k);
		const int64_t m = IHY(i, j, k);
		K1Hy[n] = (real_t)Material[m].type;
	}
	}
	}
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (    i = iMin; i <  iMax; i++) {
	for (int j = jMin; j <  jMax; j++) {
	for (int k = kMin; k <= kMax; k++) {
		const int64_t n = NA(i, j, k);
		const int64_t m = IHZ(i, j, k);
		K1Hz[n] = (real_t)Material[m].type;
	}
	}
	}
}


// correct surface index (vector)
void setupId_surface_vector(void)
{
	int i;

	// Ex
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (    i = iMin; i <  iMax; i++) {
	for (int j = jMin; j <= jMax; j++) {
	for (int k = kMin; k <= kMax; k++) {
		if ((j > 0) && (j < Ny) && (k > 0) && (k < Nz)) {
			IEX(i, j, k) = highest_vector(IEX(i, j, k),
				    IHY(i,     j,     k    ),      IHY(i,     j,     k - 1),      IHZ(i,     j,     k    ),      IHZ(i,     j - 1, k    ),
				K1Hy[NA(i,     j,     k    )], K1Hy[NA(i,     j,     k - 1)], K1Hz[NA(i,     j,     k    )], K1Hz[NA(i,     j - 1, k    )]);
		}
	}
	}
	}

	// Ey
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (    i = iMin; i <= iMax; i++) {
	for (int j = jMin; j <  jMax; j++) {
	for (int k = kMin; k <= kMax; k++) {
		if ((k > 0) && (k < Nz) && (i > 0) && (i < Nx)) {
			IEY(i, j, k) = highest_vector(IEY(i, j, k),
				    IHZ(i,     j,     k    ),      IHZ(i - 1, j,     k    ),      IHX(i,     j,     k    ),      IHX(i,     j,     k - 1),
				K1Hz[NA(i,     j,     k    )], K1Hz[NA(i - 1, j,     k    )], K1Hx[NA(i,     j,     k    )], K1Hx[NA(i,     j,     k - 1)]);
		}
	}
	}
	}

	// Ez
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (    i = iMin; i <= iMax; i++) {
	for (int j = jMin; j <= jMax; j++) {
	for (int k = kMin; k <  kMax; k++) {
		if ((i > 0) && (i < Nx) && (j > 0) && (j < Ny)) {
			IEZ(i, j, k) = highest_vector(IEZ(i, j, k),
				    IHX(i,     j,     k    ),      IHX(i,     j - 1, k    ),      IHY(i,     j,     k    ),      IHY(i - 1, j,     k    ),
				K1Hx[NA(i,     j,     k    )], K1Hx[NA(i,     j - 1, k    )], K1Hy[NA(i,     j,     k    )], K1Hy[NA(i - 1, j,     k    )]);
		}
	}
	}
	}

}


// setup material ID (vector, reuse)
void setup_material_vector(void)
{
	// clear again
	size_t size = NN * sizeof(real_t);
	memset(K1Ex, 0, size);
	memset(K1Ey, 0, size);
	memset(K1Ez, 0, size);
	memset(K1Hx, 0, size);
	memset(K1Hy, 0, size);
	memset(K1Hz, 0, size);

	int i;

#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (    i = iMin; i <  iMax; i++) {
	for (int j = jMin; j <= jMax; j++) {
	for (int k = kMin; k <= kMax; k++) {
		const int64_t n = NA(i, j, k);
		const int64_t m = IEX(i, j, k);
		K1Ex[n] = C1[m];
		K2Ex[n] = C2[m];
	}
	}
	}
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (    i = iMin; i <= iMax; i++) {
	for (int j = jMin; j <  jMax; j++) {
	for (int k = kMin; k <= kMax; k++) {
		const int64_t n = NA(i, j, k);
		const int64_t m = IEY(i, j, k);
		K1Ey[n] = C1[m];
		K2Ey[n] = C2[m];
	}
	}
	}
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (    i = iMin; i <= iMax; i++) {
	for (int j = jMin; j <= jMax; j++) {
	for (int k = kMin; k <  kMax; k++) {
		const int64_t n = NA(i, j, k);
		const int64_t m = IEZ(i, j, k);
		K1Ez[n] = C1[m];
		K2Ez[n] = C2[m];
	}
	}
	}

#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (    i = iMin; i <= iMax; i++) {
	for (int j = jMin; j <  jMax; j++) {
	for (int k = kMin; k <  kMax; k++) {
		const int64_t n = NA(i, j, k);
		const int64_t m = IHX(i, j, k);
		K1Hx[n] = D1[m];
		K2Hx[n] = D2[m];
	}
	}
	}
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (    i = iMin; i <  iMax; i++) {
	for (int j = jMin; j <= jMax; j++) {
	for (int k = kMin; k <  kMax; k++) {
		const int64_t n = NA(i, j, k);
		const int64_t m = IHY(i, j, k);
		K1Hy[n] = D1[m];
		K2Hy[n] = D2[m];
	}
	}
	}
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (    i = iMin; i <  iMax; i++) {
	for (int j = jMin; j <  jMax; j++) {
	for (int k = kMin; k <= kMax; k++) {
		const int64_t n = NA(i, j, k);
		const int64_t m = IHZ(i, j, k);
		K1Hz[n] = D1[m];
		K2Hz[n] = D2[m];
	}
	}
	}
}
