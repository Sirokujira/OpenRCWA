/*
setupId.c

setup material index
*/

#include "orcwa.h"
#include "orcwa_prototype.h"

static void boundingbox(int, const double [], double *, double *, double *, double *, double *, double *);

/*
// is air ?
static inline int isair(id_t id)
{
	material_t m = Material[id];

	return (id != PEC)
	    && (m.type == 1)
	    && (!id || (fabs(m.epsr - 1) + fabs(m.esgm) + fabs(m.amur - 1) + fabs(m.msgm) < EPS));
}
*/

/*
// higher material ?
static inline id_t higher(id_t id0, id_t id1)
{
	id_t ret = id0;

	const material_t m0 = Material[id0];
	const material_t m1 = Material[id1];

	if      ((id0 == id1) ||
	         (id0 == PEC) ||
	         (m0.type == 2)) {
		;
	}
	else if ((id1 == PEC) ||
	         (m1.type == 2) ||
	         ((m0.epsr + m0.amur) < (m1.epsr + m1.amur))) {
		ret = id1;
	}

	return ret;
}
*/


// highest material
// type = 1 or 2
static inline id_t highest(id_t id0, id_t id1, id_t id2, id_t id3, id_t id4)
{
	id_t ret = id0;

	if      ((id1 == id0) && (id2 == id0) && (id3 == id0) && (id4 == id0)) {
		;
	}
	else if ((id0 == PEC) || (id1 == PEC) || (id2 == PEC) || (id3 == PEC) || (id4 == PEC)) {
		ret = PEC;
	}
	else if (Material[id1].type == 2) {
		ret = id1;
	}
	else if (Material[id2].type == 2) {
		ret = id2;
	}
	else if (Material[id3].type == 2) {
		ret = id3;
	}
	else if (Material[id4].type == 2) {
		ret = id4;
	}

	return ret;
}


// setup Material ID
void setupId(void)
{
	const double eps = EPS * sqrt(
		(Xn[Nx] - Xn[0]) * (Xn[Nx] - Xn[0]) +
		(Yn[Ny] - Yn[0]) * (Yn[Ny] - Yn[0]) +
		(Zn[Nz] - Zn[0]) * (Zn[Nz] - Zn[0]));
	//printf("%e\n", eps);

	for (int n = 0; n < NGeometry; n++) {
		int shape = Geometry[n].shape;
		id_t m = (id_t)Geometry[n].m;
		double *g = Geometry[n].g;

		// bounding box
		//double x1, x2, y1, y2, z1, z2;
		double x1 = Xn[0];
		double x2 = Xn[Nx];
		double y1 = Yn[0];
		double y2 = Yn[Ny];
		double z1 = Zn[0];
		double z2 = Zn[Nz];
		//if (shape == 1) {
		boundingbox(shape, g, &x1, &x2, &y1, &y2, &z1, &z2);
		//}

		int i1, i2, j1, j2, k1, k2;
		int i;

		// Ex
		if (shape == 1) {
			getspan(Xc, Nx,     iMin, iMax - 1, x1, x2, &i1, &i2, eps);
			getspan(Yn, Ny + 1, jMin, jMax,     y1, y2, &j1, &j2, eps);
			getspan(Zn, Nz + 1, kMin, kMax,     z1, z2, &k1, &k2, eps);
		}
		else {
			i1 = iMin;
			i2 = iMax - 1;
			j1 = jMin;
			j2 = jMax;
			k1 = kMin;
			k2 = kMax;
		}
#ifdef _OPENMP
#pragma omp parallel for
#endif
		for (    i = i1; i <= i2; i++) {
		for (int j = j1; j <= j2; j++) {
		for (int k = k1; k <= k2; k++) {
/*
		for (    i = iMin; i <  iMax; i++) {
		for (int j = jMin; j <= jMax; j++) {
		for (int k = kMin; k <= kMax; k++) {
*/
/*
			if ((iMin <= i) && (i <  iMax) &&
			    (jMin <= j) && (j <= jMax) &&
			    (kMin <= k) && (k <= kMax) &&  // MPI
*/
			if (ingeometry(Xc[i], Yn[j], Zn[k], shape, g, eps)) {
				IEX(i, j, k) = m;
			}
		}
		}
		}

		// Ey
		if (shape == 1) {
			getspan(Yc, Ny,     jMin, jMax - 1, y1, y2, &j1, &j2, eps);
			getspan(Zn, Nz + 1, kMin, kMax,     z1, z2, &k1, &k2, eps);
			getspan(Xn, Nx + 1, iMin, iMax,     x1, x2, &i1, &i2, eps);
		}
		else {
			j1 = jMin;
			j2 = jMax - 1;
			k1 = kMin;
			k2 = kMax;
			i1 = iMin;
			i2 = iMax;
		}
#ifdef _OPENMP
#pragma omp parallel for
#endif
		for (    i = i1; i <= i2; i++) {
		for (int j = j1; j <= j2; j++) {
		for (int k = k1; k <= k2; k++) {
/*
		for (    i = iMin; i <= iMax; i++) {
		for (int j = jMin; j <  jMax; j++) {
		for (int k = kMin; k <= kMax; k++) {
*/
/*
			if ((iMin <= i) && (i <= iMax) &&
			    (jMin <= j) && (j <  jMax) &&
			    (kMin <= k) && (k <= kMax) &&  // MPI
*/
			if (ingeometry(Xn[i], Yc[j], Zn[k], shape, g, eps)) {
				IEY(i, j, k) = m;
			}
		}
		}
		}

		// Ez
		if (shape == 1) {
			getspan(Zc, Nz,     kMin, kMax - 1, z1, z2, &k1, &k2, eps);
			getspan(Xn, Nx + 1, iMin, iMax,     x1, x2, &i1, &i2, eps);
			getspan(Yn, Ny + 1, jMin, jMax,     y1, y2, &j1, &j2, eps);
		}
		else {
			k1 = kMin;
			k2 = kMax - 1;
			i1 = iMin;
			i2 = iMax;
			j1 = jMin;
			j2 = jMax;
		}
#ifdef _OPENMP
#pragma omp parallel for
#endif
		for (    i = i1; i <= i2; i++) {
		for (int j = j1; j <= j2; j++) {
		for (int k = k1; k <= k2; k++) {
/*
		for (    i = iMin; i <= iMax; i++) {
		for (int j = jMin; j <= jMax; j++) {
		for (int k = kMin; k <  kMax; k++) {
*/
/*
			if ((iMin <= i) && (i <= iMax) &&
			    (jMin <= j) && (j <= jMax) &&
			    (kMin <= k) && (k <  kMax) &&  // MPI
*/
			if (ingeometry(Xn[i], Yn[j], Zc[k], shape, g, eps)) {
				IEZ(i, j, k) = m;
			}
		}
		}
		}

		// Hx
		if (shape == 1) {
			getspan(Xn, Nx + 1, iMin,     iMax, x1, x2, &i1, &i2, eps);
			getspan(Yc, Ny,     jMin - 1, jMax, y1, y2, &j1, &j2, eps);
			getspan(Zc, Nz,     kMin - 1, kMax, z1, z2, &k1, &k2, eps);
		}
		else {
			i1 =     iMin;
			i2 =     iMax;
			j1 = MAX(jMin - 1,      0);
			j2 = MIN(jMax,     Ny - 1);
			k1 = MAX(kMin - 1,      0);
			k2 = MIN(kMax,     Nz - 1);
		}
#ifdef _OPENMP
#pragma omp parallel for
#endif
		for (    i = i1; i <= i2; i++) {
		for (int j = j1; j <= j2; j++) {
		for (int k = k1; k <= k2; k++) {
/*
		for (    i = iMin; i <= iMax; i++) {
		for (int j = jMin; j <  jMax; j++) {
		for (int k = kMin; k <  kMax; k++) {
*/
/*
			if ((iMin <= i) && (i <= iMax) &&
			    (jMin <= j) && (j <  jMax) &&
			    (kMin <= k) && (k <  kMax)) {  // MPI
*/
/*
			if ((iMin     <= i) && (i <= iMax) &&
			    (jMin - 1 <= j) && (j <= jMax) &&
			    (kMin - 1 <= k) && (k <= kMax)) {  // MPI
*/
			if (ingeometry(Xn[i], Yc[j], Zc[k], shape, g, eps)) {
				IHX(i, j, k) = m;
			}
/*
			}
*/
		}
		}
		}

		// Hy
		if (shape == 1) {
			getspan(Yn, Ny + 1, jMin,     jMax, y1, y2, &j1, &j2, eps);
			getspan(Zc, Nz,     kMin - 1, kMax, z1, z2, &k1, &k2, eps);
			getspan(Xc, Nx,     iMin - 1, iMax, x1, x2, &i1, &i2, eps);
		}
		else {
			j1 =     jMin;
			j2 =     jMax;
			k1 = MAX(kMin - 1,      0);
			k2 = MIN(kMax,     Nz - 1);
			i1 = MAX(iMin - 1,      0);
			i2 = MIN(iMax,     Nx - 1);
		}
#ifdef _OPENMP
#pragma omp parallel for
#endif
		for (    i = i1; i <= i2; i++) {
		for (int j = j1; j <= j2; j++) {
		for (int k = k1; k <= k2; k++) {
/*
		for (    i = iMin; i <  iMax; i++) {
		for (int j = jMin; j <= jMax; j++) {
		for (int k = kMin; k <  kMax; k++) {
*/
/*
			if ((iMin <= i) && (i <  iMax) &&
			    (jMin <= j) && (j <= jMax) &&
			    (kMin <= k) && (k <  kMax)) {  // MPI
*/
/*
			if ((iMin - 1 <= i) && (i <= iMax) &&
			    (jMin     <= j) && (j <= jMax) &&
			    (kMin - 1 <= k) && (k <= kMax)) {  // MPI
*/
			if (ingeometry(Xc[i], Yn[j], Zc[k], shape, g, eps)) {
				IHY(i, j, k) = m;
			}
/*
			}
*/
		}
		}
		}

		// Hz
		if (shape == 1) {
			getspan(Zn, Nz + 1, kMin,     kMax, z1, z2, &k1, &k2, eps);
			getspan(Xc, Nx,     iMin - 1, iMax, x1, x2, &i1, &i2, eps);
			getspan(Yc, Ny,     jMin - 1, jMax, y1, y2, &j1, &j2, eps);
		}
		else {
			k1 =     kMin;
			k2 =     kMax;
			i1 = MAX(iMin - 1,      0);
			i2 = MIN(iMax,     Nx - 1);
			j1 = MAX(jMin - 1,      0);
			j2 = MIN(jMax,     Ny - 1);
		}
#ifdef _OPENMP
#pragma omp parallel for
#endif
		for (    i = i1; i <= i2; i++) {
		for (int j = j1; j <= j2; j++) {
		for (int k = k1; k <= k2; k++) {
/*
		for (    i = iMin; i <  iMax; i++) {
		for (int j = jMin; j <  jMax; j++) {
		for (int k = kMin; k <= kMax; k++) {
*/
/*
			if ((iMin <= i) && (i <  iMax) &&
			    (jMin <= j) && (j <  jMax) &&
			    (kMin <= k) && (k <= kMax)) {  // MPI
*/
/*
			if ((iMin - 1 <= i) && (i <= iMax) &&
			    (jMin - 1 <= j) && (j <= jMax) &&
			    (kMin     <= k) && (k <= kMax)) {  // MPI
*/
			if (ingeometry(Xc[i], Yc[j], Zn[k], shape, g, eps)) {
				IHZ(i, j, k) = m;
			}
/*
			}
*/
		}
		}
		}

	}
}


// correct surface index
void setupId_surface(void)
{
/*
	// outer H (for the scheme on surfaces)

	// Hx
	for (int i = iMin; i <= iMax; i++) {
		for (int k = kMin; k < kMax; k++) {
			if (inside(2, i, 0,  k)) IHX(i, -1, k) = IEZ(i, 0,  k);
			if (inside(2, i, Ny, k)) IHX(i, Ny, k) = IEZ(i, Ny, k);
		}
		for (int j = jMin; j < jMax; j++) {
			if (inside(1, i, j, 0 )) IHX(i, j, -1) = IEY(i, j, 0 );
			if (inside(1, i, j, Nz)) IHX(i, j, Nz) = IEY(i, j, Nz);
		}
	}
	// Hy
	for (int j = jMin; j <= jMax; j++) {
		for (int i = iMin; i < iMax; i++) {
			if (inside(0, i, j, 0 )) IHY(i, j, -1) = IEX(i, j, 0 );
			if (inside(0, i, j, Nz)) IHY(i, j, Nz) = IEX(i, j, Nz);
		}
		for (int k = kMin; k < kMax; k++) {
			if (inside(2, 0,  j, k)) IHY(-1, j, k) = IEZ(0,  j, k);
			if (inside(2, Nx, j, k)) IHY(Nx, j, k) = IEZ(Nx, j, k);
		}
	}
	// Hz
	for (int k = kMin; k <= kMax; k++) {
		for (int j = jMin; j < jMax; j++) {
			if (inside(1, 0,  j, k)) IHZ(-1, j, k) = IEY(0,  j, k);
			if (inside(1, Nx, j, k)) IHZ(Nx, j, k) = IEY(Nx, j, k);
		}
		for (int i = iMin; i < iMax; i++) {
			if (inside(0, i, 0,  k)) IHZ(i, -1, k) = IEX(i, 0,  k);
			if (inside(0, i, Ny, k)) IHZ(i, Ny, k) = IEX(i, Ny, k);
		}
	}
*/
	// correct curved surface (E <- H)

	//printf("%d %d %d %d %d %d\n", iMin, iMax, jMin, jMax, kMin, kMax); fflush(stdout);
	int i;

	// Ex
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (    i = iMin; i <  iMax; i++) {
	for (int j = jMin; j <= jMax; j++) {
	for (int k = kMin; k <= kMax; k++) {
/*
		if ((IEX(i, j, k) != PEC) && (j > 0) && (j < Ny) && (k > 0) && (k < Nz)) {
			if ((IHY(i,     j,     k    ) == PEC) ||
			    (IHY(i,     j,     k - 1) == PEC) ||
			    (IHZ(i,     j,     k    ) == PEC) ||
			    (IHZ(i,     j - 1, k    ) == PEC)) {
				IEX(i, j, k) = PEC;
			}
		}

		if (isair(IEX(i, j, k)) && (j > 0) && (j < Ny) && (k > 0) && (k < Nz)) {
			if (!isair(IHY(i,     j,     k    ))) IEX(i, j, k) = IHY(i,     j,     k    );
			if (!isair(IHY(i,     j,     k - 1))) IEX(i, j, k) = IHY(i,     j,     k - 1);
			if (!isair(IHZ(i,     j,     k    ))) IEX(i, j, k) = IHZ(i,     j,     k    );
			if (!isair(IHZ(i,     j - 1, k    ))) IEX(i, j, k) = IHZ(i,     j - 1, k    );
		}

		if ((j > 0) && (j < Ny) && (k > 0) && (k < Nz)) {
			IEX(i, j, k) = higher(IEX(i, j, k), IHY(i,     j,     k    ));
			IEX(i, j, k) = higher(IEX(i, j, k), IHY(i,     j,     k - 1));
			IEX(i, j, k) = higher(IEX(i, j, k), IHZ(i,     j,     k    ));
			IEX(i, j, k) = higher(IEX(i, j, k), IHZ(i,     j - 1, k    ));
		}
*/
		if ((j > 0) && (j < Ny) && (k > 0) && (k < Nz)) {
			IEX(i, j, k) = highest(IEX(i, j, k), IHY(i,     j,     k    ),
			                                     IHY(i,     j,     k - 1),
			                                     IHZ(i,     j,     k    ),
			                                     IHZ(i,     j - 1, k    ));
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
/*
		if ((IEY(i, j, k) != PEC) && (k > 0) && (k < Nz) && (i > 0) && (i < Nx)) {
			if ((IHZ(i,     j,     k    ) == PEC) ||
			    (IHZ(i - 1, j,     k    ) == PEC) ||
			    (IHX(i,     j,     k    ) == PEC) ||
			    (IHX(i,     j,     k - 1) == PEC)) {
				IEY(i, j, k) = PEC;
			}
		}

		if (isair(IEY(i, j, k)) && (k > 0) && (k < Nz) && (i > 0) && (i < Nx)) {
			if (!isair(IHZ(i,     j,     k    ))) IEY(i, j, k) = IHZ(i,     j,     k    );
			if (!isair(IHZ(i - 1, j,     k    ))) IEY(i, j, k) = IHZ(i - 1, j,     k    );
			if (!isair(IHX(i,     j,     k    ))) IEY(i, j, k) = IHX(i,     j,     k    );
			if (!isair(IHX(i,     j,     k - 1))) IEY(i, j, k) = IHX(i,     j,     k - 1);
		}

		if ((k > 0) && (k < Nz) && (i > 0) && (i < Nx)) {
			IEY(i, j, k) = higher(IEY(i, j, k), IHZ(i,     j,     k    ));
			IEY(i, j, k) = higher(IEY(i, j, k), IHZ(i - 1, j,     k    ));
			IEY(i, j, k) = higher(IEY(i, j, k), IHX(i,     j,     k    ));
			IEY(i, j, k) = higher(IEY(i, j, k), IHX(i,     j,     k - 1));
		}
*/
		if ((k > 0) && (k < Nz) && (i > 0) && (i < Nx)) {
			IEY(i, j, k) = highest(IEY(i, j, k), IHZ(i,     j,     k    ),
			                                     IHZ(i - 1, j,     k    ),
			                                     IHX(i,     j,     k    ),
			                                     IHX(i,     j,     k - 1));
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
/*
		if ((IEZ(i, j, k) != PEC) && (i > 0) && (i < Nx) && (j > 0) && (j < Ny)) {
			if ((IHX(i,     j,     k    ) == PEC) ||
			    (IHX(i,     j - 1, k    ) == PEC) ||
			    (IHY(i,     j,     k    ) == PEC) ||
			    (IHY(i - 1, j,     k    ) == PEC)) {
				IEZ(i, j, k) = PEC;
			}
		}

		if (isair(IEZ(i, j, k)) && (i > 0) && (i < Nx) && (j > 0) && (j < Ny)) {
			if (!isair(IHX(i,     j,     k    ))) IEZ(i, j, k) = IHX(i,     j,     k    );
			if (!isair(IHX(i,     j - 1, k    ))) IEZ(i, j, k) = IHX(i,     j - 1, k    );
			if (!isair(IHY(i,     j,     k    ))) IEZ(i, j, k) = IHY(i,     j,     k    );
			if (!isair(IHY(i - 1, j,     k    ))) IEZ(i, j, k) = IHY(i - 1, j,     k    );

		if ((i > 0) && (i < Nx) && (j > 0) && (j < Ny)) {
			IEZ(i, j, k) = higher(IEZ(i, j, k), IHX(i,     j,     k    ));
			IEZ(i, j, k) = higher(IEZ(i, j, k), IHX(i,     j - 1, k    ));
			IEZ(i, j, k) = higher(IEZ(i, j, k), IHY(i,     j,     k    ));
			IEZ(i, j, k) = higher(IEZ(i, j, k), IHY(i - 1, j,     k    ));
		}
*/
		if ((i > 0) && (i < Nx) && (j > 0) && (j < Ny)) {
			IEZ(i, j, k) = highest(IEZ(i, j, k), IHX(i,     j,     k    ),
			                                     IHX(i,     j - 1, k    ),
			                                     IHY(i,     j,     k    ),
			                                     IHY(i - 1, j,     k    ));
		}
	}
	}
	}
}


// get bounding box
static void boundingbox(int shape, const double g[],
	double *x1, double *x2, double *y1, double *y2, double *z1, double *z2)
{
	//*x1 = *x2 = *y1 = *y2 = *z1 = *z2 = 0;

	// cube / sphere / cylinder
	if ((shape == 1) || (shape == 2) ||
	    (shape == 11) || (shape == 12) || (shape == 13)) {
		*x1 = MIN(g[0], g[1]);
		*x2 = MAX(g[0], g[1]);
		*y1 = MIN(g[2], g[3]);
		*y2 = MAX(g[2], g[3]);
		*z1 = MIN(g[4], g[5]);
		*z2 = MAX(g[4], g[5]);
	}
	// triangle pillar
	else if (shape == 31) {
		*x1 = MIN(g[0], g[1]);
		*x2 = MAX(g[0], g[1]);
		*y1 = MIN(g[2], MIN(g[3], g[4]));
		*y2 = MAX(g[2], MAX(g[3], g[4]));
		*z1 = MIN(g[5], MIN(g[6], g[7]));
		*z2 = MAX(g[5], MAX(g[6], g[7]));
	}
	else if (shape == 32) {
		*y1 = MIN(g[0], g[1]);
		*y2 = MAX(g[0], g[1]);
		*z1 = MIN(g[2], MIN(g[3], g[4]));
		*z2 = MAX(g[2], MAX(g[3], g[4]));
		*x1 = MIN(g[5], MIN(g[6], g[7]));
		*x2 = MAX(g[5], MAX(g[6], g[7]));
	}
	else if (shape == 33) {
		*z1 = MIN(g[0], g[1]);
		*z2 = MAX(g[0], g[1]);
		*x1 = MIN(g[2], MIN(g[3], g[4]));
		*x2 = MAX(g[2], MAX(g[3], g[4]));
		*y1 = MIN(g[5], MIN(g[6], g[7]));
		*y2 = MAX(g[5], MAX(g[6], g[7]));
	}
	// pyramid / cone
	else if ((shape == 41) || (shape == 51)) {
		const double hy = MAX(g[4], g[6]) / 2;
		const double hz = MAX(g[5], g[7]) / 2;
		*x1 = MIN(g[0], g[1]);
		*x2 = MAX(g[0], g[1]);
		*y1 = g[2] - hy;
		*y2 = g[2] + hy;
		*z1 = g[3] - hz;
		*z2 = g[3] + hz;
	}
	else if ((shape == 42) || (shape == 52)) {
		const double hz = MAX(g[4], g[6]) / 2;
		const double hx = MAX(g[5], g[7]) / 2;
		*y1 = MIN(g[0], g[1]);
		*y2 = MAX(g[0], g[1]);
		*z1 = g[2] - hz;
		*z2 = g[2] + hz;
		*x1 = g[3] - hx;
		*x2 = g[3] + hx;
	}
	else if ((shape == 43) || (shape == 53)) {
		const double hx = MAX(g[4], g[6]) / 2;
		const double hy = MAX(g[5], g[7]) / 2;
		*z1 = MIN(g[0], g[1]);
		*z2 = MAX(g[0], g[1]);
		*x1 = g[2] - hx;
		*x2 = g[2] + hx;
		*y1 = g[3] - hy;
		*y2 = g[3] + hy;
	}
}


// debug
void debugId(void)
{

	for (int k = 0; k <= Nz; k++) {
		printf("Ex k=%d\n", k);
		for (int j = Ny; j >= 0; j--) {
			for (int i = 0; i <  Nx; i++) {
				printf("%d", IEX(i, j, k));
			}
			printf("\n");
		}
	}

	for (int k = 0; k <= Nz; k++) {
		printf("Ey k=%d\n", k);
		for (int j = Ny - 1; j >= 0; j--) {
			for (int i = 0; i <= Nx; i++) {
				printf("%d", IEY(i, j, k));
			}
			printf("\n");
		}
	}

	for (int k = 0; k < Nz; k++) {
		printf("Ez k=%d\n", k);
		for (int j = Ny; j >= 0; j--) {
			for (int i = 0; i <= Nx; i++) {
				printf("%d", IEZ(i, j, k));
			}
			printf("\n");
		}
	}

	for (int k = 0; k < Nz; k++) {
		printf("Hx k=%d\n", k);
		for (int j = Ny - 1; j >= 0; j--) {
			for (int i = 0; i <= Nx; i++) {
				printf("%d", IHX(i, j, k));
			}
			printf("\n");
		}
	}

	for (int k = 0; k < Nz; k++) {
		printf("Hy k=%d\n", k);
		for (int j = Ny; j >= 0; j--) {
			for (int i = 0; i < Nx; i++) {
				printf("%d", IHY(i, j, k));
			}
			printf("\n");
		}
	}

	for (int k = 0; k <= Nz; k++) {
		printf("Hz k=%d\n", k);
		for (int j = Ny - 1; j >= 0; j--) {
			for (int i = 0; i < Nx; i++) {
				printf("%d", IHZ(i, j, k));
			}
			printf("\n");
		}
	}

	fflush(stdout);
}


void debugId2(void)
{
	int isum[6] = {0, 0, 0, 0, 0, 0};
	for (int n = 0; n < NN; n++) {
		isum[0] += iEx[n];
		isum[1] += iEy[n];
		isum[2] += iEz[n];
		isum[3] += iHx[n];
		isum[4] += iHy[n];
		isum[5] += iHz[n];
	}
	printf("ID: %d %d %d %d %d %d\n", isum[0], isum[1], isum[2], isum[3], isum[4], isum[5]);
}
