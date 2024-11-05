/*
setupDispersion.c

setup dispersion data
*/

#include "orcwa.h"

static void setupDispersion_id(int);


void setupDispersion(void)
{
	for (int m = 0; m < NMaterial; m++) {
		if (Material[m].type == 2) {
			const double einf = Material[m].einf;
			const double ae   = Material[m].ae;
			const double be   = Material[m].be;
			const double ce   = Material[m].ce;
			const double ke = exp(-ce * Dt);
			const double xi0 = (ae * Dt) + (be / ce) * (1 - ke);
			const double dxi0 = (be / ce) * (1 - ke) * (1 - ke);
			//printf("%d %e %e %e %e\n", m, einf, ae, be, ce);
			//printf("%d %e %e %e\n", m, ke, xi0, dxi0);
			Material[m].edisp[0] = (real_t)(1 / (einf + xi0));
			Material[m].edisp[1] = (real_t)dxi0;
			Material[m].edisp[2] = (real_t)ke;
			//printf("%d %e %e %e\n", m, Material[m].edisp[0], Material[m].edisp[1], Material[m].edisp[2]);
		}
	}

	// setup dispersion id

	numDispersionEx = numDispersionEy = numDispersionEz = 0;
	setupDispersion_id(0);
	//printf("%zd %zd %zd\n", numDispersionEx, numDispersionEy, numDispersionEz);

	if (numDispersionEx > 0) {
		mDispersionEx = (dispersion_t *)malloc(numDispersionEx * sizeof(dispersion_t));
		DispersionEx =        (real_t *)malloc(numDispersionEx * sizeof(real_t));
		memset(DispersionEx, 0,                numDispersionEx * sizeof(real_t));
	}
	if (numDispersionEy > 0) {
		mDispersionEy = (dispersion_t *)malloc(numDispersionEy * sizeof(dispersion_t));
		DispersionEy =        (real_t *)malloc(numDispersionEy * sizeof(real_t));
		memset(DispersionEy, 0,                numDispersionEy * sizeof(real_t));
	}
	if (numDispersionEz > 0) {
		mDispersionEz = (dispersion_t *)malloc(numDispersionEz * sizeof(dispersion_t));
		DispersionEz =        (real_t *)malloc(numDispersionEz * sizeof(real_t));
		memset(DispersionEz, 0,                numDispersionEz * sizeof(real_t));
	}

	setupDispersion_id(1);
}


// mode = 0/1  (OpenMP : NG for CUDA)
static void setupDispersion_id(int mode)
{
	// Ex
	int64_t nx = 0;
	for (int i = iMin; i <  iMax; i++) {
	for (int j = jMin; j <= jMax; j++) {
	for (int k = kMin; k <= kMax; k++) {
		const id_t m = IEX(i, j, k);
		if (Material[m].type == 2) {
			if (mode == 1) {
				mDispersionEx[nx].i = i;
				mDispersionEx[nx].j = j;
				mDispersionEx[nx].k = k;
				mDispersionEx[nx].f1 = Material[m].edisp[0];
				mDispersionEx[nx].f2 = Material[m].edisp[1];
				mDispersionEx[nx].f3 = Material[m].edisp[2];
			}
			nx++;
		}
	}
	}
	}
	if (mode == 0) {
		numDispersionEx = nx;
	}

	// Ey
	int64_t ny = 0;
	for (int i = iMin; i <= iMax; i++) {
	for (int j = jMin; j <  jMax; j++) {
	for (int k = kMin; k <= kMax; k++) {
		const id_t m = IEY(i, j, k);
		if (Material[m].type == 2) {
			if (mode == 1) {
				mDispersionEy[ny].i = i;
				mDispersionEy[ny].j = j;
				mDispersionEy[ny].k = k;
				mDispersionEy[ny].f1 = Material[m].edisp[0];
				mDispersionEy[ny].f2 = Material[m].edisp[1];
				mDispersionEy[ny].f3 = Material[m].edisp[2];
			}
			ny++;
		}
	}
	}
	}
	if (mode == 0) {
		numDispersionEy = ny;
	}

	// Ez
	int64_t nz = 0;
	for (int i = iMin; i <= iMax; i++) {
	for (int j = jMin; j <= jMax; j++) {
	for (int k = kMin; k <  kMax; k++) {
		const id_t m = IEZ(i, j, k);
		if (Material[m].type == 2) {
			if (mode == 1) {
				mDispersionEz[nz].i = i;
				mDispersionEz[nz].j = j;
				mDispersionEz[nz].k = k;
				mDispersionEz[nz].f1 = Material[m].edisp[0];
				mDispersionEz[nz].f2 = Material[m].edisp[1];
				mDispersionEz[nz].f3 = Material[m].edisp[2];
			}
			nz++;
		}
	}
	}
	}
	if (mode == 0) {
		numDispersionEz = nz;
	}
}
