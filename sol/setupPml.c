/*
setupPml.c

setup PML factor
*/

#include "orcwa.h"

void setupPml(void)
{
	double fc, d;

	const double kpml = (cPML.m + 1) / (2.0 * cPML.l) * log(1 / cPML.r0);

	double *f = (double *)malloc((2 * cPML.l + 1) * sizeof(double));

	f[0] = f[1] = 0;
	for (int n = 2; n <= (2 * cPML.l); n++) {
		f[n] = kpml * pow((n - 1) / (2.0 * cPML.l), cPML.m);
	}

	const int lx = cPML.l;
	const int ly = cPML.l;
	const int lz = cPML.l;

	const double cdt = C * Dt;

	// E
	for (int i = -lx + 1; i < Nx + lx; i++) {
		if      (i <= 0 ) {fc = f[- 2 * (i     )];     d = Xn[1    ] - Xn[0     ];}
		else if (i >= Nx) {fc = f[+ 2 * (i - Nx)];     d = Xn[Nx   ] - Xn[Nx - 1];}
		else              {fc = 0;                     d = Xc[i    ] - Xc[i  - 1];}
		//gPmlXn[i + lx] = (real_t)(1 / (1 + cdt / d * fc));
		gPmlXn[i + lx] = (real_t)(cdt / d * fc);
	}
	for (int j = -ly + 1; j < Ny + ly; j++) {
		if      (j <= 0 ) {fc = f[- 2 * (j     )];     d = Yn[1    ] - Yn[0     ];}
		else if (j >= Ny) {fc = f[+ 2 * (j - Ny)];     d = Yn[Ny   ] - Yn[Ny - 1];}
		else              {fc = 0;                     d = Yc[j    ] - Yc[j  - 1];}
		//gPmlYn[j + ly] = (real_t)(1 / (1 + cdt / d * fc));
		gPmlYn[j + ly] = (real_t)(cdt / d * fc);
	}
	for (int k = -lz + 1; k < Nz + lz; k++) {
		if      (k <= 0 ) {fc = f[- 2 * (k     )];     d = Zn[1    ] - Zn[0     ];}
		else if (k >= Nz) {fc = f[+ 2 * (k - Nz)];     d = Zn[Nz   ] - Zn[Nz - 1];}
		else              {fc = 0;                     d = Zc[k    ] - Zc[k  - 1];}
		//gPmlZn[k + lz] = (real_t)(1 / (1 + cdt / d * fc));
		gPmlZn[k + lz] = (real_t)(cdt / d * fc);
	}

	// H
	for (int i = -lx;     i < Nx + lx; i++) {
		if      (i <  0 ) {fc = f[- 2 * (i     ) - 1]; d = Xn[1    ] - Xn[0     ];}
		else if (i >= Nx) {fc = f[+ 2 * (i - Nx) + 1]; d = Xn[Nx   ] - Xn[Nx - 1];}
		else              {fc = 0;                     d = Xn[i + 1] - Xn[i     ];}
		//gPmlXc[i + lx] = (real_t)(1 / (1 + cdt / d * fc));
		gPmlXc[i + lx] = (real_t)(cdt / d * fc);
	}
	for (int j = -ly;     j < Ny + ly; j++) {
		if      (j <  0 ) {fc = f[- 2 * (j     ) - 1]; d = Yn[1    ] - Yn[0     ];}
		else if (j >= Ny) {fc = f[+ 2 * (j - Ny) + 1]; d = Yn[Ny   ] - Yn[Ny - 1];}
		else              {fc = 0;                     d = Yn[j + 1] - Yn[j     ];}
		//gPmlYc[j + ly] = (real_t)(1 / (1 + cdt / d * fc));
		gPmlYc[j + ly] = (real_t)(cdt / d * fc);
	}
	for (int k = -lz;     k < Nz + lz; k++) {
		if      (k <  0 ) {fc = f[- 2 * (k     ) - 1]; d = Zn[1    ] - Zn[0     ];}
		else if (k >= Nz) {fc = f[+ 2 * (k - Nz) + 1]; d = Zn[Nz   ] - Zn[Nz - 1];}
		else              {fc = 0;                     d = Zn[k + 1] - Zn[k     ];}
		//gPmlZc[k + lz] = (real_t)(1 / (1 + cdt / d * fc));
		gPmlZc[k + lz] = (real_t)(cdt / d * fc);
	}

	free(f);
/*
	// debug
	for (int i = 0; i < Nx + 2 * lx; i++) {
		printf("X : %d %e %e\n", i, gPmlXn[i], gPmlXc[i]);
	}
	for (int j = 0; j < Ny + 2 * ly; j++) {
		printf("Y : %d %e %e\n", j, gPmlYn[j], gPmlYc[j]);
	}
	for (int k = 0; k < Nz + 2 * lz; k++) {
		printf("Z : %d %e %e\n", k, gPmlZn[k], gPmlZc[k]);
	}
	fflush(stdout); exit(0);
*/
	// boundary material factor (reverse)
	for (int64_t m = 0; m < NMaterial; m++) {
		rPmlE[m] = 1 / (real_t)Material[m].epsr;
		rPmlH[m] = 1 / (real_t)Material[m].amur;
		rPml[m] = 1 / (real_t)sqrt(Material[m].epsr * Material[m].amur);
	}

	// do not delete
	rPmlE[PEC] = 0;
	rPmlH[PEC] = 0;
	rPml[PEC] = 0;
}
