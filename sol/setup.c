/*
setup.c

setup
*/

#include "orcwa.h"
#include "orcwa_prototype.h"

extern void setupId_vector(void);
extern void setupId_surface_vector(void);
extern void setup_material_vector(void);
extern void setupId_surface(void);
//extern void debugId(void);
//extern void debugId2(void);

// time step
static void setup_timestep(void)
{
	double dxmin = fabs(Xn[1] - Xn[0]);
	for (int i = 1; i < Nx; i++) {
		dxmin = MIN(dxmin, fabs(Xn[i + 1] - Xn[i]));
	}

	double dymin = fabs(Yn[1] - Yn[0]);
	for (int j = 1; j < Ny; j++) {
		dymin = MIN(dymin, fabs(Yn[j + 1] - Yn[j]));
	}

	double dzmin = fabs(Zn[1] - Zn[0]);
	for (int k = 1; k < Nz; k++) {
		dzmin = MIN(dzmin, fabs(Zn[k + 1] - Zn[k]));
	}

	Dt = 1 / sqrt(1 / (dxmin * dxmin) + 1 / (dymin * dymin) + 1 / (dzmin * dzmin)) / C;
	//printf("%e %e %e %e\n", dxmin, dymin, dzmin, Dt);
}


// pulse width
static void setup_pulsewidth(void)
{
	double f0;
	if      (NFreq1 > 0) {
		f0 = (Freq1[0] + Freq1[NFreq1 - 1]) / 2;
	}
	else if (NFreq2 > 0) {
		f0 = (Freq2[0] + Freq2[NFreq2 - 1]) / 2;
	}
	else {
		f0 = 1 / (20 * Dt);
	}

	Tw = 1.27 / f0;
	//printf("f0=%e Tw=%e\n", f0, Tw);
}


// material parameters
static void setup_material(void)
{
	// air
	C1[0] = 1;
	C2[0] = 1;
	D1[0] = 1;
	D2[0] = 1;

	// PEC
	C1[1] = 0;
	C2[1] = 0;
	D1[1] = 0;
	D2[1] = 0;

	for (int m = 2; m < NMaterial; m++) {
		if      (Material[m].type == 1) {
			const double denom = Material[m].epsr + (Material[m].esgm * ETA0 * C * Dt);
			C1[m] = (real_t)(Material[m].epsr / denom);
			C2[m] = (real_t)(1 / denom);
		}
		else if (Material[m].type == 2) {
			const double einf = Material[m].einf;
			const double ae   = Material[m].ae;
			const double be   = Material[m].be;
			const double ce   = Material[m].ce;
			const double ke = exp(-ce * Dt);
			const double xi0 = (ae * Dt) + (be / ce) * (1 - ke);
			const double denom = einf + xi0;
			C1[m] = (real_t)(einf / denom);
			C2[m] = (real_t)(1    / denom);
		}
	}

	for (int m = 2; m < NMaterial; m++) {
		const double denom = Material[m].amur + (Material[m].msgm / ETA0 * C * Dt);
		D1[m] = (real_t)(Material[m].amur / denom);
		D2[m] = (real_t)(1 / denom);
	}
/*
	// debug
	for (int m = 0; m < NMaterial; m++) {
		printf("%d %e %e %e %e\n", m, C1[m], C2[m], D1[m], D2[m]);
	}
*/
}


// mesh factor : c * dt / d
static void setup_mesh(void)
{
	double cdt = C * Dt;

	for (int i = 1; i < Nx; i++) {
		RXn[i] = (real_t)(cdt / ((Xn[i + 1] - Xn[i - 1]) / 2));
	}
	for (int j = 1; j < Ny; j++) {
		RYn[j] = (real_t)(cdt / ((Yn[j + 1] - Yn[j - 1]) / 2));
	}
	for (int k = 1; k < Nz; k++) {
		RZn[k] = (real_t)(cdt / ((Zn[k + 1] - Zn[k - 1]) / 2));
	}

	RXn[0]  = (real_t)(cdt / (Xn[1]  - Xn[0]));
	RYn[0]  = (real_t)(cdt / (Yn[1]  - Yn[0]));
	RZn[0]  = (real_t)(cdt / (Zn[1]  - Zn[0]));
	RXn[Nx] = (real_t)(cdt / (Xn[Nx] - Xn[Nx - 1]));
	RYn[Ny] = (real_t)(cdt / (Yn[Ny] - Yn[Ny - 1]));
	RZn[Nz] = (real_t)(cdt / (Zn[Nz] - Zn[Nz - 1]));
/*
	// debug
	for (int i = 0; i <= Nx; i++) {
		printf("RXn[%d]=%.5f\n", i, RXn[i]);
	}
	for (int j = 0; j <= Ny; j++) {
		printf("RYn[%d]=%.5f\n", j, RYn[j]);
	}
	for (int k = 0; k <= Nz; k++) {
		printf("RZn[%d]=%.5f\n", k, RZn[k]);
	}
*/
	for (int i = 0; i < Nx; i++) {
		RXc[i] = (real_t)(cdt / (Xn[i + 1] - Xn[i]));
	}
	for (int j = 0; j < Ny; j++) {
		RYc[j] = (real_t)(cdt / (Yn[j + 1] - Yn[j]));
	}
	for (int k = 0; k < Nz; k++) {
		RZc[k] = (real_t)(cdt / (Zn[k + 1] - Zn[k]));
	}
/*
	// debug
	for (int i = 0; i < Nx; i++) {
		printf("RXc[%d]=%.5f\n", i, RXc[i]);
	}
	for (int j = 0; j < Ny; j++) {
		printf("RYc[%d]=%.5f\n", j, RYc[j]);
	}
	for (int k = 0; k < Nz; k++) {
		printf("RZc[%d]=%.5f\n", k, RZc[k]);
	}
*/
}


// Mur factor
double factorMur(double d, id_t m)
{
	if (m != PEC) {
		const double vdt = (C * Dt) / sqrt(Material[m].epsr * Material[m].amur);
		return (vdt - d) / (vdt + d);
	}
	else {
		return -1;
	}
}


// setup
void setup(void)
{
	// time step (default value)
	if (Dt < 1e-20) {
		setup_timestep();
	}

	// pulse width (default value)
	if (Tw < 1e-20) {
		setup_pulsewidth();
	}

	// material ID
	setupId();
	//debugId2();
	if (VECTOR) {
		setupId_vector();
		setupId_surface_vector();
	}
	else {
		setupId_surface();
	}
	//debugId2();
	//debugId();

	// material factor
	setup_material();
	if (VECTOR) {
		setup_material_vector();
	}

	// dispersion
	setupDispersion();

	// mesh factor
	setup_mesh();

	// ABC arrays
	if      (iABC == 0) {
		setupMurHx(1);
		setupMurHy(1);
		setupMurHz(1);
	}
	else if (iABC == 1) {
		setupPmlEx(1);
		setupPmlEy(1);
		setupPmlEz(1);
		setupPmlHx(1);
		setupPmlHy(1);
		setupPmlHz(1);
		setupPml();
	}

	// setup DFT factor
	setupDft();
}
