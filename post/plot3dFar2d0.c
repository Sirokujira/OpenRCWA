/*
plot3dFar2d0.c (OpenMOM/OpenRCWA/OpenTHFD)

plot a far2d pattern (3D)
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "ev.h"

static void plot3dwire(int, double (*)[2][3], double);

void plot3dFar2d0(
	int divtheta, int divphi, double **d,
	int scaledb, int scaleuser, double scalemin, double scalemax,
	int nwire, double (*wire)[2][3], double rscale,
	int ncomment, char **comment)
{
	if ((divtheta <= 0) || (divphi <= 0)) return;

	const double dbauto = 40;
	const double eps = 1e-10;

	// to dB
	if (scaledb) {
		for (int itheta = 0; itheta <= divtheta; itheta++) {
		for (int iphi   = 0; iphi   <= divphi;   iphi++  ) {
			d[itheta][iphi] = 20 * log10((d[itheta][iphi] > eps) ? d[itheta][iphi] : eps);
		}
		}
	}

	// max
	double dmax = d[0][0];
	for (int itheta = 0; itheta <= divtheta; itheta++) {
	for (int iphi   = 0; iphi   <= divphi;   iphi++  ) {
		if (d[itheta][iphi] > dmax) {
			dmax = d[itheta][iphi];
		}
	}
	}

	// min, max
	double fmin = 0, fmax = 0;
	if (scaleuser) {
		// user scale
		fmin = scalemin;
		fmax = scalemax;
	}
	else {
		// auto scale
		fmax = dmax;
		if (scaledb) {
			// dB
			fmin = fmax - dbauto;
		}
		else {
			// linear
			fmin = 0;
		}
	}

	// scaled data
	for (int itheta = 0; itheta <= divtheta; itheta++) {
	for (int iphi   = 0; iphi   <= divphi;   iphi++  ) {
		double f = d[itheta][iphi];
		if (f < fmin) f = fmin;
		d[itheta][iphi] = (f - fmin) / (fmax - fmin);
	}
	}

	// new page
	ev3d_newPage();

	// objects
	if (nwire && (rscale > 0)) {
		plot3dwire(nwire, wire, rscale);
	}

	// plot
	ev3dlib_func(divtheta, 0, 180, divphi, 0, 360, d);

	// comment
	ev3d_setColor(0, 0, 0);
	for (int n = 0; n < ncomment; n++) {
		ev3d_drawTitle(comment[n]);
	}
}


static void plot3dwire(int nwire, double (*wire)[2][3], double rscale)
{
	double gmin[3], gmax[3];
	for (int m = 0; m < 3; m++) {
		gmin[m] = +1e10;
		gmax[m] = -1e10;
	}
	for (int n = 0; n < nwire; n++) {
		for (int v = 0; v < 2; v++) {
			for (int m = 0; m < 3; m++) {
				const double g = wire[n][v][m];
				if (g < gmin[m]) gmin[m] = g;
				if (g > gmax[m]) gmax[m] = g;
			}
		}
	}

	const double dx = gmax[0] - gmin[0];
	const double dy = gmax[1] - gmin[1];
	const double dz = gmax[2] - gmin[2];
	const double rr = sqrt((dx * dx) + (dy * dy) + (dz * dz)) / 2;
	const double x0 = (gmin[0] + gmax[0]) / 2;
	const double y0 = (gmin[1] + gmax[1]) / 2;
	const double z0 = (gmin[2] + gmax[2]) / 2;
	const double rf = rscale / rr;
	ev3d_setColor(0, 0, 0);
	for (int n = 0; n < nwire; n++) {
		const double x1 = rf * (wire[n][0][0] - x0);
		const double y1 = rf * (wire[n][0][1] - y0);
		const double z1 = rf * (wire[n][0][2] - z0);
		const double x2 = rf * (wire[n][1][0] - x0);
		const double y2 = rf * (wire[n][1][1] - y0);
		const double z2 = rf * (wire[n][1][2] - z0);
		ev3d_drawLine(x1, y1, z1, x2, y2, z2);
	}
}
