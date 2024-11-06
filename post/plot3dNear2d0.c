/*
plot3dNear2d0.c (OpenMOM/OpenRCWA/OpenTHFD)

plot near field on a plane (3D)
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "ev.h"

void plot3dNear2d0(
	int div1, int div2, double **d, const double pos1[], const double pos2[],
	int phase, int scaledb, int scaleuser, double scalemin, double scalemax, int contour,
	const char title[], double frequency,
	char dir, double pos0, const char component[],
	int nwire, double (*wire)[2][3])
{
	// new page
	ev3d_newPage();

	if ((div1 <= 0) || (div2 <= 0)) return;
	if (phase) scaledb = 0;

	const double dbauto = 40;

	// max
	double dmax = 0;
	if (!phase) {
		dmax = d[0][0];
		for (int i = 0; i <= div1; i++) {
		for (int j = 0; j <= div2; j++) {
			if (d[i][j] > dmax) dmax = d[i][j];
		}
		}
	}

	// min, max, division
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

	// contour
	for (int i = 0; i < div1; i++) {
	for (int j = 0; j < div2; j++) {
		double p1 = 0, q1 = 0, p2 = 0, q2 = 0;
		if      (dir == 'X') {
			p1 = pos1[i];
			q1 = pos2[j];
			p2 = pos1[i + 1];
			q2 = pos2[j + 1];
		}
		else if (dir == 'Y') {
			q1 = pos1[i];
			p1 = pos2[j];
			q2 = pos1[i + 1];
			p2 = pos2[j + 1];
		}
		else if (dir == 'Z') {
			p1 = pos1[i];
			q1 = pos2[j];
			p2 = pos1[i + 1];
			q2 = pos2[j + 1];
		}
		double v = (d[i][j] + d[i + 1][j] + d[i][j + 1] + d[i + 1][j + 1]) / 4;
		v = (v - fmin) / (fmax - fmin);
		ev3d_setColorV(v, (contour < 2));
		ev3d_fillRectangle(dir, pos0, p1, q1, p2, q2);
	}
	}

	// objects
	ev3d_setColor(0, 0, 0);
	ev3d_index(1);
	for (int i = 0; i < nwire; i++) {
		ev3d_drawLine(
			wire[i][0][0], wire[i][0][1], wire[i][0][2],
			wire[i][1][0], wire[i][1][1], wire[i][1][2]);
	}

	// title
	char str[4 * BUFSIZ], strmax[BUFSIZ], strunit[BUFSIZ];
	sprintf(strunit, "[%s%s]", (scaledb ? "dB" : ""), (phase ? "deg" : (component[0] == 'E' ? "V/m" : "A/m")));
	if (!phase) {
		sprintf(strmax, (scaledb ? "max=%.3f%s" : "max=%.3e%s"), dmax, strunit);
	}
	else{
		strcpy(strmax, "");
	}
	sprintf(str, "%s%s %c=%.3e[m] f=%.3e[Hz] %s", component, strunit, dir, pos0, frequency, strmax);
	ev3d_setColor(0, 0, 0);
	ev3d_drawTitle(str);
	ev3d_drawTitle(title);
}
