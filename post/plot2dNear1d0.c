/*
plot2dNear1d0.c (OpenMOM/OpenRCWA/OpenTHFD)

plot a near field the line (2D)
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "ev.h"

void plot2dNear1d0(
	const char comp[], int div, double *amp[3], double *deg[3], const double lng[],
	int scaledb, int scaleuser, double scalemin, double scalemax, int scalediv,
	const char title[], double frequency, double pos[2][3],
	double width, double height, double font)
{
	const int dbdiv = 5;
	const int dbauto = 40;
	const unsigned char rgbamp[][3] = {
		{  0,   0,   0},
		{255,   0,   0},
		{  0, 255,   0},
		{  0,   0, 255}};
	const unsigned char rgbdeg[][3] = {
		{255,   0,   0}};
	const char scmp[][2] = {"", "x", "y", "z"};
	const double eps = 1e-10;
	char str[3 * BUFSIZ], str1[BUFSIZ], str2[BUFSIZ], stru[BUFSIZ];

	// layout
	const double x1 = 0.15 * width;
	const double x2 = 0.92 * width;
	const double y1 = 0.12 * height;
	const double y2 = 0.90 * height;
	const double h = font;

	// alloc
	int namp = (!strcmp(comp, "E") || !strcmp(comp, "H") || !strcmp(comp, "I")) ? 4 : 1;
	int ndeg = (!strcmp(comp, "E") || !strcmp(comp, "H") || !strcmp(comp, "I")) ? 0 : 1;
	double **damp = (double **)malloc(namp * sizeof(double *));
	for (int n = 0; n < namp; n++) {
		damp[n] = (double *)malloc((div + 1) * sizeof(double));
	}
	double **ddeg = (double **)malloc(ndeg * sizeof(double *));
	for (int n = 0; n < ndeg; n++) {
		ddeg[n] = (double *)malloc((div + 1) * sizeof(double));
	}

	// field
	for (int i = 0; i <= div; i++) {
		const int c = toupper(comp[1]);
		if      (c == 'X') {
			damp[0][i] = amp[0][i];
			ddeg[0][i] = deg[0][i];
		}
		else if (c == 'Y') {
			damp[0][i] = amp[1][i];
			ddeg[0][i] = deg[1][i];
		}
		else if (c == 'Z') {
			damp[0][i] = amp[2][i];
			ddeg[0][i] = deg[2][i];
		}
		else {
			damp[0][i] = sqrt((amp[0][i] * amp[0][i])
			                + (amp[1][i] * amp[1][i])
			                + (amp[2][i] * amp[2][i]));
			damp[1][i] = amp[0][i];
			damp[2][i] = amp[1][i];
			damp[3][i] = amp[2][i];
		}
	}

	// max
	double dmax = damp[0][0];
	for (int i = 1; i <= div; i++) {
		if (damp[0][i] > dmax) dmax = damp[0][i];
	}
	if (dmax < eps) dmax = eps;

	// min, max, division
	double ymin = 0, ymax = 0;
	int    ydiv = 0;
	if (scaleuser) {
		// user scale
		ymin = scalemin;
		ymax = scalemax;
		ydiv = scalediv;
	}
	else {
		// auto scale
		if (scaledb) {
			// dB
			double dbmax = 20 * log10(dmax);
			ydiv = dbauto / dbdiv;
			ymax = dbdiv * ceil(dbmax / dbdiv);
			ymin = ymax - dbauto;
		}
		else {
			// linear
			double exp = pow(10, floor(log10(dmax)));
			ydiv = (int)((ceil)(dmax / exp));
			ymax = ydiv * exp;
			ymin = 0;
		}
	}

	// to dB
	if (scaledb) {
		for (int n = 0; n < namp; n++) {
			for (int i = 0; i <= div; i++) {
				damp[n][i] = 20 * log10((damp[n][i] > eps) ? damp[n][i] : eps);
			}
		}
	}

	// new page
	ev2d_newPage();

	// grid
	ev2dlib_grid(x1, y1, x2, y2, 10, ydiv);

	// plot amplitude
	for (int n = 0; n < namp; n++) {
		ev2d_setColorA(rgbamp[n]);
		ev2dlib_func2(div, lng, damp[n], ymin, ymax, x1, y1, x2, y2);
	}

	// plot phase
	for (int n = 0; n < ndeg; n++) {
		ev2d_setColorA(rgbdeg[n]);
		ev2dlib_func2(div, lng, ddeg[n], -180, +180, x1, y1, x2, y2);
	}

	ev2d_setColor(0, 0, 0);

	// X-axis
	for (int k = 0; k < 2; k++) {
		double x = (k == 0) ? (x1 - 3.0 * h) : (x2 - 7.0 * h);
		sprintf(str, "X=%.3e[m]", pos[k][0]);
		ev2d_drawString(x, y1 - 1.2 * h, 0.8 * h, str);
		sprintf(str, "Y=%.3e[m]", pos[k][1]);
		ev2d_drawString(x, y1 - 2.2 * h, 0.8 * h, str);
		sprintf(str, "Z=%.3e[m]", pos[k][2]);
		ev2d_drawString(x, y1 - 3.2 * h, 0.8 * h, str);
	}

	// Y-axis
	sprintf(stru, "[%s%s]", (scaledb ? "dB" : ""), ((comp[0] == 'E') ? "V/m" : (comp[0] == 'H') ? "A/m" : "A"));
	sprintf(str1, (scaledb ? "%.1f" : "%.3e"), ymin);
	sprintf(str2, (scaledb ? "%.1f" : "%.3e"), ymax);
	ev2dlib_Yaxis(y1, y2, x1, 0.8 * h, str1, str2, stru);

	if (ndeg > 0) {
		ev2d_setColorA(rgbdeg[0]);
		for (int m = 1; m < 12; m++) {
			double y = y1 + (m / 12.0) * (y2 - y1);
			ev2d_drawLine(x2 - 0.5 * h, y, x2, y);
		}
		double x = x2 + 0.5 * h;
		ev2d_drawString(x, y1            - 0.3 * h, 0.8 * h, "-180");
		ev2d_drawString(x, (y1 + y2) / 2 - 0.3 * h, 0.8 * h, "0");
		ev2d_drawString(x, y2            - 0.3 * h, 0.8 * h, "+180");
		ev2d_drawString(x, y2            - 1.5 * h, 0.8 * h, "[deg]");
		ev2d_setColor(0, 0, 0);
	}

	// title
	ev2d_drawString(x1, y2 + 0.4 * h, h, title);

	// comment
	sprintf(str1, (scaledb ? "%.3f" : "%.3e"), (scaledb ? 20 * log10(dmax) : dmax));
	sprintf(str, "f=%.3e[Hz] max=%s%s", frequency, str1, stru);
	ev2d_drawString(x1, y2 + 1.8 * h, h, str);

	// component
	if (!strcmp(comp, "E") || !strcmp(comp, "H") || !strcmp(comp, "I")) {
		for (int n = 0; n < namp; n++) {
			sprintf(str, "%c%s", comp[0], scmp[n]);
			ev2d_setColorA(rgbamp[n]);
			double x = x2 - 8.0 * h + (n * 2 * h);
			ev2d_drawString(x, y2 + 0.4 * h, h, str);
		}
		ev2d_setColor(0, 0, 0);
	}
	else {
		ev2d_drawString(x2 - 8.0 * h, y2 + 0.4 * h, h, comp);
	}

	// free
	for (int n = 0; n < namp; n++) {
		free(damp[n]);
	}
	for (int n = 0; n < ndeg; n++) {
		free(ddeg[n]);
	}
}
