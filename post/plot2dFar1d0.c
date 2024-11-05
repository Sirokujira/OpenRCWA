/*
plot2dFar1d0.c (OpenMOM/OpenRCWA/OpenTHFD)

plot far-1d (2D)
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "ev.h"

void plot2dFar1d0(
	int div, double (*ef)[7],
	int component[3], char dir, double angle, int figurestyle,
	int scaledb, int scaleuser, double scalemin, double scalemax, int scalediv,
	const char title[], const char strunit[], double frequency,
	double width, double height, double font)
{
	const double dbdiv = 5;  // dB division
	const char estr[][2][BUFSIZ] = {{"E-theta", "E-phi"}, {"E-major", "E-minor"}, {"E-RHCP", "E-LHCP"}};
	const double eps = 1e-15;
	const double pi = 4 * atan(1);
	char str[8][BUFSIZ], fmt[BUFSIZ];

	// layout
	const double x1 = 0.12 * width;
	const double x2 = 0.92 * width;
	const double y1 = 0.10 * height;
	const double y2 = 0.90 * height;
	const double h = font;

	const double x0 = (x1 + x2) / 2;
	const double y0 = (y1 + y2) / 2;
	const double r0 = fabs(y2 - y1) / 2;

	// max
	double emax = 0;
	for (int i = 0; i <= div; i++) {
		double e = ef[i][0];
		if (e > emax) emax = e;
	}
	if (emax < eps) emax = eps;

	// min, max, division
	double fmin = 0, fmax = 0;
	int    fdiv = 0;
	if (scaleuser) {
		// user scale
		fmin = scalemin;
		fmax = scalemax;
		fdiv = scalediv;
	}
	else {
		// auto scale
		if (scaledb) {
			// dB
			double dbmax = 20 * log10(emax);
			fdiv = 8;
			fmax = dbdiv * ceil(dbmax / dbdiv);
			fmin = fmax - (dbdiv * fdiv);
		}
		else {
			// linear
			double exp = pow(10, floor(log10(emax)));
			fdiv = (int)((ceil)(emax / exp));
			fmax = fdiv * exp;
			fmin = 0;
		}
	}
	//printf("%f %f %d\n", fmin, fmax, fdiv);

	// alloc
	double *ef0 = (double *)malloc((div + 1) * sizeof(double));
	double *ef1 = (double *)malloc((div + 1) * sizeof(double));
	double *ef2 = (double *)malloc((div + 1) * sizeof(double));

	// plot parameter
	double p0 = 0;
	double sgn = 0;
	if      (dir == 'X') {
		strcpy(str[0], "+Y");
		strcpy(str[1], "+Z");
		strcpy(str[2], " -Y");
		strcpy(str[3], "-Z");
		strcpy(str[4], "theta[deg]  (phi[deg]=90)");
		p0 = pi / 2;
		sgn = -1;
	}
	else if (dir == 'Y') {
		strcpy(str[0], "+X");
		strcpy(str[1], "+Z");
		strcpy(str[2], " -X");
		strcpy(str[3], "-Z");
		strcpy(str[4], "theta[deg]  (phi[deg]=0)");
		p0 = pi / 2;
		sgn = -1;
	}
	else if (dir == 'Z') {
		strcpy(str[0], "+X");
		strcpy(str[1], "+Y");
		strcpy(str[2], " -X");
		strcpy(str[3], "-Y");
		strcpy(str[4], "phi[deg]  (theta[deg]=90)");
		p0 = 0;
		sgn = +1;
	}
	else if (dir == 'V') {
		strcpy(str[0], "90");
		strcpy(str[1], "0");
		strcpy(str[2], "270");
		strcpy(str[3], "180");
		sprintf(str[4], "theta[deg]  (phi[deg]=%.2f)", angle);
		p0 = pi / 2;
		sgn = -1;
	}
	else if (dir == 'H') {
		strcpy(str[0], "0");
		strcpy(str[1], "90");
		strcpy(str[2], "180");
		strcpy(str[3], "270");
		sprintf(str[4], "phi[deg]  (theta[deg]=%.2f)", angle);
		p0 = 0;
		sgn = +1;
	}

	for (int icomponent = 0; icomponent < 3; icomponent++) {
		if (!component[icomponent]) continue;

		// E field
		for (int i = 0; i <= div; i++) {
			ef0[i] = ef[i][0];
			if      (icomponent == 0) {
				ef1[i] = ef[i][1];
				ef2[i] = ef[i][2];
			}
			else if (icomponent == 1) {
				ef1[i] = ef[i][3];
				ef2[i] = ef[i][4];
			}
			else if (icomponent == 2) {
				ef1[i] = ef[i][5];
				ef2[i] = ef[i][6];
			}

			// to dB
			if (scaledb) {
				ef0[i] = 20 * log10((ef0[i] > eps) ? ef0[i] : eps);
				ef1[i] = 20 * log10((ef1[i] > eps) ? ef1[i] : eps);
				ef2[i] = 20 * log10((ef2[i] > eps) ? ef2[i] : eps);
			}
		}

		// new page
		ev2d_newPage();

		// grid
		if      (figurestyle == 0) {
			ev2dlib_CircleMesh(x1, y1, x2, y2, fdiv, 18);
		}
		else if (figurestyle == 1) {
			ev2dlib_grid(x1, y1, x2, y2, 12, fdiv);
		}

		// (1) abs (black)
		ev2d_setColor(  0,   0,   0);
		if      (figurestyle == 0) {
			ev2dlib_CircleFunc(x1, y1, x2, y2, ef0, div, fmin, fmax, p0, sgn);
		}
		else if (figurestyle == 1) {
			ev2dlib_func1(div, ef0, fmin, fmax, x1, y1, x2, y2);
		}
		ev2d_drawString(x1 + 24.0 * h, y2 + 1.8 * h, h, "E-abs");

		// (2) component-1 (red)
		ev2d_setColor(255,   0,   0);
		if      (figurestyle == 0) {
			ev2dlib_CircleFunc(x1, y1, x2, y2, ef1, div, fmin, fmax, p0, sgn);
		}
		else if (figurestyle == 1) {
			ev2dlib_func1(div, ef1, fmin, fmax, x1, y1, x2, y2);
		}
		ev2d_drawString(x1 + 29.0 * h, y2 + 1.8 * h, h, estr[icomponent][0]);

		// (3) component-2 (blue)
		ev2d_setColor(  0,   0, 255);
		if      (figurestyle == 0) {
			ev2dlib_CircleFunc(x1, y1, x2, y2, ef2, div, fmin, fmax, p0, sgn);
		}
		else if (figurestyle == 1) {
			ev2dlib_func1(div, ef2, fmin, fmax, x1, y1, x2, y2);
		}
		ev2d_drawString(x1 + 35.0 * h, y2 + 1.8 * h, h, estr[icomponent][1]);

		// angle
		ev2d_setColor(  0,   0,   0);
		if      (figurestyle == 0) {
			ev2d_drawString(x0 + r0 + 0.2 * h, y0      - 0.2 * h, h, str[0]);
			ev2d_drawString(x0      - 0.5 * h, y0 + r0 + 0.2 * h, h, str[1]);
			ev2d_drawString(x0 - r0 - 2.0 * h, y0      - 0.2 * h, h, str[2]);
			ev2d_drawString(x0      - 0.5 * h, y0 - r0 - 1.0 * h, h, str[3]);
			// constant angle
			if ((dir == 'V') || (dir == 'H')) {
				ev2d_drawString(x0 + 4.0 * h, y2 + 0.6 * h, h, str[4]);
			}
		}
		else if (figurestyle == 1) {
			ev2dlib_Xaxis(x1, x2, y1, h, "0", "360", str[4]);
		}

		// level
		ev2d_setColor(  0,   0,   0);
		if      (figurestyle == 0) {
			strcpy(fmt,  (scaledb ? "%.0f%s" : "%.1e%s"));
			sprintf(str[5], "%.0f",   fmin);
			sprintf(str[6], fmt, fmax, strunit);
			ev2d_drawString(x0      - 0.8 * h, y0 - 1.0 * h, h, str[5]);
			ev2d_drawString(x0 + r0 - 0.8 * h, y0 - 1.0 * h, h, str[6]);
		}
		else if (figurestyle == 1) {
			strcpy(fmt, (scaledb ? "%.0f" : "%.1e"));
			sprintf(str[5], "%.0f", fmin);
			sprintf(str[6], fmt, fmax);
			ev2dlib_Yaxis(y1, y2, x1, h, str[5], str[6], strunit);
		}

		// title
		ev2d_drawString(x1, y2 + 0.4 * h, h, title);

		// comment
		double pmax = scaledb ? 20 * log10(emax) : emax;
		sprintf(fmt, "f=%%.3e[Hz] max=%s%%s", (scaledb ? "%.3f" : "%.3e"));
		sprintf(str[7], fmt, frequency, pmax, strunit);
		ev2d_drawString(x1, y2 + 1.8 * h, h, str[7]);
	}

	// free
	free(ef0);
	free(ef1);
	free(ef2);
}
