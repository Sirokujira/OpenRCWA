/*
plot2dFar0d0.c (OpenMOM/OpenRCWA/OpenTHFD)

plot far0d F-char (2D)
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "ev.h"

void plot2dFar0d0(
	int nfreq, double *ef[7],
	int scaleuser, double scalemin, double scalemax, int scalediv,
	int fdiv, double fmin, double fmax,
	const char title[], const char comment[], const char unit[],
	double width, double height, double font)
{
	char str1[BUFSIZ], str2[BUFSIZ];
	const unsigned char rgb[][3] = {{0, 0, 0}, {255, 0, 0}, {0, 0, 255}};
	const char *ecomp[][3] = {{"E", "E-theta", "E-phi"}, {"E", "E-major", "E-minor"}, {"E", "E-RHCP", "E-LHCP"}};
	const double eps = 1e-10;

	// layout
	const double x1 = 0.12 * width;
	const double x2 = 0.90 * width;
	const double y1 = 0.10 * height;
	const double y2 = 0.90 * height;
	const double h = font;

	// alloc
	double **f = (double **)malloc(7 * sizeof(double *));
	for (int comp = 0; comp < 7; comp++) {
		f[comp] = (double *)malloc(nfreq * sizeof(double));
	}

	// to dB
	for (int comp = 0; comp < 7; comp++) {
		for (int ifreq = 0; ifreq < nfreq; ifreq++) {
			const double e = (ef[comp][ifreq] > eps) ? ef[comp][ifreq] : eps;
			f[comp][ifreq] = 20 * log10(e);
		}
	}

	// data : max
	double dmax = f[0][0];
	for (int ifreq = 1; ifreq < nfreq; ifreq++) {
		const double d = f[0][ifreq];
		if (d > dmax) dmax = d;
	}

	// scale : min, max, div
	double ymin, ymax;
	int ydiv;
	if (scaleuser == 0) {
		const double dbmax = +100;
		const double dbunit = 5;
		ymax = dbunit * (ceil)(dmax / dbunit);
		if (ymax > dbmax) ymax = dbmax;
		ydiv = 4;
		ymin = ymax - dbunit * ydiv;
	}
	else {
		ymin = scalemin;
		ymax = scalemax;
		ydiv = scalediv;
	}

	// plot (3 pages)
	for (int page = 0; page < 3; page++) {
		// new page
		ev2d_newPage();

		// grid
		ev2dlib_grid(x1, y1, x2, y2, fdiv, ydiv);

		// 3 components
		for (int comp = 0; comp < 3; comp++) {
			ev2d_setColor(rgb[comp][0], rgb[comp][1], rgb[comp][2]);

			// function
			const int icomp = (comp == 0) ? 0 : (comp == 1) ? (2 * page + 1) : (2 * page + 2);
			ev2dlib_func1(nfreq - 1, f[icomp], ymin, ymax, x1, y1, x2, y2);

			// component
			const double dx = (comp == 0) ? 11 : (comp == 1) ? 9 : 3;
			ev2d_drawString(x2 - dx * h, y2 + 1.8 * h, h, ecomp[page][comp]);
		}
		ev2d_setColor(0, 0, 0);

		// max, min
		sprintf(str1, "max=%.3f%s", dmax, unit);
		ev2d_drawString(x2 - 10.0 * h, y2 + 0.4 * h, h, str1);

		// X-axis
		sprintf(str1, "%.3e", fmin);
		sprintf(str2, "%.3e", fmax);
		ev2dlib_Xaxis(x1, x2, y1, h, str1, str2, "frequency[Hz]");

		// Y-axis
		sprintf(str1, "%.0f", ymin);
		sprintf(str2, "%.0f", ymax);
		ev2dlib_Yaxis(y1, y2, x1, h, str1, str2, unit);

		// title
		ev2d_drawString(x1, y2 + 1.8 * h, h, comment);
		ev2d_drawString(x1, y2 + 0.4 * h, h, title);
	}

	// free
	for (int comp = 0; comp < 7; comp++) {
		free(f[comp]);
	}
	free(f);
}
