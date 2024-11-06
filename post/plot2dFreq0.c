/*
plot2dFreq0.c (OpenMOM/OpenRCWA/OpenTHFD)

plot frequency char.s (2D)
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

typedef struct {double r, i;} d_complex_t;
typedef struct {float  r, i;} f_complex_t;

#include "complex.h"
#include "ev.h"

typedef struct {
	int    db;
	int    user;
	double min, max;
	int    div;
} scale_t;

#define MAX(x, y) ((x) > (y) ? (x) : (y))
#define MIN(x, y) ((x) < (y) ? (x) : (y))

static void smithChart(int, double *, double *, double, double, double, double);

// smith chart
void plot2dSmith(
	int nfeed, int nfreq, const d_complex_t zin[], const double z0[], const double freq[], const char title[],
	double width, double height, double font)
{
	char str[BUFSIZ];

	// layout
	const double r0 = 0.47 * MIN(width, height);
	const double y0 = 0.5 * height;
	const double x0 = width - 0.53 * MIN(width, height);
	const double h = font;

	// alloc
	double *rin = (double *)malloc(nfreq * sizeof(double));
	double *xin = (double *)malloc(nfreq * sizeof(double));

	for (int ifeed = 0; ifeed < nfeed; ifeed++) {

		ev2d_newPage();

		for (int ifreq = 0; ifreq < nfreq; ifreq++) {
			d_complex_t z = zin[(ifeed * nfreq) + ifreq];
			rin[ifreq] = z.r;
			xin[ifreq] = z.i;
		}

		smithChart(nfreq, rin, xin, z0[ifeed], x0, y0, r0);

		// comment
		double x = 1.0 * h;
		ev2d_drawString(x, height - 2.0 * h, h, "SMITH CHART");

		ev2d_drawString(x, height - 3.5 * h, h, title);

		ev2d_fillRectangle(x, height - 5.0 * h, x + 0.5 * h, height - 4.5 * h);
		sprintf(str, "start=%.3e[Hz]", freq[0]);
		ev2d_drawString(x + 0.8 * h, height - 5.0 * h, h, str);

		ev2d_drawRectangle(x, height - 6.5 * h, x + 0.5 * h, height - 6.0 * h);
		sprintf(str, "stop=%.3e[Hz]", freq[nfreq - 1]);
		ev2d_drawString(x + 0.8 * h, height - 6.5 * h, h, str);

		sprintf(str, "Z0=%.1f[ohm]", z0[ifeed]);
		ev2d_drawString(x, height - 8.0 * h, h, str);

		if (nfeed > 1) {
			sprintf(str, "feed #%d", ifeed + 1);
			ev2d_drawString(x, height - 9.5 * h, h, str);
		}
	}

	// free
	free(rin);
	free(xin);
}


// Zin
void plot2dZin(
	int nfeed, int nfreq, const d_complex_t zin[],
	scale_t scale, int fdiv, const double freq[], const char title[],
	double width, double height, double font)
{
	char str1[BUFSIZ], str2[BUFSIZ];
	const unsigned char rgb[][3] = {{255,   0,   0},
	                                {  0,   0, 255}};

	// layout
	const double x1 = 0.12 * width;
	const double x2 = 0.90 * width;
	const double y1 = 0.10 * height;
	const double y2 = 0.90 * height;
	const double h = font;

	// alloc
	double *f = (double *)malloc(nfreq * sizeof(double));

	// data : max
	double dmax = 0;
	for (int ifeed = 0; ifeed < nfeed; ifeed++) {
		for (int ifreq = 0; ifreq < nfreq; ifreq++) {
			d_complex_t z = zin[(ifeed * nfreq) + ifreq];
			dmax = MAX(dmax, MAX(fabs(z.r), fabs(z.i)));
		}
	}

	// scale : min, max, div
	double ymin, ymax;
	int ydiv;
	if (scale.user == 0) {
		const double unit = 100;
		ydiv = MIN((int)(dmax / unit) + 1, 50);
		ymax = unit * ydiv;
		ymin = -ymax;
		ydiv *= 2;
	}
	else {
		ymin = scale.min;
		ymax = scale.max;
		ydiv = scale.div;
	}

	// plot
	ev2d_newPage();

	// grid
	ev2dlib_grid(x1, y1, x2, y2, fdiv, ydiv);

	ev2d_setColor(0, 0, 0);

	for (int ifeed = 0; ifeed < nfeed; ifeed++) {
		for (int ic = 0; ic < 2; ic++) {
			for (int ifreq = 0; ifreq < nfreq; ifreq++) {
				f[ifreq] = (ic == 0) ? zin[(ifeed * nfreq) + ifreq].r
				                     : zin[(ifeed * nfreq) + ifreq].i;
			}

			ev2d_setColor(rgb[ic][0], rgb[ic][1], rgb[ic][2]);
			ev2dlib_func1(nfreq - 1, f, ymin, ymax, x1, y1, x2, y2);
			if (nfeed > 1) {
				sprintf(str1, "#%d", ifeed + 1);
				const double y = y1 + (y2 - y1) * (f[nfreq - 1] - ymin) / (ymax - ymin);
				ev2d_drawString(x2, y - 0.4 * h, h, str1);
			}
		}
	}

	ev2d_setColor(0, 0, 0);

	// X-axis
	sprintf(str1, "%.3e", freq[0]);
	sprintf(str2, "%.3e", freq[nfreq - 1]);
	ev2dlib_Xaxis(x1, x2, y1, h, str1, str2, "frequency[Hz]");

	// Y-axis
	sprintf(str1, "%.0f", ymin);
	sprintf(str2, "%.0f", ymax);
	ev2dlib_Yaxis(y1, y2, x1, h, str1, str2, "[ohm]");

	// Y=0
	const double y0 = y1 + (y2 - y1) * (0 - ymin) / (ymax - ymin);
	if ((y0 > y1) && (y0 < y2)) {
		ev2d_drawLine(x1, y0, x2, y0);
		ev2d_drawString(x1 - 1.5 * h, y0 - 0.3 * h, h, "0");
	}

	// title
	ev2d_drawString(x1, y2 + 0.4 * h, h, title);
	ev2d_drawString(x1, y2 + 1.8 * h, h, "input impedance");
	ev2d_setColor(rgb[0][0], rgb[0][1], rgb[0][2]);
	ev2d_drawString(x1 + 16.0 * h, y2 + 1.8 * h, h, "Rin");
	ev2d_setColor(rgb[1][0], rgb[1][1], rgb[1][2]);
	ev2d_drawString(x1 + 19.0 * h, y2 + 1.8 * h, h, "Xin");

	// free
	free(f);
}


// Yin
void plot2dYin(
	int nfeed, int nfreq, const d_complex_t zin[],
	scale_t scale, int fdiv, const double freq[], const char title[],
	double width, double height, double font)
{
	char str1[BUFSIZ], str2[BUFSIZ];
	const unsigned char rgb[][3] = {{255,   0,   0},
	                                {  0,   0, 255}};

	// layout
	const double x1 = 0.12 * width;
	const double x2 = 0.90 * width;
	const double y1 = 0.10 * height;
	const double y2 = 0.90 * height;
	const double h = font;

	// alloc
	double *f = (double *)malloc(nfreq * sizeof(double));

	// data : max
	double dmax = 0;
	for (int ifeed = 0; ifeed < nfeed; ifeed++) {
		for (int ifreq = 0; ifreq < nfreq; ifreq++) {
			const d_complex_t z = d_inv(zin[(ifeed * nfreq) + ifreq]);
			dmax = MAX(dmax, MAX(fabs(z.r), fabs(z.i)));
		}
	}

	// scale : min, max, div
	double ymin, ymax;
	int ydiv;
	if (scale.user == 0) {
		const double unit = 10e-3;
		ydiv = MIN((int)(dmax / unit) + 1, 50);
		ymax = unit * ydiv;
		ymin = -ymax;
		ydiv *= 2;
	}
	else {
		ymin = scale.min;
		ymax = scale.max;
		ydiv = scale.div;
	}

	// plot
	ev2d_newPage();

	// grid
	ev2dlib_grid(x1, y1, x2, y2, fdiv, ydiv);

	ev2d_setColor(0, 0, 0);

	for (int ifeed = 0; ifeed < nfeed; ifeed++) {
		for (int ic = 0; ic < 2; ic++) {
			for (int ifreq = 0; ifreq < nfreq; ifreq++) {
				f[ifreq] = (ic == 0) ? d_inv(zin[(ifeed * nfreq) + ifreq]).r
				                     : d_inv(zin[(ifeed * nfreq) + ifreq]).i;
			}

			ev2d_setColor(rgb[ic][0], rgb[ic][1], rgb[ic][2]);
			ev2dlib_func1(nfreq - 1, f, ymin, ymax, x1, y1, x2, y2);
			if (nfeed > 1) {
				sprintf(str1, "#%d", ifeed + 1);
				const double y = y1 + (y2 - y1) * (f[nfreq - 1] - ymin) / (ymax - ymin);
				ev2d_drawString(x2, y - 0.4 * h, h, str1);
			}
		}
	}

	ev2d_setColor(0, 0, 0);

	// X-axis
	sprintf(str1, "%.3e", freq[0]);
	sprintf(str2, "%.3e", freq[nfreq - 1]);
	ev2dlib_Xaxis(x1, x2, y1, h, str1, str2, "frequency[Hz]");

	// Y-axis
	sprintf(str1, "%.0f", ymin * 1e3);
	sprintf(str2, "%.0f", ymax * 1e3);
	ev2dlib_Yaxis(y1, y2, x1, h, str1, str2, "[mS]");

	// Y=0
	const double y0 = y1 + (y2 - y1) * (0 - ymin) / (ymax - ymin);
	if ((y0 > y1) && (y0 < y2)) {
		ev2d_drawLine(x1, y0, x2, y0);
		ev2d_drawString(x1 - 1.5 * h, y0 - 0.3 * h, h, "0");
	}

	// title
	ev2d_drawString(x1, y2 + 0.4 * h, h, title);
	ev2d_drawString(x1, y2 + 1.8 * h, h, "input admittance");
	ev2d_setColor(rgb[0][0], rgb[0][1], rgb[0][2]);
	ev2d_drawString(x1 + 16.0 * h, y2 + 1.8 * h, h, "Gin");
	ev2d_setColor(rgb[1][0], rgb[1][1], rgb[1][2]);
	ev2d_drawString(x1 + 19.0 * h, y2 + 1.8 * h, h, "Bin");

	// free
	free(f);
}


// reflection
void plot2dRef(
	int nfeed, int nfreq, const d_complex_t zin[], const double z0[],
	scale_t scale, int fdiv, const double freq[], const char title[],
	double width, double height, double font)
{
	const double vswr[] = {1.1, 1.2, 1.5, 2.0, 3.0, 5.0, 10.0};
	char str1[BUFSIZ], str2[BUFSIZ];

	// layout
	const double x1 = 0.12 * width;
	const double x2 = 0.90 * width;
	const double y1 = 0.10 * height;
	const double y2 = 0.90 * height;
	const double h = font;

	// alloc
	double *f = (double *)malloc(nfreq * sizeof(double));
	double *ref = (double *)malloc(nfeed * nfreq * sizeof(double));

	// reflection [dB]
	for (int ifeed = 0; ifeed < nfeed; ifeed++) {
		const d_complex_t z0_ = d_complex(z0[ifeed], 0);
		for (int ifreq = 0; ifreq < nfreq; ifreq++) {
			const d_complex_t zin_ = zin[(ifeed * nfreq) + ifreq];
			const double gamma = d_abs(d_div(d_sub(zin_, z0_), d_add(zin_, z0_)));
			ref[(ifeed * nfreq) + ifreq] = 20 * log10(MAX(gamma, 1e-6));
		}
	}

	// data : min
	double dmin = ref[0];
	for (int ifeed = 0; ifeed < nfeed; ifeed++) {
		for (int ifreq = 0; ifreq < nfreq; ifreq++) {
			dmin = MIN(dmin, ref[(ifeed * nfreq) + ifreq]);
		}
	}

	// scale : min, max, div
	double ymin, ymax;
	int ydiv;
	if (scale.user == 0) {
		const double dbmin = -80;
		const double unit = 5;
		ydiv = (int)(-dmin / unit) + 1;
		ymin = -unit * ydiv;
		ymin = MAX(ymin, dbmin);
		ymax = 0;
	}
	else {
		ymin = scale.min;
		ymax = scale.max;
		ydiv = scale.div;
	}

	// plot
	ev2d_newPage();

	// grid
	ev2dlib_grid(x1, y1, x2, y2, fdiv, ydiv);

	ev2d_setColor(0, 0, 0);

	for (int ifeed = 0; ifeed < nfeed; ifeed++) {
		for (int ifreq = 0; ifreq < nfreq; ifreq++) {
			f[ifreq] = ref[(ifeed * nfreq) + ifreq];
		}

		ev2dlib_func1(nfreq - 1, f, ymin, ymax, x1, y1, x2, y2);
		if (nfeed > 1) {
			sprintf(str1, "#%d", ifeed + 1);
			double y = y1 + (y2 - y1) * (f[nfreq - 1] - ymin) / (ymax - ymin);
			ev2d_drawString(x2, y - 0.4 * h, h, str1);
		}

		// Z0 : constant!
		sprintf(str1, "Z0=%.1f[ohm]", z0[ifeed]);
		ev2d_drawString(x2 - 10.0 * h, y2 + 0.6 * h, h, str1);
	}

	ev2d_setColor(0, 0, 0);

	// VSWR (right side)
	const int nvswr = sizeof(vswr) / sizeof(double);
	ev2d_drawString(x2, y2 + 0.6 * h, h, "VSWR");
	for (int i = 0; i < nvswr; i++) {
		double y = 20.0 * log10((vswr[i] - 1) / (vswr[i] + 1));
		y = y1 + (y2 - y1) * (y - ymin) / (ymax - ymin);
		if ((y >= y1) && (y <= y2)) {
			ev2d_drawLine(x2, y, x2 - 0.7 * h, y);
			sprintf(str1, "%g", vswr[i]);
			ev2d_drawString(x2 + 0.2 * h, y - 0.3 * h, h, str1);
		}
	}

	// X-axis
	sprintf(str1, "%.3e", freq[0]);
	sprintf(str2, "%.3e", freq[nfreq - 1]);
	ev2dlib_Xaxis(x1, x2, y1, h, str1, str2, "frequency[Hz]");

	// Y-axis
	sprintf(str1, "%.0f", ymin);
	sprintf(str2, "%.0f", ymax);
	ev2dlib_Yaxis(y1, y2, x1, h, str1, str2, "[dB]");

	// title
	ev2d_drawString(x1, y2 + 0.4 * h, h, title);
	ev2d_drawString(x1, y2 + 1.8 * h, h, "reflection");

	// min
	sprintf(str1, "min=%.3f[dB]", dmin);
	ev2d_drawString(x2 - 10.0 * h, y2 + 1.8 * h, h, str1);

	// free
	free(f);
	free(ref);
}


// S-parameter
void plot2dSpara(
	int npoint, int nfreq, const d_complex_t spara[],
	scale_t scale, int fdiv, const double freq[], const char title[],
	double width, double height, double font)
{
	char str1[BUFSIZ], str2[BUFSIZ];
	const double eps = 1e-10;

	// alloc
	double *f = (double *)malloc(nfreq * sizeof(double));

	// layout
	const double x1 = 0.12 * width;
	const double x2 = 0.90 * width;
	const double y1 = 0.10 * height;
	const double y2 = 0.90 * height;
	const double h = font;

	// data : min, max
	double dmin = +1000;
	double dmax = -1000;
	for (int ipoint = 0; ipoint < npoint; ipoint++) {
		for (int ifreq = 0; ifreq < nfreq; ifreq++) {
			const double d = 20 * log10(MAX(d_abs(spara[(ipoint * nfreq) + ifreq]), eps));
			dmin = MIN(dmin, d);
			dmax = MAX(dmax, d);
		}
	}

	// scale : min, max, div
	double ymin, ymax;
	int ydiv;
	if (scale.user == 0) {
		const double dbmin = -80;
		const double unit = 5;
		ydiv = (int)(-dmin / unit) + 1;
		ymin = -unit * ydiv;
		ymin = MAX(ymin, dbmin);
		ymax = 0;
	}
	else {
		ymin = scale.min;
		ymax = scale.max;
		ydiv = scale.div;
	}

	// plot
	ev2d_newPage();

	// grid
	ev2dlib_grid(x1, y1, x2, y2, fdiv, ydiv);

	ev2d_setColor(0, 0, 0);

	for (int ipoint = 0; ipoint < npoint; ipoint++) {
		for (int ifreq = 0; ifreq < nfreq; ifreq++) {
			f[ifreq] = 20 * log10(MAX(d_abs(spara[(ipoint * nfreq) + ifreq]), eps));
		}

		if      (ipoint == 0) {
			ev2d_setColor(255,   0,   0);
		}
		else if (ipoint == 1) {
			ev2d_setColor(  0,   0, 255);
		}
		else if (ipoint == 2) {
			ev2d_setColor(  0, 255,   0);
		}
		else {
			ev2d_setColor(rand() & 0xff, rand() & 0xff, rand() & 0xff); 
		}

		ev2dlib_func1(nfreq - 1, f, ymin, ymax, x1, y1, x2, y2);

		sprintf(str1, "S%d1", ipoint + 1);
		double y = y1 + (y2 - y1) * (f[nfreq - 1] - ymin) / (ymax - ymin);
		ev2d_drawString(x2, y - 0.4 * h, h, str1);
	}

	ev2d_setColor(0, 0, 0);

	// min, max
	sprintf(str1, "max=%.3f[dB]", dmax);
	sprintf(str2, "min=%.3f[dB]", dmin);
	ev2d_drawString(x2 - 10.0 * h, y2 + 1.8 * h, h, str1);
	ev2d_drawString(x2 - 10.0 * h, y2 + 0.4 * h, h, str2);

	// X-axis
	sprintf(str1, "%.3e", freq[0]);
	sprintf(str2, "%.3e", freq[nfreq - 1]);
	ev2dlib_Xaxis(x1, x2, y1, h, str1, str2, "frequency[Hz]");

	// Y-axis
	sprintf(str1, "%.0f", ymin);
	sprintf(str2, "%.0f", ymax);
	ev2dlib_Yaxis(y1, y2, x1, h, str1, str2, "[dB]");

	// title
	ev2d_drawString(x1, y2 + 0.4 * h, h, title);
	ev2d_drawString(x1, y2 + 1.8 * h, h, "S-parameter");

	// free
	free(f);
}


// Coupling
void plot2dCoupling(
	int nfeed, int npoint, int nfreq, double ***couple,
	scale_t scale, int fdiv, const double freq[], const char title[],
	double width, double height, double font)
{
	char str1[BUFSIZ], str2[BUFSIZ];
	const double eps = 1e-10;

	// alloc
	double *f = (double *)malloc(nfreq * sizeof(double));

	// layout
	const double x1 = 0.12 * width;
	const double x2 = 0.90 * width;
	const double y1 = 0.10 * height;
	const double y2 = 0.90 * height;
	const double h = font;

	// data : min, max
	double dmin = +1000;
	double dmax = -1000;
	for (int ifeed = 0; ifeed < nfeed; ifeed++) {
		for (int ipoint = 0; ipoint < npoint; ipoint++) {
			for (int ifreq = 0; ifreq < nfreq; ifreq++) {
				const double d = 20 * log10(MAX(couple[ifeed][ipoint][ifreq], eps));
				dmin = MIN(dmin, d);
				dmax = MAX(dmax, d);
			}
		}
	}

	// scale : min, max, div
	double ymin, ymax;
	int ydiv;
	if (scale.user == 0) {
		const double dbmax = +100;
		const double dbmin = -300;
		const double dbunit = 5;
		ymax = dbunit *  (ceil)(dmax / dbunit);
		ymin = dbunit * (floor)(dmin / dbunit);
		ymax = MIN(ymax, dbmax);
		ymin = MAX(ymin, dbmin);
		ydiv = (int)((ymax - ymin) / dbunit + 1e-6);
	}
	else {
		ymin = scale.min;
		ymax = scale.max;
		ydiv = scale.div;
	}

	// plot
	ev2d_newPage();

	// grid
	ev2dlib_grid(x1, y1, x2, y2, fdiv, ydiv);

	ev2d_setColor(0, 0, 0);

	for (int ifeed = 0; ifeed < nfeed; ifeed++) {
		for (int ipoint = 0; ipoint < npoint; ipoint++) {
			for (int ifreq = 0; ifreq < nfreq; ifreq++) {
				f[ifreq] = 20 * log10(MAX(couple[ifeed][ipoint][ifreq], eps));
			}

			if      (ipoint == 0) {
				ev2d_setColor(255,   0,   0);
			}
			else if (ipoint == 1) {
				ev2d_setColor(  0,   0, 255);
			}
			else if (ipoint == 2) {
				ev2d_setColor(  0, 255,   0);
			}
			else {
				ev2d_setColor(rand() & 0xff, rand() & 0xff, rand() & 0xff); 
			}

			ev2dlib_func1(nfreq - 1, f, ymin, ymax, x1, y1, x2, y2);

			sprintf(str1, "C%d%d", ipoint + 1, ifeed + 1);
			double y = y1 + (y2 - y1) * (f[nfreq - 1] - ymin) / (ymax - ymin);
			ev2d_drawString(x2, y - 0.4 * h, h, str1);
		}
	}

	ev2d_setColor(0, 0, 0);

	// min, max
	sprintf(str1, "max=%.3f[dB]", dmax);
	sprintf(str2, "min=%.3f[dB]", dmin);
	ev2d_drawString(x2 - 10.0 * h, y2 + 1.8 * h, h, str1);
	ev2d_drawString(x2 - 10.0 * h, y2 + 0.4 * h, h, str2);

	// X-axis
	sprintf(str1, "%.3e", freq[0]);
	sprintf(str2, "%.3e", freq[nfreq - 1]);
	ev2dlib_Xaxis(x1, x2, y1, h, str1, str2, "frequency[Hz]");

	// Y-axis
	sprintf(str1, "%.0f", ymin);
	sprintf(str2, "%.0f", ymax);
	ev2dlib_Yaxis(y1, y2, x1, h, str1, str2, "[dB]");

	// title
	ev2d_drawString(x1, y2 + 0.4 * h, h, title);
	ev2d_drawString(x1, y2 + 1.8 * h, h, "Coupling");

	// free
	free(f);
}


// f-char
void plot2dFchar(
	int nfreq, const double fchar[], scale_t scale, int fdiv,
	const double freq[], const char title[], const char figure[], const char sunit[],
	double width, double height, double font)
{
	char str1[BUFSIZ], str2[BUFSIZ], str3[BUFSIZ];
	const double eps = 1e-30;

	// layout
	const double x1 = 0.12 * width;
	const double x2 = 0.90 * width;
	const double y1 = 0.10 * height;
	const double y2 = 0.90 * height;
	const double h = font;

	// alloc
	double *f = (double *)malloc(nfreq * sizeof(double));

	// to dB
	for (int ifreq = 0; ifreq < nfreq; ifreq++) {
		f[ifreq] = 10 * log10(MAX(fchar[ifreq], eps));
	}

	// data : min, max
	double dmin = f[0];
	double dmax = dmin;
	for (int ifreq = 1; ifreq < nfreq; ifreq++) {
		const double d = f[ifreq];
		dmin = MIN(dmin, d);
		dmax = MAX(dmax, d);
	}

	// scale : min, max, div
	double ymin, ymax;
	int ydiv;
	if (scale.user == 0) {
		const double dbmax = +100;
		const double dbmin = -300;
		const double dbunit = 5;
		ymax = dbunit *  (ceil)(dmax / dbunit);
		ymin = dbunit * (floor)(dmin / dbunit);
		ymax = MIN(ymax, dbmax);
		ymin = MAX(ymin, dbmin);
		ydiv = (int)((ymax - ymin) / dbunit + 1e-6);
	}
	else {
		ymin = scale.min;
		ymax = scale.max;
		ydiv = scale.div;
	}

	// plot
	ev2d_newPage();

	// grid
	ev2dlib_grid(x1, y1, x2, y2, fdiv, ydiv);

	ev2d_setColor(0, 0, 0);

	// function
	ev2dlib_func1(nfreq - 1, f, ymin, ymax, x1, y1, x2, y2);

	// max, min
	sprintf(str1, "max=%.3f[%s]", dmax, sunit);
	sprintf(str2, "min=%.3f[%s]", dmin, sunit);
	ev2d_drawString(x2 - 10.0 * h, y2 + 1.8 * h, h, str1);
	ev2d_drawString(x2 - 10.0 * h, y2 + 0.4 * h, h, str2);

	// X-axis
	sprintf(str1, "%.3e", freq[0]);
	sprintf(str2, "%.3e", freq[nfreq - 1]);
	ev2dlib_Xaxis(x1, x2, y1, h, str1, str2, "frequency[Hz]");

	// Y-axis
	sprintf(str1, "%.0f", ymin);
	sprintf(str2, "%.0f", ymax);
	sprintf(str3, "[%s]", sunit);
	ev2dlib_Yaxis(y1, y2, x1, h, str1, str2, str3);

	// title
	ev2d_drawString(x1, y2 + 1.8 * h, h, figure);
	ev2d_drawString(x1, y2 + 0.4 * h, h, title);

	// free
	free(f);
}


static void smithChart(int nfrequency, double *rin, double *xin, double z0, double x0, double y0, double r0)
{
	char str[BUFSIZ];
	const int cdiv = 90, nr = 6, nx = 6;
	const double h = 0.05 * r0;
	const double rdat[] = {0.2, 0.5, 1.0, 2.0, 5.0, 10.0};
	const double xdat[] = {0.2, 0.5, 1.0, 2.0, 5.0, 10.0};
	const double pdat[] = {-0.923, -0.596, 0.0, 0.596, 0.923, 0.981};
	const double pi = 4 * atan(1);

	// plot constant r lines

	ev2d_setColor(0, 0, 0);
	ev2d_drawEllipse(x0 - r0, y0 - r0, x0 + r0, y0 + r0);
	ev2d_drawString(x0 - r0 - 1.0 * h, y0 - 0.5 * h, h, "0");

	for (int ir = 0; ir < nr; ir++) {
		const double xc = x0 + (r0 * rdat[ir]) / (rdat[ir] + 1);
		const double rc = r0 / (rdat[ir] + 1);
		ev2d_setColor(150, 150, 150);
		ev2d_drawEllipse(xc - rc, y0 - rc, xc + rc, y0 + rc);
		ev2d_setColor(0, 0, 0);
		sprintf(str, "%d", (int)(z0 * rdat[ir]));
		ev2d_drawString(xc - rc - 0.7 * h, y0 - 0.4 * h, h, str);
	}

	// plot constant x lines

	ev2d_setColor(150, 150, 150);
	ev2d_drawLine(x0 - r0, y0, x0 + r0, y0);
	for (int ix = 0; ix < nx; ix++) {
		const int n = cdiv * (nx + 1) / (ix + 1);
		double *xa = (double *)malloc(n * sizeof(double));
		double *ya = (double *)malloc(n * sizeof(double));
		const double r = 1 / xdat[ix];
		int nl = 0;
		for (int i = 0; i < n; i++) {
			const double arg = (1.5 * pi) - (pi * i / n);
			const double x = 1 + r * cos(arg);
			const double y = + r + r * sin(arg);
			if ((x * x) + (y * y) <= 1) {
				xa[i] = x0 + (r0 * x);
				ya[i] = y0 + (r0 * y);
				nl = i;
			}
		}
		ev2d_drawPolyline(nl + 1, xa, ya);
		nl = 0;
		for (int i = 0; i < n; i++) {
			const double arg = (0.5 * pi) + (pi * i / n);
			const double x = 1 + r * cos(arg);
			const double y = - r + r * sin(arg);
			if ((x * x) + (y * y) <= 1) {
				xa[i] = x0 + (r0 * x);
				ya[i] = y0 + (r0 * y);
				nl = i;
			}
		}
		ev2d_drawPolyline(nl + 1, xa, ya);
		free(xa);
		free(ya);
	}

	ev2d_setColor(0, 0, 0);
	for (int ix = 0; ix < nx; ix++) {
		const double tmp = sqrt(1 - (pdat[ix] * pdat[ix]));
		const double x = x0 + (r0 * pdat[ix]);
		sprintf(str, "%d", (int)(+z0 * xdat[ix]));
		ev2d_drawString(x - 0.8 * h, y0 + r0 * tmp - 0.4 * h, h, str);
		sprintf(str, "%d", (int)(-z0 * xdat[ix]));
		ev2d_drawString(x - 1.3 * h, y0 - r0 * tmp - 0.4 * h, h, str);
	}

	// plot

	ev2d_setColor(0, 0, 0);
	double *xb = (double *)malloc(nfrequency * sizeof(double));
	double *yb = (double *)malloc(nfrequency * sizeof(double));
	double c1[2], c2[2], c3[2];

	for (int i = 0; i < nfrequency; i++) {
		c1[0] = rin[i] - z0;
		c1[1] = xin[i];
		c2[0] = rin[i] + z0;
		c2[1] = xin[i];
		c3[0] = (c1[0] * c2[0] + c1[1] * c2[1]) / (c2[0] * c2[0] + c2[1] * c2[1]);
		c3[1] = (c1[1] * c2[0] - c1[0] * c2[1]) / (c2[0] * c2[0] + c2[1] * c2[1]);
		xb[i] = x0 + (r0 * c3[0]);
		yb[i] = y0 + (r0 * c3[1]);
		if      (i == 0) {
			const double dm = 0.4 * h;
			ev2d_fillRectangle(xb[i] - dm, yb[i] - dm, xb[i] + dm, yb[i] + dm);
		}
		else if (i == nfrequency / 2) {
			const double dm = 0.5 * h;
			ev2d_fillEllipse(xb[i] - dm, yb[i] - dm, xb[i] + dm, yb[i] + dm);
		}
		else if (i == nfrequency - 1) {
			const double dm = 0.4 * h;
			ev2d_drawRectangle(xb[i] - dm, yb[i] - dm, xb[i] + dm, yb[i] + dm);
		}
		else {
			const double dm = 0.15 * h;
			ev2d_fillRectangle(xb[i] - dm, yb[i] - dm, xb[i] + dm, yb[i] + dm);
		}
	}

	ev2d_setColor(255,   0,   0);
	ev2d_drawPolyline(nfrequency, xb, yb);
	ev2d_setColor(0, 0, 0);

	free(xb);
	free(yb);
}
