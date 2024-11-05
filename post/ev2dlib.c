/*
ev2dlib.c

ev2d utilities (2D)
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "ev.h"


// plot a grid
void ev2dlib_grid(double x1, double y1, double x2, double y2, int xdiv, int ydiv)
{
	ev2d_setColor(  0,   0,   0);
	ev2d_drawRectangle(x1, y1, x2, y2);

	ev2d_setColor(192, 192, 192);
	for (int i = 1; i < xdiv; i++) {
		double x = x1 + (x2 - x1) * i / xdiv;
		ev2d_drawLine(x, y1, x, y2);
	}
	for (int j = 1; j < ydiv; j++) {
		double y = y1 + (y2 - y1) * j / ydiv;
		ev2d_drawLine(x1, y, x2, y);
	}

	ev2d_setColor(  0,   0,   0);
}


// plot a function : y[x[0]]...y[x[div]]
void ev2dlib_func2(
	int div, const double x[], const double y[], double ymin, double ymax,
	double x1, double y1, double x2, double y2)
{
	if (div < 0) return;
	if (ymin == ymax) return;

	double yf = (y2 - y1) / (ymax - ymin);
	if (div > 0) {
		double xf = (x2 - x1) / (x[div] - x[0]);
		for (int i = 0; i < div; i++) {
			double xa = x1 + xf * (x[i]     - x[0]);
			double xb = x1 + xf * (x[i + 1] - x[0]);
			double ya = y1 + yf * (y[i]     - ymin);
			double yb = y1 + yf * (y[i + 1] - ymin);
			if (ya < y1) ya = y1;
			if (ya > y2) ya = y2;
			if (yb < y1) yb = y1;
			if (yb > y2) yb = y2;
			ev2d_drawLine(xa, ya, xb, yb);
		}
	}
	else {
		double y0 = y1 + yf * (y[0] - ymin);
		if ((y0 >= y1) && (y0 <= y2)) {
			ev2d_drawLine(x1, y0, x2, y0);
		}
	}
}


// plot a function : f[0]...f[div]
void ev2dlib_func1(
	int div, const double f[], double fmin, double fmax,
	double x1, double y1, double x2, double y2)
{
	if (div < 0) return;
	if (fmin == fmax) return;

	double *x = (double *)malloc((div + 1) * sizeof(double));

	if (div > 0) {
		for (int i = 0; i <= div; i++) {
			x[i] = (double)i / div;
		}
	}

	ev2dlib_func2(div, x, f, fmin, fmax, x1, y1, x2, y2);

	free(x);
}


// plot color sample
void ev2dlib_sample(double x1, double y1, double x2, double y2, int color)
{
	const int div = 20;

	for (int n = 0; n < div; n++) {
		ev2d_setColorV((n + 0.5) / div, color);
		ev2d_fillRectangle(x1, y1 + (y2 - y1) * (n + 0.0) / div,
		                   x2, y1 + (y2 - y1) * (n + 1.0) / div);
	}

	ev2d_setColor(0, 0, 0);
}


// (static) plot contour in a triangle
static void triangleContour(double d[3][3], double zmin, double zmax, int color)
{
	const double eps = 1e-10;
	const int ncolor = 20;

	// check

	if (zmin >= zmax) return;
	if (fabs((d[1][0] - d[0][0]) * (d[2][1] - d[0][1])
	       - (d[2][0] - d[0][0]) * (d[1][1] - d[0][1])) < eps * eps ) return;

	// ordering by z : 0 < 1 < 2

	double z0 = d[0][2];
	double z1 = d[1][2];
	double z2 = d[2][2];

	int iv0, iv1, iv2;
	if      ((z0 <= z1) && (z1 <= z2)) {iv0 = 0; iv1 = 1; iv2 = 2;}
	else if ((z0 <= z2) && (z2 <= z1)) {iv0 = 0; iv1 = 2; iv2 = 1;}
	else if ((z1 <= z2) && (z2 <= z0)) {iv0 = 1; iv1 = 2; iv2 = 0;}
	else if ((z1 <= z0) && (z0 <= z2)) {iv0 = 1; iv1 = 0; iv2 = 2;}
	else if ((z2 <= z0) && (z0 <= z1)) {iv0 = 2; iv1 = 0; iv2 = 1;}
	else if ((z2 <= z1) && (z1 <= z0)) {iv0 = 2; iv1 = 1; iv2 = 0;}
	else                               {iv0 = 0; iv1 = 1; iv2 = 2;}
	double *v0 = d[iv0];
	double *v1 = d[iv1];
	double *v2 = d[iv2];

	// height index of three corners

	const double dz = (zmax - zmin) / ncolor;

	int iz0 = (int)floor((v0[2] - zmin) / dz);
	int iz1 = (int)floor((v1[2] - zmin) / dz);
	int iz2 = (int)floor((v2[2] - zmin) / dz);

	// contour

	if ((iz0 == iz1) && (iz1 == iz2)) {
		ev2d_setColorV((double)iz0 / ncolor, color);
		ev2d_fillTriangle(v0[0], v0[1], v1[0], v1[1], v2[0], v2[1]);
	}
	else {
		// 0 -> 1 & 0 -> 2
		if (iz0 < iz1) {
			for (int iz = iz0; iz <= iz1; iz++) {
				double za = zmin + (iz + 0) * dz;
				double zb = zmin + (iz + 1) * dz;
				double a01 = (fabs(v1[2] - v0[2]) > eps) ? (za - v0[2]) / (v1[2] - v0[2]) : 0;
				double a02 = (fabs(v2[2] - v0[2]) > eps) ? (za - v0[2]) / (v2[2] - v0[2]) : 0;
				double b01 = (fabs(v1[2] - v0[2]) > eps) ? (zb - v0[2]) / (v1[2] - v0[2]) : 1;
				double b02 = (fabs(v2[2] - v0[2]) > eps) ? (zb - v0[2]) / (v2[2] - v0[2]) : 1;
				a01 = (a01 < 0) ? 0 : (a01 > 1) ? 1 : a01;
				a02 = (a02 < 0) ? 0 : (a02 > 1) ? 1 : a02;
				b01 = (b01 < 0) ? 0 : (b01 > 1) ? 1 : b01;
				b02 = (b02 < 0) ? 0 : (b02 > 1) ? 1 : b02;
				double x1 = v0[0] + a01 * (v1[0] - v0[0]);
				double y1 = v0[1] + a01 * (v1[1] - v0[1]);
				double x2 = v0[0] + b01 * (v1[0] - v0[0]);
				double y2 = v0[1] + b01 * (v1[1] - v0[1]);
				double x3 = v0[0] + b02 * (v2[0] - v0[0]);
				double y3 = v0[1] + b02 * (v2[1] - v0[1]);
				double x4 = v0[0] + a02 * (v2[0] - v0[0]);
				double y4 = v0[1] + a02 * (v2[1] - v0[1]);
				ev2d_setColorV((double)iz / ncolor, color);
				ev2d_fillQuadrangle(x1, y1, x2, y2, x3, y3, x4, y4);
			}
		}

		// 1 -> 2 & 0 -> 2
		if (iz1 < iz2) {
			for (int iz = iz1; iz <= iz2; iz++) {
				double za = zmin + (iz + 0) * dz;
				double zb = zmin + (iz + 1) * dz;
				double a02 = (fabs(v2[2] - v0[2]) > eps) ? (za - v0[2]) / (v2[2] - v0[2]) : 0;
				double a12 = (fabs(v2[2] - v1[2]) > eps) ? (za - v1[2]) / (v2[2] - v1[2]) : 0;
				double b02 = (fabs(v2[2] - v0[2]) > eps) ? (zb - v0[2]) / (v2[2] - v0[2]) : 1;
				double b12 = (fabs(v2[2] - v1[2]) > eps) ? (zb - v1[2]) / (v2[2] - v1[2]) : 1;
				a02 = (a02 < 0) ? 0 : (a02 > 1) ? 1 : a02;
				a12 = (a12 < 0) ? 0 : (a12 > 1) ? 1 : a12;
				b02 = (b02 < 0) ? 0 : (b02 > 1) ? 1 : b02;
				b12 = (b12 < 0) ? 0 : (b12 > 1) ? 1 : b12;
				double x1 = v0[0] + a02 * (v2[0] - v0[0]);
				double y1 = v0[1] + a02 * (v2[1] - v0[1]);
				double x2 = v0[0] + b02 * (v2[0] - v0[0]);
				double y2 = v0[1] + b02 * (v2[1] - v0[1]);
				double x3 = v1[0] + b12 * (v2[0] - v1[0]);
				double y3 = v1[1] + b12 * (v2[1] - v1[1]);
				double x4 = v1[0] + a12 * (v2[0] - v1[0]);
				double y4 = v1[1] + a12 * (v2[1] - v1[1]);
				ev2d_setColorV((double)iz / ncolor, color);
				ev2d_fillQuadrangle(x1, y1, x2, y2, x3, y3, x4, y4);
			}
		}
	}
}


// (static) plot contour in a rectangle
static void rectangleContour(double d[4][3], double zmin, double zmax, int color)
{
	double t[3][3];

	// triangle 0-1-2
	for (int n = 0; n < 3; n++) {
		t[0][n] = d[0][n];
		t[1][n] = d[1][n];
		t[2][n] = d[2][n];
	}
	triangleContour(t, zmin, zmax, color);

	// triangle 0-2-3
	for (int n = 0; n < 3; n++) {
		t[0][n] = d[0][n];
		t[1][n] = d[2][n];
		t[2][n] = d[3][n];
	}
	triangleContour(t, zmin, zmax, color);
}


/*
plot contour in a rectangle region
non-uniform mesh
mode = 0/1 : color/gray fine
       2/3 : color/gray coarse
*/
void ev2dlib_contour(
	int nxdiv, int nydiv, const double x[], const double y[], double **z,
	double zmin, double zmax, int mode)
{
	if ((nxdiv <= 0) || (nydiv <= 0)) return;
	if ((mode < 0) || (mode > 3)) return;
	if (zmin >= zmax) {
		ev2d_setColor(0, 0, 255);
		ev2d_fillRectangle(x[0], y[0], x[nxdiv], y[nydiv]);
		ev2d_setColor(0, 0, 0);
		ev2d_drawRectangle(x[0], y[0], x[nxdiv], y[nydiv]);
		return;
	}

	double quad[4][3];

	for (int i = 0; i < nxdiv; i++) {
	for (int j = 0; j < nydiv; j++) {
		quad[0][0] = x[i    ];
		quad[1][0] = x[i + 1];
		quad[2][0] = x[i + 1];
		quad[3][0] = x[i    ];
		quad[0][1] = y[j    ];
		quad[1][1] = y[j    ];
		quad[2][1] = y[j + 1];
		quad[3][1] = y[j + 1];
		quad[0][2] = z[i    ][j    ];
		quad[1][2] = z[i + 1][j    ];
		quad[2][2] = z[i + 1][j + 1];
		quad[3][2] = z[i    ][j + 1];

		if ((mode == 0) || (mode == 1)) {
			// fine (color or gray)
			rectangleContour(quad, zmin, zmax, (mode == 0));
		}
		else if ((mode == 2) || (mode == 3)) {
			// coarse (color or gray)
			double zav = (quad[0][2] + quad[1][2] + quad[2][2] + quad[3][2]) / 4;
			double f = (zav - zmin) / (zmax - zmin);
			ev2d_setColorV(f, (mode == 2));
			ev2d_fillRectangle(quad[0][0], quad[0][1], quad[2][0], quad[2][1]);
		}
	}
	}

	ev2d_setColor(0, 0, 0);
	ev2d_drawRectangle(x[0], y[0], x[nxdiv], y[nydiv]);
}


// draw an arrow
void ev2dlib_arrow(double x1, double y1, double x2, double y2)
{
	double dx = x2 - x1;
	double dy = y2 - y1;

	ev2d_drawLine(x1, y1, x2, y2);
	ev2d_drawLine(x2, y2, x1 + 0.7 * dx - 0.15 * dy, y1 + 0.7 * dy + 0.15 * dx);
	ev2d_drawLine(x2, y2, x1 + 0.7 * dx + 0.15 * dy, y1 + 0.7 * dy - 0.15 * dx);
}


// plot X-axis
void ev2dlib_Xaxis(double x1, double x2, double y1, double fs, const char *smin, const char *smax, const char *caption)
{
	ev2d_drawString(x1 - 0.3 * strlen(smin) * fs, y1 - 1.5 * fs, fs, smin);
	ev2d_drawString(x2 - 0.3 * strlen(smax) * fs, y1 - 1.5 * fs, fs, smax);
	ev2d_drawString((x1 + x2) / 2 - 0.3 * strlen(caption) * fs, y1 - 1.5 * fs, fs, caption);
}


// plot Y-axis
void ev2dlib_Yaxis(double y1, double y2, double x1, double fs, const char *smin, const char *smax, const char *sdim)
{
	ev2d_drawString(x1 - 0.75 * strlen(smin) * fs - 0.5 * fs, y1 - 0.4 * fs, fs, smin);
	ev2d_drawString(x1 - 0.75 * strlen(smax) * fs - 0.5 * fs, y2 - 0.4 * fs, fs, smax);
	ev2d_drawString(x1 - 0.75 * strlen(sdim) * fs - 0.5 * fs, y2 - 2.0 * fs, fs, sdim);
}


// plot a function in the circle
void ev2dlib_CircleFunc(double x1, double y1, double x2, double y2, const double f[], int div, double fmin, double fmax, double p0, double sgn)
{
	const double pi = 4 * atan(1);

	double x0 = (x1 + x2) / 2;
	double y0 = (y1 + y2) / 2;
	double r0 = fabs(y2 - y1) / 2;
	double rf = r0 / (fmax - fmin);

	for (int i = 0; i < div; i++) {
		double pa = p0 + sgn * (2 * pi * (i + 0)) / div;
		double pb = p0 + sgn * (2 * pi * (i + 1)) / div;
		double ra = rf * (f[i + 0] - fmin);
		double rb = rf * (f[i + 1] - fmin);
		if (ra <  0) ra =  0;
		if (ra > r0) ra = r0;
		if (rb <  0) rb =  0;
		if (rb > r0) rb = r0;
		double xa = x0 + ra * cos(pa);
		double ya = y0 + ra * sin(pa);
		double xb = x0 + rb * cos(pb);
		double yb = y0 + rb * sin(pb);
		ev2d_drawLine(xa, ya, xb, yb);
	}
}


// plot circle mesh
void ev2dlib_CircleMesh(double x1, double y1, double x2, double y2, int rdiv, int adiv)
{
	const unsigned char gray = 200;
	const unsigned char gray2 = 130;
	const double pi = 4 * atan(1);

	double x0 = (x1 + x2) / 2;
	double y0 = (y1 + y2) / 2;
	double r0 = fabs(y2 - y1) / 2;

	// circles
	for (int ir = 1; ir <= rdiv; ir++) {
		unsigned char ic = (ir < rdiv) ? gray : 0;
		ev2d_setColor(ic, ic, ic);
		double r = r0 * ir / rdiv;
		ev2d_drawEllipse(x0 - r, y0 - r, x0 + r, y0 + r);
	}

	// radial lines
	for (int ia = 0; ia < adiv; ia++) {
		unsigned char ic = (ia % (adiv / 2)) ? gray : gray2;
		ev2d_setColor(ic, ic, ic);
		double xa = x0 + r0 * cos(pi * ia / adiv);
		double ya = y0 + r0 * sin(pi * ia / adiv);
		double xb = 2 * x0 - xa;
		double yb = 2 * y0 - ya;
		ev2d_drawLine(xa, ya, xb, yb);
	}

	ev2d_setColor(0, 0, 0);
}


// array argument
void ev2dlib_drawLineA(double g[2][2])
{
	ev2d_drawLine(g[0][0], g[0][1],
	              g[1][0], g[1][1]);
}
void ev2dlib_drawTriangleA(double g[3][2])
{
	ev2d_drawTriangle(g[0][0], g[0][1],
	                  g[1][0], g[1][1],
	                  g[2][0], g[2][1]);
}
void ev2dlib_fillTriangleA(double g[3][2])
{
	ev2d_fillTriangle(g[0][0], g[0][1],
	                  g[1][0], g[1][1],
	                  g[2][0], g[2][1]);
}
void ev2dlib_drawQuadrangleA(double g[4][2])
{
	ev2d_drawQuadrangle(g[0][0], g[0][1],
	                    g[1][0], g[1][1],
	                    g[2][0], g[2][1],
	                    g[3][0], g[3][1]);
}
void ev2dlib_fillQuadrangleA(double g[4][2])
{
	ev2d_fillQuadrangle(g[0][0], g[0][1],
	                    g[1][0], g[1][1],
	                    g[2][0], g[2][1],
	                    g[3][0], g[3][1]);
}
