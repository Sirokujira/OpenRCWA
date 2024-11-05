/*
plot2dNear2d0.c (OpenMOM/OpenRCWA/OpenTHFD)

plot near field on a plane (2D)
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "ev.h"

void plot2dNear2d0(
	int div1, int div2, double **d, const double pos1[], const double pos2[],
	int phase, int scaledb, int scaleuser, double scalemin, double scalemax, int contour,
	const char title[], double frequency,
	char dir, double pos0, const char component[],
	int ngline, double (*gline)[2][3],
	double width, double height, double font)
{
	if ((div1 <= 0) || (div2 <= 0)) return;
	if (phase) scaledb = 0;

	const double dbauto = 40;
	const double fw = 0.7;  // width margin
	const double fh = 0.8;  // height margin
	char str[4 * BUFSIZ], strmin[BUFSIZ], strmax[BUFSIZ], strunit[BUFSIZ];

	// layout
	const double w0 = width / 2.0;
	const double h0 = height / 2.0;
	const double h = font;

	const double x1 = pos1[0];
	const double x2 = pos1[div1];
	const double y1 = pos2[0];
	const double y2 = pos2[div2];

	const double x0 = (x1 + x2) / 2;
	const double y0 = (y1 + y2) / 2;

	const double xfctr = fw * width  / fabs(x2 - x1);
	const double yfctr = fh * height / fabs(y2 - y1);
	const double fctr = (xfctr < yfctr) ? xfctr : yfctr;

	const double px0 = w0 - (fctr * x0);
	const double py0 = h0 - (fctr * y0);

	const double rx1 = px0 + (fctr * x1);
	const double rx2 = px0 + (fctr * x2);
	const double ry1 = py0 + (fctr * y1);
	const double ry2 = py0 + (fctr * y2);

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

	// node
	double *xn = (double *)malloc((div1 + 1) * sizeof(double));
	double *yn = (double *)malloc((div2 + 1) * sizeof(double));
	for (int i = 0; i <= div1; i++) {
		xn[i] = px0 + (fctr * pos1[i]);
	}
	for (int j = 0; j <= div2; j++) {
		yn[j] = py0 + (fctr * pos2[j]);
	}

	// new page
	ev2d_newPage();

	// contour
	int mode = (contour == 0) ? 0 : (contour == 1) ? 2 : (contour == 2) ? 1 : 3;
	ev2dlib_contour(div1, div2, xn, yn, d, fmin, fmax, mode);

	ev2d_setColor(0, 0, 0);

	// border
	ev2d_drawRectangle(rx1, ry1, rx2, ry2);

	// comment
	sprintf(strunit, "[%s%s]", (scaledb ? "dB" : ""), (phase ? "deg" : (component[0] == 'E') ? "V/m" : (component[0] == 'H') ? "A/m" : "A"));
	if (!phase) {
		sprintf(strmax, (scaledb ? "max=%.3f%s" : "max=%.3e%s"), dmax, strunit);
	}
	else{
		strcpy(strmax, "");
	}
	sprintf(str, "%s%s %c=%.3e[m] f=%.3e[Hz] %s", component, strunit, dir, pos0, frequency, strmax);
	ev2d_drawString(2.0 * h, ry2 + 1.8 * h, h, str);

	// axis
	char cx[BUFSIZ] = "";
	char cy[BUFSIZ] = "";
	if      (dir == 'X') {
		// Y-Z
		strcpy(cx, "Y[m]");
		strcpy(cy, "Z[m]");
	}
	else if (dir == 'Y') {
		// X-Z
		strcpy(cx, "X[m]");
		strcpy(cy, "Z[m]");
	}
	else if (dir == 'Z') {
		// X-Y
		strcpy(cx, "X[m]");
		strcpy(cy, "Y[m]");
	}
	ev2d_drawString(w0 - 2.0 * h, ry1 - 1.5 * h, h, cx);
	ev2d_drawString(rx1 - 3.5 * h, h0 - 0.5 * h, h, cy);

	sprintf(str, "%.3e", x1);
	ev2d_drawString(rx1 - 2.0 * h, ry1 - 1.5 * h, 0.8 * h, str);
	sprintf(str, "%.3e", x2);
	ev2d_drawString(rx2 - 2.0 * h, ry1 - 1.5 * h, 0.8 * h, str);

	sprintf(str, "%.3e", y1);
	ev2d_drawString(rx1 - 7.0 * h, ry1 - 0.3 * h, 0.8 * h, str);
	sprintf(str, "%.3e", y2);
	ev2d_drawString(rx1 - 7.0 * h, ry2 - 0.3 * h, 0.8 * h, str);

	// title
	ev2d_drawString(rx1, ry2 + 0.4 * h, h, title);

	// color sample
	ev2dlib_sample(rx2 + 0.5 * h, ry2 - 5.0 * h, rx2 + 1.5 * h, ry2, ((contour == 0) || (contour == 1)));
	ev2d_setColor(0, 0, 0);
	if (!phase) {
		sprintf(strmax, (scaledb ? "%.3f" : "%.3e"), fmax);
		sprintf(strmin, (scaledb ? "%.3f" : "%.3e"), fmin);
	}
	else {
		sprintf(strmax, "%+d", (int)(fmax + 0.5));
		sprintf(strmin, "%+d", (int)(fmin - 0.5));
	}
	ev2d_drawString(rx2 + 1.5 * h, ry2 - 0.3 * h, 0.8 * h, strmax);
	ev2d_drawString(rx2 + 1.5 * h, ry2 - 5.3 * h, 0.8 * h, strmin);
	ev2d_drawString(rx2 + 1.7 * h, ry2 - 1.3 * h, 0.8 * h, strunit);

	// border of objects
	if (ngline > 0) {
		double px1 = 0, py1 = 0, px2 = 0, py2 = 0;
		for (int i = 0; i < ngline; i++) {
			if      (dir == 'X') {
				// Y-Z
				px1 = px0 + (fctr * gline[i][0][1]);
				px2 = px0 + (fctr * gline[i][1][1]);
				py1 = py0 + (fctr * gline[i][0][2]);
				py2 = py0 + (fctr * gline[i][1][2]);
			}
			else if (dir == 'Y') {
				// X-Z
				px1 = px0 + (fctr * gline[i][0][0]);
				px2 = px0 + (fctr * gline[i][1][0]);
				py1 = py0 + (fctr * gline[i][0][2]);
				py2 = py0 + (fctr * gline[i][1][2]);
			}
			else if (dir == 'Z') {
				// X-Y
				px1 = px0 + (fctr * gline[i][0][0]);
				px2 = px0 + (fctr * gline[i][1][0]);
				py1 = py0 + (fctr * gline[i][0][1]);
				py2 = py0 + (fctr * gline[i][1][1]);
			}
			if ((px1 >= rx1) && (px1 <= rx2) &&
			    (px2 >= rx1) && (px2 <= rx2) &&
			    (py1 >= ry1) && (py1 <= ry2) &&
			    (py2 >= ry1) && (py2 <= ry2)) {
				ev2d_drawLine(px1, py1, px2, py2);
			}
		}
	}

	// free
	free(xn);
	free(yn);
}
