/*
plot2dIter.c

plot E/H average of iteration (2D)
*/

#include "orcwa.h"
#include "ev.h"

void plot2dIter(void)
{
	if (Niter < 1) return;

	char str1[BUFSIZ], str2[BUFSIZ];
	const int ylog = 6;

	// layout
	const double x1 = 0.12 * Width2d;
	const double x2 = 0.90 * Width2d;
	const double y1 = 0.10 * Height2d;
	const double y2 = 0.90 * Height2d;
	const double h = Fontsize2d;

	// min, max
	double fmax = 0;
	for (int n = 0; n < Niter; n++) {
		fmax = MAX(fmax, MAX(Eiter[n], Hiter[n]));
	}
	double ymax = ceil(log10(fmax));
	double ymin = ymax - ylog;
	//printf("%f %f\n", ymax, ymin);

	// plot
	ev2d_newPage();

	ev2dlib_grid(x1, y1, x2, y2, 10, ylog);

	// <E>
	ev2d_setColor(255,   0,   0);
	for (int n = 0; n < Niter - 1; n++) {
		double px1 = x1 + (x2 - x1) * (n + 0) / (Niter - 1);
		double px2 = x1 + (x2 - x1) * (n + 1) / (Niter - 1);
		double py1 = MAX(log10(MAX(Eiter[n + 0], EPS2)), ymin);
		double py2 = MAX(log10(MAX(Eiter[n + 1], EPS2)), ymin);
		py1 = y1 + (y2 - y1) * (py1 - ymin) / (ymax - ymin);
		py2 = y1 + (y2 - y1) * (py2 - ymin) / (ymax - ymin);
		ev2d_drawLine(px1, py1, px2, py2);
	}
	ev2d_drawString(x2 - 6.0 * h, y2 + 0.5 * h, h, "<E>");

	// <H>
	ev2d_setColor(  0,   0, 255);
	for (int n = 0; n < Niter - 1; n++) {
		double px1 = x1 + (x2 - x1) * (n + 0) / (Niter - 1);
		double px2 = x1 + (x2 - x1) * (n + 1) / (Niter - 1);
		double py1 = MAX(log10(MAX(Hiter[n + 0], EPS2)), ymin);
		double py2 = MAX(log10(MAX(Hiter[n + 1], EPS2)), ymin);
		py1 = y1 + (y2 - y1) * (py1 - ymin) / (ymax - ymin);
		py2 = y1 + (y2 - y1) * (py2 - ymin) / (ymax - ymin);
		ev2d_drawLine(px1, py1, px2, py2);
	}
	ev2d_drawString(x2 - 3.0 * h, y2 + 0.5 * h, h, "<H>");

	ev2d_setColor(  0,   0,   0);

	// title
	ev2d_drawString(x1, y2 + 0.4 * h, h, Title);
	ev2d_drawString(x1, y2 + 1.8 * h, h, "iteration");

	// X-axis
	sprintf(str2, "%d", (Niter - 1) * Solver.nout);
	ev2dlib_Xaxis(x1, x2, y1, h, "0", str2, "time step");
	sprintf(str1, "%.3e[sec]", (Niter - 1) * Solver.nout * Dt);
	ev2d_drawString(x2 - 7.0 * h, y1 - 2.5 * h, h, str1);

	// Y-axis
	sprintf(str1, "%.0e", pow(10, ymin));
	sprintf(str2, "%.0e", pow(10, ymax));
	ev2dlib_Yaxis(y1, y2, x1, h, str1, str2, "");
}
