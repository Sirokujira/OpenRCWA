/*
plot2dPoint.c

plot waveform and spectrum on points (2D)
*/

#include "orcwa.h"
#include "complex.h"
#include "orcwa_prototype.h"
#include "ev.h"

void plot2dPoint(void)
{
	if (!NPoint || (Ntime < 2)) return;

	char str1[BUFSIZ], str2[BUFSIZ];
	const int ylog = 8;

	// layout
	const double x1 = 0.12 * Width2d;
	const double x2 = 0.90 * Width2d;
	const double y1 = 0.10 * Height2d;
	const double y2 = 0.90 * Height2d;
	const double h = Fontsize2d;

	// open
	FILE *fp;
	if ((fp = fopen(FN_point, "w")) == NULL) {
		fprintf(stderr, "*** %s open error.\n", FN_point);
		return;
	}

	// === waveform ===

	for (int ipoint = 0; ipoint < NPoint; ipoint++) {
		// min, max
		double fmax, ymax, ymin;
		fmax = 0;
		for (int itime = 0; itime < Ntime; itime++) {
			int id = ipoint * (Solver.maxiter + 1) + itime;
			//printf("%d %d %e\n", ipoint, itime, VPoint[id]);
			fmax = MAX(fmax, fabs(VPoint[id]));
		}
		//printf("%f\n", fmax);
		ymax = ceil(log10(fmax));
		ymin = ymax - ylog;
		//printf("%f %f\n", ymax, ymin);

		// plot
		ev2d_newPage();

		// grid
		ev2dlib_grid(x1, y1, x2, y2, 10, ylog);

		ev2d_setColor(0, 0, 0);

		// E (dB)
		for (int itime = 0; itime < Ntime - 1; itime++) {
			int id = ipoint * (Solver.maxiter + 1) + itime;
			double px1 = x1 + (x2 - x1) * (itime + 0) / (Ntime - 1);
			double px2 = x1 + (x2 - x1) * (itime + 1) / (Ntime - 1);
			double py1 = log10(MAX(fabs(VPoint[id + 0]), EPS2));
			double py2 = log10(MAX(fabs(VPoint[id + 1]), EPS2));
			py1 = y1 + (y2 - y1) * (py1 - ymin) / (ymax - ymin);
			py2 = y1 + (y2 - y1) * (py2 - ymin) / (ymax - ymin);
			py1 = MAX(py1, y1);
			py2 = MAX(py2, y1);
			ev2d_drawLine(px1, py1, px2, py2);
		}

		// title
		ev2d_drawString(x1, y2 + 0.4 * h, h, Title);
		ev2d_drawString(x1, y2 + 1.8 * h, h, "E waveform on points");

		// Y-axis
		sprintf(str1, "%.0e", pow(10, ymin));
		sprintf(str2, "%.0e", pow(10, ymax));
		ev2dlib_Yaxis(y1, y2, x1, h, str1, str2, "[V/m]");

		// X-axis
		sprintf(str2, "%d", Ntime - 1);
		ev2dlib_Xaxis(x1, x2, y1, h, "0", str2, "time step");

		sprintf(str1, "%.3e[sec]", (Ntime - 1) * Dt);
		ev2d_drawString(x2 - 8.0 * h, y1 - 2.5 * h, h, str1);

		// point #
		sprintf(str1, "point #%d", ipoint + 1);
		ev2d_drawString(x2 - 6.0 * h, y2 + 0.5 * h, h, str1);

		// log
		fprintf(fp, "point #%d (waveform)\n", ipoint + 1);
		fprintf(fp, "%s\n", "    No.    time[sec]      E[V/m]");
		for (int itime = 0; itime < Ntime; itime++) {
			int id = ipoint * (Solver.maxiter + 1) + itime;
			fprintf(fp, "%7d %13.5e %13.5e\n", itime, itime * Dt, fabs(VPoint[id]));
		}
	}

	// === spectrum ===

	d_complex_t *esp = (d_complex_t *)malloc(NFreq1 * sizeof(d_complex_t));

	for (int ipoint = 0; ipoint < NPoint; ipoint++) {
		// DFT
		for (int ifreq = 0; ifreq < NFreq1; ifreq++) {
			esp[ifreq] = calcdft(Ntime, &VPoint[ipoint * (Solver.maxiter + 1)], Freq1[ifreq], Dt, 0);
		}

		// normalize amplitude
		double emax = 0;
		for (int ifreq = 0; ifreq < NFreq1; ifreq++) {
			emax = MAX(emax, d_abs(esp[ifreq]));
		}
		for (int ifreq = 0; ifreq < NFreq1; ifreq++) {
			esp[ifreq] = d_rmul(1 / emax, esp[ifreq]);
		}

		// plot
		ev2d_newPage();

		// grid
		ev2dlib_grid(x1, y1, x2, y2, 10, ylog);

		ev2d_setColor(0, 0, 0);

		// E (linear)
		const double ymax = 1;
		const double ymin = 0;
		for (int ifreq = 0; ifreq < NFreq1 - 1; ifreq++) {
			double px1 = x1 + (x2 - x1) * (ifreq + 0) / (NFreq1 - 1);
			double px2 = x1 + (x2 - x1) * (ifreq + 1) / (NFreq1 - 1);
			double py1 = d_abs(esp[ifreq + 0]);
			double py2 = d_abs(esp[ifreq + 1]);
			py1 = y1 + (y2 - y1) * (py1 - ymin) / (ymax - ymin);
			py2 = y1 + (y2 - y1) * (py2 - ymin) / (ymax - ymin);
			ev2d_drawLine(px1, py1, px2, py2);
		}

		// title
		ev2d_drawString(x1, y2 + 0.4 * h, h, Title);
		ev2d_drawString(x1, y2 + 1.8 * h, h, "E spectrum on points");

		// Y-axis
		sprintf(str1, "%.0f", ymin);
		sprintf(str2, "%.0f", ymax);
		ev2dlib_Yaxis(y1, y2, x1, h, str1, str2, "[V/m]");

		// X-axis
		sprintf(str1, "%.3e", Freq1[0]);
		sprintf(str2, "%.3e", Freq1[NFreq1 - 1]);
		ev2dlib_Xaxis(x1, x2, y1, h, str1, str2, "frequency[Hz]");

		// point #
		sprintf(str1, "point #%d", ipoint + 1);
		ev2d_drawString(x2 - 6.0 * h, y2 + 0.5 * h, h, str1);

		// log
		fprintf(fp, "point #%d (spectrum)\n", ipoint + 1);
		fprintf(fp, "%s\n", " No. frequency[Hz]  amplitude   degree");
		for (int ifreq = 0; ifreq < NFreq1; ifreq++) {
			fprintf(fp, "%4d %13.5e %9.5f %9.3f\n", ifreq, Freq1[ifreq], d_abs(esp[ifreq]), d_deg(esp[ifreq]));
		}
	}

	// free
	free(esp);

	// close
	fclose(fp);
}
