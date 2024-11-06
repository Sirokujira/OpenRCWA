/*
plot2dFeed.c

plot waveform and spectrum on feeds (2D)
*/

#include "orcwa.h"
#include "complex.h"
#include "orcwa_prototype.h"
#include "ev.h"

void plot2dFeed(void)
{
	if (!NFeed || (Ntime < 2)) return;

	char str1[BUFSIZ], str2[BUFSIZ];
	const char fig1[2][BUFSIZ] = {"V[V]", "I[A]"};
	const char fig2[2][BUFSIZ] = {"V", "I"};
	const unsigned char rgb[2][3] = {{255, 0, 0}, {0, 0, 255}};
	const int ylog = 8;

	// layout
	const double x1 = 0.12 * Width2d;
	const double x2 = 0.90 * Width2d;
	const double y1 = 0.10 * Height2d;
	const double y2 = 0.90 * Height2d;
	const double h = Fontsize2d;

	// open
	FILE *fp;
	if ((fp = fopen(FN_feed, "w")) == NULL) {
		fprintf(stderr, "*** %s open error.\n", FN_feed);
		return;
	}

	// === waveform ===

	for (int ifeed = 0; ifeed < NFeed; ifeed++) {
		// min, max
		double fmax[2], ymax[2], ymin[2];
		fmax[0] = fmax[1] = 0;
		for (int itime = 0; itime < Ntime; itime++) {
			int id = ifeed * (Solver.maxiter + 1) + itime;
			fmax[0] = MAX(fmax[0], fabs(VFeed[id]));
			fmax[1] = MAX(fmax[1], fabs(IFeed[id]));
		}
		for (int m = 0; m < 2; m++) {
			ymax[m] = ceil(log10(fmax[m]));
			ymin[m] = ymax[m] - ylog;
		}

		// plot
		ev2d_newPage();

		// grid
		ev2dlib_grid(x1, y1, x2, y2, 10, ylog);

		// V, I [dB]
		for (int m = 0; m < 2; m++) {
			ev2d_setColor(rgb[m][0], rgb[m][1], rgb[m][2]);
			double *ptr = (m == 0) ? VFeed : IFeed;
			for (int itime = 0; itime < Ntime - 1; itime++) {
				int id = ifeed * (Solver.maxiter + 1) + itime;
				double px1 = x1 + (x2 - x1) * (itime + 0) / (Ntime - 1);
				double px2 = x1 + (x2 - x1) * (itime + 1) / (Ntime - 1);
				double py1 = MAX(log10(MAX(fabs(ptr[id + 0]), EPS2)), ymin[m]);
				double py2 = MAX(log10(MAX(fabs(ptr[id + 1]), EPS2)), ymin[m]);
				py1 = y1 + (y2 - y1) * (py1 - ymin[m]) / (ymax[m] - ymin[m]);
				py2 = y1 + (y2 - y1) * (py2 - ymin[m]) / (ymax[m] - ymin[m]);
				ev2d_drawLine(px1, py1, px2, py2);
			}

			// Y-axis
			double x = (m == 0) ? x1 - 5.0 * h : x2 + 0.3 * h;
			sprintf(str1, "%.0e", pow(10, ymax[m]));
			sprintf(str2, "%.0e", pow(10, ymin[m]));
			ev2d_drawString(x, y2 - 0.2 * h, h, str1);
			ev2d_drawString(x, y1 - 0.2 * h, h, str2);
			ev2d_drawString(x, y2 - 1.5 * h, h, fig1[m]);
		}

		ev2d_setColor(0, 0, 0);

		// title
		ev2d_drawString(x1, y2 + 0.4 * h, h, Title);
		ev2d_drawString(x1, y2 + 1.8 * h, h, "V and I waveform on feeds");

		// X-axis
		sprintf(str2, "%d", Ntime -1);
		ev2dlib_Xaxis(x1, x2, y1, h, "0", str2, "time step");

		sprintf(str1, "%.3e[sec]", (Ntime - 1) * Dt);
		ev2d_drawString(x2 - 8.0 * h, y1 - 2.5 * h, h, str1);

		// feed #
		sprintf(str1, "feed #%d", ifeed + 1);
		ev2d_drawString(x2 - 6.0 * h, y2 + 0.5 * h, h, str1);

		// log
		fprintf(fp, "feed #%d (waveform)\n", ifeed + 1);
		fprintf(fp, "%s\n", "    No.    time[sec]      V[V]          I[A]");
		for (int itime = 0; itime < Ntime; itime++) {
			int id = ifeed * (Solver.maxiter + 1) + itime;
			fprintf(fp, "%7d %13.5e %13.5e %13.5e\n", itime, itime * Dt, fabs(VFeed[id]), fabs(IFeed[id]));
		}
	}

	// === spectrum ===

	double *vsp = (double *)malloc(NFreq1 * sizeof(double));
	double *isp = (double *)malloc(NFreq1 * sizeof(double));

	for (int ifeed = 0; ifeed < NFeed; ifeed++) {
		// DFT
		for (int ifreq = 0; ifreq < NFreq1; ifreq++) {
			const d_complex_t vsum = calcdft(Ntime, &VFeed[ifeed * (Solver.maxiter + 1)], Freq1[ifreq], Dt, 0);
			const d_complex_t isum = calcdft(Ntime, &IFeed[ifeed * (Solver.maxiter + 1)], Freq1[ifreq], Dt, 0);
			vsp[ifreq] = d_abs(vsum);
			isp[ifreq] = d_abs(isum);
			//printf("%d %d %e %e\n", ifeed, ifreq, vsp[ifreq], isp[ifreq]);
		}

		// normalize
		double vmax = 0, imax = 0;
		for (int ifreq = 0; ifreq < NFreq1; ifreq++) {
			vmax = MAX(vmax, vsp[ifreq]);
			imax = MAX(imax, isp[ifreq]);
		}
		for (int ifreq = 0; ifreq < NFreq1; ifreq++) {
			vsp[ifreq] /= vmax;
			isp[ifreq] /= imax;
		}
		const double ymax = 1;
		const double ymin = 0;

		// plot
		ev2d_newPage();

		// grid
		ev2dlib_grid(x1, y1, x2, y2, 10, ylog);

		// V, I (linear)
		for (int m = 0; m < 2; m++) {
			ev2d_setColor(rgb[m][0], rgb[m][1], rgb[m][2]);
			double *ptr = (m == 0) ? vsp : isp;
			for (int ifreq = 0; ifreq < NFreq1 - 1; ifreq++) {
				double px1 = x1 + (x2 - x1) * (ifreq + 0) / (NFreq1 - 1);
				double px2 = x1 + (x2 - x1) * (ifreq + 1) / (NFreq1 - 1);
				double py1 = ptr[ifreq + 0];
				double py2 = ptr[ifreq + 1];
				py1 = y1 + (y2 - y1) * (py1 - ymin) / (ymax - ymin);
				py2 = y1 + (y2 - y1) * (py2 - ymin) / (ymax - ymin);
				ev2d_drawLine(px1, py1, px2, py2);
			}

			// V or I
			double x = x2 - 4.0 * h + 1.5 * m * h;
			ev2d_drawString(x, y2 + 1.7 * h, h, fig2[m]);
		}

		ev2d_setColor(0, 0, 0);

		// title
		ev2d_drawString(x1, y2 + 0.4 * h, h, Title);
		ev2d_drawString(x1, y2 + 1.8 * h, h, "V and I spectrum on feeds");

		// Y-axis
		sprintf(str1, "%.0f", ymin);
		sprintf(str2, "%.0f", ymax);
		ev2dlib_Yaxis(y1, y2, x1, h, str1, str2, "");

		// X-axis
		sprintf(str1, "%.3e", Freq1[0]);
		sprintf(str2, "%.3e", Freq1[NFreq1 - 1]);
		ev2dlib_Xaxis(x1, x2, y1, h, str1, str2, "frequency[Hz]");

		// feed #
		sprintf(str1, "feed #%d", ifeed + 1);
		ev2d_drawString(x2 - 6.0 * h, y2 + 0.5 * h, h, str1);

		// log
		fprintf(fp, "feed #%d (spectrum)\n", ifeed + 1);
		fprintf(fp, "%s\n", " No. frequency[Hz]       V          I");
		for (int ifreq = 0; ifreq < NFreq1; ifreq++) {
			fprintf(fp, "%4d %13.5e %10.5f %10.5f\n", ifreq, Freq1[ifreq], vsp[ifreq], isp[ifreq]);
		}
	}

	// free
	free(vsp);
	free(isp);

	// close
	fclose(fp);
}
