/*
outputNear2d.c

plot and log of near2d (2D+3D)
*/

#include "orcwa.h"
#include "complex.h"
#include "orcwa_prototype.h"
#include "ev.h"

static void field(int, int, int, int, int64_t, d_complex_t *, double []);
static void fcomp(const char [], const d_complex_t [], double *, double *);
static void fcomp2(double, const char [], const d_complex_t [], double *);


void outputNear2d(void)
{
	const int ngline = Near2dObj ? (int)NGline : 0;
	const double eps = 1e-10;
	const char fmt[] = "%5d%5d%12.3e%12.3e%12.3e%11.3e%11.3e%9.3f%11.3e%9.3f%11.3e%9.3f\n";

	// open log file
	FILE *fp;
	if ((fp = fopen(FN_near2d, "w")) == NULL) {
		fprintf(stderr, "*** %s open error.\n", FN_near2d);
		return;
	}

	int64_t adr = 0;
	for (int n = 0; n < NNear2d; n++) {
		int    div1 = 0, div2 = 0;
		int    start1 = 0, stop1 = 0;
		int    start2 = 0, stop2 = 0;
		double pos0 = 0, *pos1 = NULL, *pos2 = NULL;
		if      (Near2d[n].dir == 'X') {
			// Y-Z
			start1 = !Near2dIzoom ?  0 : nearest(Near2dHzoom[0], 0, Ny, Yn);
			stop1  = !Near2dIzoom ? Ny : nearest(Near2dHzoom[1], 0, Ny, Yn);
			start2 = !Near2dIzoom ?  0 : nearest(Near2dVzoom[0], 0, Nz, Zn);
			stop2  = !Near2dIzoom ? Nz : nearest(Near2dVzoom[1], 0, Nz, Zn);
			pos0 = Xn[Near2d[n].id0];
			pos1 = &Yn[start1];
			pos2 = &Zn[start2];
			div1 = Ny;
			div2 = Nz;
		}
		else if (Near2d[n].dir == 'Y') {
			// X-Z
			start1 = !Near2dIzoom ?  0 : nearest(Near2dHzoom[0], 0, Nx, Xn);
			stop1  = !Near2dIzoom ? Nx : nearest(Near2dHzoom[1], 0, Nx, Xn);
			start2 = !Near2dIzoom ?  0 : nearest(Near2dVzoom[0], 0, Nz, Zn);
			stop2  = !Near2dIzoom ? Nz : nearest(Near2dVzoom[1], 0, Nz, Zn);
			pos0 = Yn[Near2d[n].id0];
			pos1 = &Xn[start1];
			pos2 = &Zn[start2];
			div1 = Nx;
			div2 = Nz;
		}
		else if (Near2d[n].dir == 'Z') {
			// X-Y
			start1 = !Near2dIzoom ?  0 : nearest(Near2dHzoom[0], 0, Nx, Xn);
			stop1  = !Near2dIzoom ? Nx : nearest(Near2dHzoom[1], 0, Nx, Xn);
			start2 = !Near2dIzoom ?  0 : nearest(Near2dVzoom[0], 0, Ny, Yn);
			stop2  = !Near2dIzoom ? Ny : nearest(Near2dVzoom[1], 0, Ny, Yn);
			pos0 = Zn[Near2d[n].id0];
			pos1 = &Xn[start1];
			pos2 = &Yn[start2];
			div1 = Nx;
			div2 = Ny;
		}

		// alloc (zoom)
		const int adiv1 = stop1 - start1;
		const int adiv2 = stop2 - start2;
		double **amp = (double **)malloc((adiv1 + 1) * sizeof(double *));
		double **deg = (double **)malloc((adiv1 + 1) * sizeof(double *));
		for (int n1 = 0; n1 <= adiv1; n1++) {
			amp[n1] = (double *)malloc((adiv2 + 1) * sizeof(double));
			deg[n1] = (double *)malloc((adiv2 + 1) * sizeof(double));
		}
		//printf("%d %d\n", adiv1, adiv2);

		for (int ifreq = 0; ifreq < NFreq2; ifreq++) {
			// log
			fprintf(fp, "#%d : frequency[Hz] = %.5e\n", (n + 1), Freq2[ifreq]);
			fprintf(fp, "%s", "  No.  No.     X[m]        Y[m]        Z[m]      ");
			fprintf(fp, "%s\n", (Near2d[n].cmp[0] == 'E')
				? "E[V/m]     Ex[V/m]   Ex[deg]   Ey[V/m]   Ey[deg]   Ez[V/m]   Ez[deg]"
				: "H[A/m]     Hx[A/m]   Hx[deg]   Hy[A/m]   Hy[deg]   Hz[A/m]   Hz[deg]");

			// field and log
			for (int n1 = 0; n1 <= div1; n1++) {
			for (int n2 = 0; n2 <= div2; n2++) {
				if ((n1 >= start1) && (n1 <= stop1) &&
				    (n2 >= start2) && (n2 <= stop2)) {  // zoom
					// field
					d_complex_t f[3];
					double pos[3];
					field(n, ifreq, n1, n2, adr, f, pos);
					fcomp(Near2d[n].cmp, f, &amp[n1 - start1][n2 - start2], &deg[n1 - start1][n2 - start2]);

					// log
					fprintf(fp, fmt,
						n1, n2, pos[0], pos[1], pos[2],
						sqrt(d_norm(f[0]) + d_norm(f[1]) + d_norm(f[2])),
						d_abs(f[0]), d_deg(f[0]),
						d_abs(f[1]), d_deg(f[1]),
						d_abs(f[2]), d_deg(f[2]));
				}
				// counter++
				adr++;
			}
			}

			// to dB
			if (Near2dScale.db) {
				for (int n1 = 0; n1 <= adiv1; n1++) {
				for (int n2 = 0; n2 <= adiv2; n2++) {
					amp[n1][n2] = 20 * log10(MAX(amp[n1][n2], EPS2));
				}
				}
			}

			// plot (2D)
			if (Near2dDim[0]) {
				// harmonic
				if (!Near2dFrame) {
					plot2dNear2d0(
						adiv1, adiv2, amp, pos1, pos2,
						0, Near2dScale.db, Near2dScale.user, Near2dScale.min, Near2dScale.max, Near2dContour,
						Title, Freq2[ifreq],
						Near2d[n].dir, pos0, Near2d[n].cmp,
						ngline, Gline,
						Width2d, Height2d, Fontsize2d);

					// phase
					if (strcmp(Near2d[n].cmp, "E") && strcmp(Near2d[n].cmp, "H")) {
						plot2dNear2d0(
							adiv1, adiv2, deg, pos1, pos2,
							1, 0, 1, -180, +180, Near2dContour,
							Title, Freq2[ifreq],
							Near2d[n].dir, pos0, Near2d[n].cmp,
							ngline, Gline,
							Width2d, Height2d, Fontsize2d);
					}
				}
				// animation
				else if ((n == 0) && (ifreq == 0)) {
					// alloc (zoom)
					double **mag = (double **)malloc((adiv1 + 1) * sizeof(double *));
					for (int n1 = 0; n1 <= adiv1; n1++) {
						mag[n1] = (double *)malloc((adiv2 + 1) * sizeof(double));
					}
					// auto scale -> user scale
					scale_t scale = Near2dScale;
					if (scale.user == 0) {
						// data : max
						double dmax = 0;
						for (int frame = 0; frame < Near2dFrame; frame++) {
							int64_t adr0 = 0;  // 1st plane and 1st frequency
							for (int n1 = 0; n1 <= div1; n1++) {
							for (int n2 = 0; n2 <= div2; n2++) {
							if ((n1 >= start1) && (n1 <= stop1) &&
							    (n2 >= start2) && (n2 <= stop2)) {  // zoom
									d_complex_t f[3];
									double pos[3], wf;
									field(n, ifreq, n1, n2, adr0, f, pos);
									const double wt = 2 * PI * frame / Near2dFrame;
									fcomp2(wt, Near2d[n].cmp, f, &wf);
									dmax = MAX(dmax, wf);
								}
								adr0++;
							}
							}
						}
						scale.user = 1;
						if (scale.db) {
							// dB
							scale.max = 20 * log10(MAX(dmax, eps));
							scale.min = scale.max - 40;
						}
						else {
							// linear
							scale.max = dmax;
							scale.min = 0;
						}
					}
					// frame
					for (int frame = 0; frame < Near2dFrame; frame++) {
						// data
						int64_t adr0 = 0;  // 1st plane and 1st frequency
						for (int n1 = 0; n1 <= div1; n1++) {
						for (int n2 = 0; n2 <= div2; n2++) {
							if ((n1 >= start1) && (n1 <= stop1) &&
							    (n2 >= start2) && (n2 <= stop2)) {  // zoom
								d_complex_t f[3];
								double pos[3];
								field(n, ifreq, n1, n2, adr0, f, pos);
								const double wt = 2 * PI * frame / Near2dFrame;
								fcomp2(wt, Near2d[n].cmp, f, &mag[n1 - start1][n2 - start2]);
							}
							adr0++;
						}
						}
						// to dB
						if (scale.db) {
							for (int n1 = 0; n1 <= adiv1; n1++) {
							for (int n2 = 0; n2 <= adiv2; n2++) {
								mag[n1][n2] = 20 * log10(MAX(mag[n1][n2], eps));
							}
							}
						}
						// plot 2D
						plot2dNear2d0(
							adiv1, adiv2, mag, pos1, pos2,
							0, scale.db, scale.user, scale.min, scale.max, Near2dContour,
							Title, Freq2[ifreq],
							Near2d[n].dir, Near2d[n].pos0, Near2d[n].cmp,
							ngline, Gline,
							Width2d, Height2d, Fontsize2d);
					}
					// free
					for (int n1 = 0; n1 <= adiv1; n1++) {
						free(mag[n1]);
					}
					free(mag);
				}
			}

			// plot (3D)
			if (Near2dDim[1]) {
				plot3dNear2d0(
					adiv1, adiv2, amp, pos1, pos2,
					0, Near2dScale.db, Near2dScale.user, Near2dScale.min, Near2dScale.max, Near2dContour,
					Title, Freq2[ifreq],
					Near2d[n].dir, pos0, Near2d[n].cmp,
					ngline, Gline);
				ev3d_setColor(160, 160, 160);
				ev3d_index(2);
				ev3d_drawBox(Xn[0], Yn[0], Zn[0], Xn[Nx], Yn[Ny], Zn[Nz]);

				// phase
				if (strcmp(Near2d[n].cmp, "E") && strcmp(Near2d[n].cmp, "H")) {
					plot3dNear2d0(
						adiv1, adiv2, deg, pos1, pos2,
						1, 0, 1, -180, +180, Near2dContour,
						Title, Freq2[ifreq],
						Near2d[n].dir, pos0, Near2d[n].cmp,
						ngline, Gline);
					ev3d_setColor(160, 160, 160);
					ev3d_index(2);
					ev3d_drawBox(Xn[0], Yn[0], Zn[0], Xn[Nx], Yn[Ny], Zn[Nz]);
				}
			}
		}

		// free
		for (int n1 = 0; n1 <= adiv1; n1++) {
			free(amp[n1]);
			free(deg[n1]);
		}
		free(amp);
		free(deg);
	}

	// close
	fclose(fp);
}


static void field(int n, int ifreq, int n1, int n2, int64_t adr, d_complex_t *f, double pos[])
{
	// position
	if      (Near2d[n].dir == 'X') {
		pos[0] = Xn[Near2d[n].id0];
		pos[1] = Yn[n1];
		pos[2] = Zn[n2];
	}
	else if (Near2d[n].dir == 'Y') {
		pos[0] = Xn[n1];
		pos[1] = Yn[Near2d[n].id0];
		pos[2] = Zn[n2];
	}
	else if (Near2d[n].dir == 'Z') {
		pos[0] = Xn[n1];
		pos[1] = Yn[n2];
		pos[2] = Zn[Near2d[n].id0];
	}

	// E or H (complex)
	d_complex_t e[3], h[3];
	if      (Near2d[n].cmp[0] == 'E') {
		f[0] = Near2dEx[adr];
		f[1] = Near2dEy[adr];
		f[2] = Near2dEz[adr];
		if (IPlanewave && !Near2dNoinc) {
			planewave(Freq2[ifreq], pos[0], pos[1], pos[2], e, h);
			for (int m = 0; m < 3; m++) {
				f[m] = d_add(f[m], e[m]);
			}
		}
	}
	else if (Near2d[n].cmp[0] == 'H') {
		f[0] = Near2dHx[adr];
		f[1] = Near2dHy[adr];
		f[2] = Near2dHz[adr];
		if (IPlanewave && !Near2dNoinc) {
			planewave(Freq2[ifreq], pos[0], pos[1], pos[2], e, h);
			for (int m = 0; m < 3; m++) {
				f[m] = d_add(f[m], h[m]);
			}
		}
		for (int m = 0; m < 3; m++) {
			f[m] = d_rmul(1 / ETA0, f[m]);
		}
	}
}


// amplitude and phase
static void fcomp(const char comp[], const d_complex_t f[], double *amp, double *deg)
{
	const int c = toupper(comp[1]);

	if      (c == 'X') {
		*amp = d_abs(f[0]);
		*deg = d_deg(f[0]);
	}
	else if (c == 'Y') {
		*amp = d_abs(f[1]);
		*deg = d_deg(f[1]);
	}
	else if (c == 'Z') {
		*amp = d_abs(f[2]);
		*deg = d_deg(f[2]);
	}
	else {
		*amp = sqrt(d_norm(f[0])
		          + d_norm(f[1])
		          + d_norm(f[2]));
		*deg = 0;
	}
}


// waveform
static void fcomp2(double wt, const char comp[], const d_complex_t f[], double *mag)
{
	const int c = toupper(comp[1]);

	if      (c == 'X') {
		*mag = d_abs(f[0]) * cos(d_rad(f[0]) + wt);
	}
	else if (c == 'Y') {
		*mag = d_abs(f[1]) * cos(d_rad(f[1]) + wt);
	}
	else if (c == 'Z') {
		*mag = d_abs(f[2]) * cos(d_rad(f[2]) + wt);
	}
	else {
		*mag = sqrt(d_norm(f[0]) * cos(d_rad(f[0]) + wt) * cos(d_rad(f[0]) + wt)
		          + d_norm(f[1]) * cos(d_rad(f[1]) + wt) * cos(d_rad(f[1]) + wt)
		          + d_norm(f[2]) * cos(d_rad(f[2]) + wt) * cos(d_rad(f[2]) + wt));
	}
}
