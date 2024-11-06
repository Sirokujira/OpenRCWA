/*
outputNear1d.c

near1d field (2D plot and log)
*/

#include "orcwa.h"
#include "complex.h"
#include "orcwa_prototype.h"


void outputNear1d(void)
{
	double *amp[3], *deg[3], pos[2][3];
	double *lng = NULL;
	const char fmt[] = "%4d%12.3e%12.3e%12.3e%11.3e%11.3e%9.3f%11.3e%9.3f%11.3e%9.3f\n";

	// open log file
	FILE *fp;
	if ((fp = fopen(FN_near1d, "w")) == NULL) {
		fprintf(stderr, "*** %s open error.\n", FN_near1d);
		return;
	}

	int adr = 0;
	for (int n = 0; n < NNear1d; n++) {

		// division and position
		int div = 0;
		if      (Near1d[n].dir == 'X') {
			lng = Xn;
			div = Nx;
			pos[0][0] = lng[0];
			pos[1][0] = lng[div];
			pos[0][1] = pos[1][1] = Yn[Near1d[n].id1];
			pos[0][2] = pos[1][2] = Zn[Near1d[n].id2];
		}
		else if (Near1d[n].dir == 'Y') {
			lng = Yn;
			div = Ny;
			pos[0][1] = lng[0];
			pos[1][1] = lng[div];
			pos[0][2] = pos[1][2] = Zn[Near1d[n].id1];
			pos[0][0] = pos[1][0] = Xn[Near1d[n].id2];
		}
		else if (Near1d[n].dir == 'Z') {
			lng = Zn;
			div = Nz;
			pos[0][2] = lng[0];
			pos[1][2] = lng[div];
			pos[0][0] = pos[1][0] = Xn[Near1d[n].id1];
			pos[0][1] = pos[1][1] = Yn[Near1d[n].id2];
		}

		// alloc
		for (int m = 0; m < 3; m++) {
			amp[m] = (double *)malloc((div + 1) * sizeof(double));
			deg[m] = (double *)malloc((div + 1) * sizeof(double));
		}

		for (int ifreq = 0; ifreq < NFreq2; ifreq++) {
			// log
			fprintf(fp, "#%d : frequency[Hz] = %.5e\n", (n + 1), Freq2[ifreq]);
			fprintf(fp, "%s", " No.     X[m]        Y[m]        Z[m]      ");
			fprintf(fp, "%s\n", (Near1d[n].cmp[0] == 'E')
				? "E[V/m]     Ex[V/m]   Ex[deg]   Ey[V/m]   Ey[deg]   Ez[V/m]   Ez[deg]"
				: "H[A/m]     Hx[A/m]   Hx[deg]   Hy[A/m]   Hy[deg]   Hz[A/m]   Hz[deg]");

			// field
			for (int i = 0; i <= div; i++) {
				// position
				double x = 0, y = 0, z = 0;
				if      (Near1d[n].dir == 'X') {
					x = Xn[i];
					y = Yn[Near1d[n].id1];
					z = Zn[Near1d[n].id2];
				}
				else if (Near1d[n].dir == 'Y') {
					x = Xn[Near1d[n].id2];
					y = Yn[i];
					z = Zn[Near1d[n].id1];
				}
				else if (Near1d[n].dir == 'Z') {
					x = Xn[Near1d[n].id1];
					y = Yn[Near1d[n].id2];
					z = Zn[i];
				}

				// E or H (complex)
				d_complex_t c[3], e[3], h[3];
				if      (Near1d[n].cmp[0] == 'E') {
					c[0] = Near1dEx[adr];
					c[1] = Near1dEy[adr];
					c[2] = Near1dEz[adr];
					if (IPlanewave && !Near1dNoinc) {
						planewave(Freq2[ifreq], x, y, z, e, h);
						for (int m = 0; m < 3; m++) {
							c[m] = d_add(c[m], e[m]);
						}
					}
				}
				else if (Near1d[n].cmp[0] == 'H') {
					c[0] = Near1dHx[adr];
					c[1] = Near1dHy[adr];
					c[2] = Near1dHz[adr];
					if (IPlanewave && !Near1dNoinc) {
						planewave(Freq2[ifreq], x, y, z, e, h);
						for (int m = 0; m < 3; m++) {
							c[m] = d_add(c[m], h[m]);
						}
					}
					for (int m = 0; m < 3; m++) {
						c[m] = d_rmul(1 / ETA0, c[m]);
					}
				}
				adr++;

				// amplitude and phase
				for (int m = 0; m < 3; m++) {
					amp[m][i] = d_abs(c[m]);
					deg[m][i] = d_deg(c[m]);
				}

				// log
				fprintf(fp, fmt, i, x, y, z,
					sqrt(d_norm(c[0]) + d_norm(c[1]) + d_norm(c[2])),
					d_abs(c[0]), d_deg(c[0]),
					d_abs(c[1]), d_deg(c[1]),
					d_abs(c[2]), d_deg(c[2]));
			}

			// plot
			plot2dNear1d0(
				Near1d[n].cmp, div, amp, deg, lng,
				Near1dScale.db, Near1dScale.user, Near1dScale.min, Near1dScale.max, Near1dScale.div,
				Title, Freq2[ifreq], pos,
				Width2d, Height2d, Fontsize2d);
		}

		// free
		for (int m = 0; m < 3; m++) {
			free(amp[m]);
			free(deg[m]);
		}
	}

	// close
	fclose(fp);
}
