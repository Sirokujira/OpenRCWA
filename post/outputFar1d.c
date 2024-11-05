/*
outputFar1d.c

output far1d field
*/

#include "orcwa.h"
#include "complex.h"
#include "orcwa_prototype.h"


static void plot2dFar1d(d_complex_t ***etheta, d_complex_t ***ephi)
{
	double (*ef)[7];
	char str[BUFSIZ];

	for (int n = 0; n < NFar1d; n++) {
		// alloc
		ef = (double (*)[7])malloc((Far1d[n].div + 1) * 7 * sizeof(double));

		for (int ifreq = 0; ifreq < NFreq2; ifreq++) {
			// E
			for (int i = 0; i <= Far1d[n].div; i++) {
				farComponent(etheta[n][ifreq][i], ephi[n][ifreq][i], ef[i]);
			}

			// unit
			if (Far1dScale.db) {
				strcpy(str, (NFeed ? "[dB]" : "[dBsm]"));
			}
			else {
				strcpy(str, (NFeed ? "" : "[sm]"));
			}

			// plot
			plot2dFar1d0(
				Far1d[n].div, ef,
				Far1dComp, Far1d[n].dir, Far1d[n].angle, Far1dStyle,
				Far1dScale.db, Far1dScale.user, Far1dScale.min, Far1dScale.max, Far1dScale.div,
				Title, str, Freq2[ifreq],
				Width2d, Height2d, Fontsize2d);
		}

		// free
		free(ef);
	}
}


static void logFar1d(d_complex_t ***etheta, d_complex_t ***ephi)
{
	double e[7];

	FILE *fp;
	if ((fp = fopen(FN_far1d, "w")) == NULL) {
		fprintf(stderr, "*** %s open error.\n", FN_far1d);
		return;
	}

	for (int n = 0; n < NFar1d; n++) {
		for (int ifreq = 0; ifreq < NFreq2; ifreq++) {
			fprintf(fp, "#%d : %c-plane", (n + 1), Far1d[n].dir);
			if ((Far1d[n].dir == 'V') || (Far1d[n].dir == 'H')) {
				fprintf(fp, " (%s = %.2f[deg])", (Far1d[n].dir == 'V' ? "phi" : "theta"), Far1d[n].angle);
			}
			fprintf(fp, ", frequency[Hz] = %.5e\n", Freq2[ifreq]);
			fprintf(fp, "  No.   deg    E-abs[dB]  E-theta[dB] E-theta[deg]    E-phi[dB]   E-phi[deg]  E-major[dB]  E-minor[dB]   E-RHCP[dB]   E-LHCP[dB] AxialRatio[dB]\n");
			for (int i = 0; i <= Far1d[n].div; i++) {
				farComponent(etheta[n][ifreq][i], ephi[n][ifreq][i], e);
				// to dB
				for (int k = 0; k < 7; k++) {
					e[k] = 20 * log10(MAX(e[k], EPS2));
				}
				// output
				double angle = (360.0 * i) / Far1d[n].div;
				fprintf(fp, "%4d%7.1f%13.4f%13.4f%13.4f%13.4f%13.4f%13.4f%13.4f%13.4f%13.4f%13.4f\n",
					i, angle, e[0], e[1], d_deg(etheta[n][ifreq][i]), e[2], d_deg(ephi[n][ifreq][i]), e[3], e[4], e[5], e[6], e[3] - e[4]);
			}
		}
	}

	fclose(fp);
}


void outputFar1d(void)
{
	// alloc
	size_t size = NFar1d * sizeof(d_complex_t **);
	d_complex_t ***etheta = (d_complex_t ***)malloc(size);
	d_complex_t ***ephi   = (d_complex_t ***)malloc(size);
	for (int n = 0; n < NFar1d; n++) {
		size = NFreq2 * sizeof(d_complex_t *);
		etheta[n] = (d_complex_t **)malloc(size);
		ephi[n]   = (d_complex_t **)malloc(size);
		for (int ifreq = 0; ifreq < NFreq2; ifreq++) {
			size = (Far1d[n].div + 1) * sizeof(d_complex_t);
			etheta[n][ifreq] = (d_complex_t *)malloc(size);
			ephi[n][ifreq]   = (d_complex_t *)malloc(size);
			memset(etheta[n][ifreq], 0, size);
			memset(  ephi[n][ifreq], 0, size);
		}
	}
	double *ffctr = (double *)malloc(NFreq2 * sizeof(double));

	// setup
	setup_farfield();

	// factor
	for (int ifreq = 0; ifreq < NFreq2; ifreq++) {
		ffctr[ifreq] = farfactor(ifreq);
	}

	// calculation
	calcFar1d(ffctr, etheta, ephi);

	// plot
	plot2dFar1d(etheta, ephi);

	// write log
	logFar1d(etheta, ephi);

	// free
	for (int n = 0; n < NFar1d; n++) {
		for (int ifreq = 0; ifreq < NFreq2; ifreq++) {
			free(etheta[n][ifreq]);
			free(ephi[n][ifreq]);
		}
		free(etheta[n]);
		free(ephi[n]);
	}
	free(etheta);
	free(ephi);

	free(ffctr);
}
