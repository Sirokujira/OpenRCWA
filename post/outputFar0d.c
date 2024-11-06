/*
outputFar0d.c

output far0d field
*/

#include "orcwa.h"
#include "complex.h"
#include "orcwa_prototype.h"

static void calcFar0d(const double ffctr[], d_complex_t etheta[], d_complex_t ephi[], double theta, double phi)
{
	for (int ifreq = 0; ifreq < NFreq2; ifreq++) {
		farfield(ifreq, theta, phi, ffctr[ifreq], &etheta[ifreq], &ephi[ifreq]);
	}
}


static void plot2dFar0d(const d_complex_t etheta[], const d_complex_t ephi[], double theta, double phi)
{
	char str1[BUFSIZ], str2[BUFSIZ];
	double e[7];

	// alloc
	double **ef = (double **)malloc(7 * sizeof(double *));
	for (int comp = 0; comp < 7; comp++) {
		ef[comp] = (double *)malloc(NFreq2 * sizeof(double));
	}

	// data
	for (int ifreq = 0; ifreq < NFreq2; ifreq++) {
		farComponent(etheta[ifreq], ephi[ifreq], e);
		for (int comp = 0; comp < 7; comp++) {
			ef[comp][ifreq] = e[comp];
		}
	}

	// title
	sprintf(str1, "far field (theta=%gdeg, phi=%gdeg)", theta, phi);
	strcpy(str2, (NFeed ? "[dB]" : "[dBsm]"));

	// plot
	plot2dFar0d0(
		NFreq2, ef,
		Far0dScale.user, Far0dScale.min, Far0dScale.max, Far0dScale.div,
		Freqdiv, Freq2[0], Freq2[NFreq2 - 1],
		Title, str1, str2,
		Width2d, Height2d, Fontsize2d);

	// free
	for (int comp = 0; comp < 7; comp++) {
		free(ef[comp]);
	}
	free(ef);
}


static void logFar0d(const d_complex_t etheta[], const d_complex_t ephi[], double theta, double phi)
{
	double e[7];

	FILE *fp;
	if ((fp = fopen(FN_far0d, "w")) == NULL) {
		fprintf(stderr, "*** %s open error.\n", FN_far0d);
		return;
	}

	fprintf(fp, "theta=%.3f[deg] phi=%.3f[deg]\n", theta, phi);
	fprintf(fp, "  No. frequency[Hz]    E-abs[dB]  E-theta[dB] E-theta[deg]    E-phi[dB]   E-phi[deg]  E-major[dB]  E-minor[dB]   E-RHCP[dB]   E-LHCP[dB] AxialRatio[dB]\n");
	for (int ifreq = 0; ifreq < NFreq2; ifreq++) {
		farComponent(etheta[ifreq], ephi[ifreq], e);
		// to dB
		for (int k = 0; k < 7; k++) {
			e[k] = 20 * log10(MAX(e[k], EPS2));
		}
		// output
		fprintf(fp, "%4d%15.5e%13.4f%13.4f%13.4f%13.4f%13.4f%13.4f%13.4f%13.4f%13.4f%13.4f\n",
			ifreq + 1, Freq2[ifreq], e[0], e[1], d_deg(etheta[ifreq]), e[2], d_deg(ephi[ifreq]), e[3], e[4], e[5], e[6], e[3] - e[4]);
	}

	fclose(fp);
}


void outputFar0d(void)
{
	const double theta = Far0d[0];
	const double phi   = Far0d[1];

	// alloc
	double *ffctr       =      (double *)malloc(NFreq2 * sizeof(double));
	d_complex_t *etheta = (d_complex_t *)malloc(NFreq2 * sizeof(d_complex_t));
	d_complex_t *ephi   = (d_complex_t *)malloc(NFreq2 * sizeof(d_complex_t));

	// setup
	setup_farfield();

	// factor
	for (int ifreq = 0; ifreq < NFreq2; ifreq++) {
		ffctr[ifreq] = farfactor(ifreq);
	}

	// calculation
	calcFar0d(ffctr, etheta, ephi, theta, phi);

	// plot
	plot2dFar0d(etheta, ephi, theta, phi);

	// write log
	logFar0d(etheta, ephi, theta, phi);

	// free
	free(ffctr);
	free(etheta);
	free(ephi);
}
