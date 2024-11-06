/*
outputFar2d.c

output far2d field
*/

#include "orcwa.h"
#include "complex.h"
#include "orcwa_prototype.h"

static void logFar2d(d_complex_t ***etheta, d_complex_t ***ephi)
{
	double e[7];

	FILE *fp;
	if ((fp = fopen(FN_far2d, "w")) == NULL) {
		fprintf(stderr, "*** %s open error.\n", FN_far2d);
		return;
	}

	for (int ifreq = 0; ifreq < NFreq2; ifreq++) {
		fprintf(fp, "frequency[Hz] = %.5e\n", Freq2[ifreq]);
		fprintf(fp, " No. No. theta[deg] phi[deg]   E-abs[dB]  E-theta[dB] E-theta[deg]    E-phi[dB]   E-phi[deg]  E-major[dB]  E-minor[dB]   E-RHCP[dB]   E-LHCP[dB] AxialRatio[dB]\n");
		for (int itheta = 0; itheta <= Far2d.divtheta; itheta++) {
		for (int iphi   = 0; iphi   <= Far2d.divphi;   iphi++  ) {
			farComponent(etheta[ifreq][itheta][iphi], ephi[ifreq][itheta][iphi], e);
			// to dB
			for (int k = 0; k < 7; k++) {
				e[k] = 20 * log10(MAX(e[k], EPS2));
			}
			// output
			double theta = (180.0 * itheta) / Far2d.divtheta;
			double phi   = (360.0 * iphi)   / Far2d.divphi;
			fprintf(fp, "%4d%4d %9.1f%9.1f%13.4f%13.4f%13.4f%13.4f%13.4f%13.4f%13.4f%13.4f%13.4f%13.4f\n",
				itheta, iphi, theta, phi, e[0], e[1], d_deg(etheta[ifreq][itheta][iphi]), e[2], d_deg(ephi[ifreq][itheta][iphi]), e[3], e[4], e[5], e[6], e[3] - e[4]);
		}
		}
	}

	fclose(fp);
}


void outputFar2d(void)
{
	// alloc
	d_complex_t ***etheta = (d_complex_t ***)malloc(NFreq2 * sizeof(d_complex_t **));
	d_complex_t ***ephi   = (d_complex_t ***)malloc(NFreq2 * sizeof(d_complex_t **));
	for (int ifreq = 0; ifreq < NFreq2; ifreq++) {
		etheta[ifreq] = (d_complex_t **)malloc((Far2d.divtheta + 1) * sizeof(d_complex_t *));
		ephi[ifreq]   = (d_complex_t **)malloc((Far2d.divtheta + 1) * sizeof(d_complex_t *));
		for (int itheta = 0; itheta <= Far2d.divtheta; itheta++) {
			etheta[ifreq][itheta] = (d_complex_t *)malloc((Far2d.divphi + 1) * sizeof(d_complex_t));
			ephi[ifreq][itheta]   = (d_complex_t *)malloc((Far2d.divphi + 1) * sizeof(d_complex_t));
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
	calcFar2d(ffctr, etheta, ephi);

	// plot
	plot3dFar2d(etheta, ephi);

	// write log
	logFar2d(etheta, ephi);

	// free
	for (int ifreq = 0; ifreq < NFreq2; ifreq++) {
		for (int itheta = 0; itheta <= Far2d.divtheta; itheta++) {
			free(etheta[ifreq][itheta]);
			free(ephi[ifreq][itheta]);
		}
		free(etheta[ifreq]);
		free(ephi[ifreq]);
	}
	free(etheta);
	free(ephi);

	free(ffctr);
}
