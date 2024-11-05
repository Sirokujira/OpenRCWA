/*
plot3dFar2d.c

plot far2d field (3D)
*/

#include "orcwa.h"
#include "complex.h"
#include "orcwa_prototype.h"

static void statistics(d_complex_t **, d_complex_t **, char [], char []);


void plot3dFar2d(d_complex_t ***etheta, d_complex_t ***ephi)
{
	const char scomp[][BUFSIZ] = {"E-abs", "E-theta", "E-phi", "E-major", "E-minor", "E-RHCP", "E-LHCP"};
	char unit[BUFSIZ], stat1[BUFSIZ], stat2[BUFSIZ];
	double e[7];

	strcpy(unit, (Far2dScale.db ? "[dB]" : ""));

	// alloc
	double **ef = (double **)malloc((Far2d.divtheta + 1) * sizeof(double *));
	for (int itheta = 0; itheta <= Far2d.divtheta; itheta++) {
		ef[itheta] = (double *)malloc((Far2d.divphi + 1) * sizeof(double));
	}
	char *comment[5];
	for (int n = 0; n < 5; n++) {
		comment[n] = (char *)malloc(BUFSIZ * sizeof(char));
	}

	for (int ifreq = 0; ifreq < NFreq2; ifreq++) {
		// calculate statistics
		statistics(etheta[ifreq], ephi[ifreq], stat1, stat2);

		for (int icomp = 0; icomp < 7; icomp++) {
			if (Far2dComp[icomp]) {
				double emax = 0;
				double tmax = 0;
				double pmax = 0;
				for (int itheta = 0; itheta <= Far2d.divtheta; itheta++) {
				for (int iphi   = 0; iphi   <= Far2d.divphi;   iphi++  ) {
					farComponent(etheta[ifreq][itheta][iphi], ephi[ifreq][itheta][iphi], e);
					ef[itheta][iphi] = e[icomp];
					if (e[icomp] > emax) {
						emax = e[icomp];
						tmax = 180.0 * itheta / Far2d.divtheta;
						pmax = 360.0 * iphi   / Far2d.divphi;
					}
				}
				}

				// comment
				strcpy(comment[0], Title);
				sprintf(comment[1], "%s%s f=%.3e[Hz]", scomp[icomp], unit, Freq2[ifreq]);
				strcpy(comment[2], stat1);
				strcpy(comment[3], stat2);
				sprintf(comment[4], "max[deg] @ theta=%.1f, phi=%.1f", tmax, pmax);

				// plot
				plot3dFar2d0(
					Far2d.divtheta, Far2d.divphi, ef,
					Far2dScale.db, Far2dScale.user, Far2dScale.min, Far2dScale.max,
					(int)NGline, Gline, Far2dObj,
					5, comment);
			}
		}
	}

	// free
	for (int itheta = 0; itheta <= Far2d.divtheta; itheta++) {
		free(ef[itheta]);
	}
	free(ef);
	for (int n = 0; n < 5; n++) {
		free(comment[n]);
	}
}


static void statistics(d_complex_t **etheta, d_complex_t **ephi, char stat1[], char stat2[])
{
	const double eps = 1e-20;

	double pmax = 0;
	double sum = 0;
	const double sumfactor = (PI / Far2d.divtheta) * (2 * PI / Far2d.divphi) / (4 * PI);
	for (int itheta = 0; itheta <= Far2d.divtheta; itheta++) {
	for (int iphi   = 0; iphi   <  Far2d.divphi;   iphi++  ) {
		const double pow = d_norm(etheta[itheta][iphi])
		                 + d_norm(  ephi[itheta][iphi]);
		sum += sumfactor * sin(PI * itheta / Far2d.divtheta) * pow;
		pmax = MAX(pow, pmax);
	}
	}

	if (NFeed) {
		double gain = pmax / sum;
		if (Far2dScale.db) {
			sprintf(stat1, "directive gain = %.3f[dBi]", 10 * log10(MAX(gain, eps)));
		}
		else {
			sprintf(stat1, "directive gain = %.3f", gain);
		}
		sprintf(stat2, "efficiency = %.3f[%%]", sum * 100);
	}
	else if (IPlanewave) {
		if (Far2dScale.db) {
			sprintf(stat1, "total cross section = %.3f[dBsm]", 10 * log10(MAX(sum, eps)));
			sprintf(stat2, "max = %.3f[dBsm]", 10 * log10(MAX(pmax, eps)));
		}
		else {
			sprintf(stat1, "total cross section = %.3e[sm]", sum);
			sprintf(stat2, "max = %.3e[sm]", pmax);
		}
	}
}
