/*
calcFar1d.c

calculate far1d field
*/

#include "orcwa.h"
#include "complex.h"
#include "orcwa_prototype.h"

void calcFar1d(const double ffctr[], d_complex_t ***etheta, d_complex_t ***ephi)
{
	for (int n = 0; n < NFar1d; n++) {
		for (int ifreq = 0; ifreq < NFreq2; ifreq++) {
			int i;
#ifdef _OPENMP
#pragma omp parallel for
#endif
			for (i = 0; i <= Far1d[n].div; i++) {
				// direction
				double angle = (360.0 * i) / Far1d[n].div;
				double theta = 0;
				double phi = 0;
				if      (Far1d[n].dir == 'X') {
					theta = angle;
					phi   = 90;
				}
				else if (Far1d[n].dir == 'Y') {
					theta = angle;
					phi   = 0;
				}
				else if (Far1d[n].dir == 'Z') {
					theta = 90;
					phi   = angle;
				}
				else if (Far1d[n].dir == 'V') {
					theta = angle;
					phi   = Far1d[n].angle;
				}
				else if (Far1d[n].dir == 'H') {
					theta = Far1d[n].angle;
					phi   = angle;
				}

				// far field
				farfield(ifreq, theta, phi, ffctr[ifreq], &etheta[n][ifreq][i], &ephi[n][ifreq][i]);
			}

		 	// normalization
			if (Far1dNorm) {
				double pmax = 0;
				for (i = 0; i <= Far1d[n].div; i++) {
					pmax = MAX(pmax, d_norm(etheta[n][ifreq][i]) + d_norm(ephi[n][ifreq][i]));
				}
				const double fctr = 1 / sqrt(pmax);
				for (i = 0; i <= Far1d[n].div; i++) {
					etheta[n][ifreq][i] = d_rmul(fctr, etheta[n][ifreq][i]);
					ephi[n][ifreq][i]   = d_rmul(fctr, ephi[n][ifreq][i]);
				}
			}
		}
	}
}
