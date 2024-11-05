/*
calcFar2d.c

calculate far2d field
*/

#include "orcwa.h"
#include "orcwa_prototype.h"

void calcFar2d(const double ffctr[], d_complex_t ***etheta, d_complex_t ***ephi)
{
	for (int ifreq = 0; ifreq < NFreq2; ifreq++) {
		int itheta;
#ifdef _OPENMP
#pragma omp parallel for
#endif
		for (    itheta = 0; itheta <= Far2d.divtheta; itheta++) {
		for (int iphi   = 0; iphi   <= Far2d.divphi;   iphi++  ) {
			const double theta = 180.0 * itheta / Far2d.divtheta;
			const double phi   = 360.0 * iphi   / Far2d.divphi;
			farfield(ifreq, theta, phi, ffctr[ifreq], &etheta[ifreq][itheta][iphi], &ephi[ifreq][itheta][iphi]);
		}
		}
	}
}
