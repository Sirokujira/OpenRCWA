/*
setupNear.c

setup near field factor
*/

#include "orcwa.h"
#include "complex.h"
#include "orcwa_prototype.h"
#include "finc.h"

void setupNear(void)
{
	if (NFreq2 <= 0) return;
	if (Solver.maxiter <= 0) return;

	Fnorm = (d_complex_t *)malloc(NFreq2 * sizeof(d_complex_t));

	// normalizing factor
	for (int ifreq = 0; ifreq < NFreq2; ifreq++) {
		double omega = 2 * PI * Freq2[ifreq];
		d_complex_t csum = d_complex(0, 0);
		for (int itime = 0; itime <= Solver.maxiter; itime++) {
			double t = itime * Dt;
			double fi = 0;
			if (NFeed) {
				fi = vfeed(t, Tw, 0);
			}
			else if (IPlanewave) {
				const double x0 = (Xn[0] + Xn[Nx]) / 2;
				const double y0 = (Yn[0] + Yn[Ny]) / 2;
				const double z0 = (Zn[0] + Zn[Nz]) / 2;
				real_t fi_, dfi;
				finc(x0, y0, z0, t, Planewave.r0, Planewave.ri, 1, Planewave.ai, Dt, &fi_, &dfi);
				fi = fi_;
			}
			const double phase = omega * t;
			csum = d_add(csum, d_rmul(fi, d_exp(-phase)));
		}
		Fnorm[ifreq] = csum;
	}

	// DFT factor
	for (int itime = 0; itime <= Solver.maxiter; itime++) {
		for (int ifreq = 0; ifreq < NFreq2; ifreq++) {
			const double omega = 2 * PI * Freq2[ifreq];
			const double hphase = omega * (itime + 0.5) * Dt;
			const double ephase = omega * (itime + 1.0) * Dt;
			const int id = (itime * NFreq2) + ifreq;
			cHdft[id] = d_div(d_exp(-hphase), Fnorm[ifreq]);
			cEdft[id] = d_div(d_exp(-ephase), Fnorm[ifreq]);
		}
	}
}
