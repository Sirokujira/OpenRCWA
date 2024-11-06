/*
plot2dFreq.c

plot frequency char.s (2D)
*/

#include "orcwa.h"
#include "complex.h"
#include "orcwa_prototype.h"

void plot2dFreq(void)
{
	const int nfreq = NFreq1;
	const double *freq = Freq1;
	const int nfeed = NFeed;
	const int npoint = NPoint;
	int id;

	// alloc
	double *z0 = NULL;
	if (nfeed) {
		z0 = (double *)malloc(nfeed * sizeof(double));
		for (int ifeed = 0; ifeed < nfeed; ifeed++) {
			z0[ifeed] = Feed[ifeed].z0;
		}
	}

	// (1) Smith Chart

	id = 0;
	if (IFreq[id] && nfeed && nfreq) {
		plot2dSmith(nfeed, nfreq, Zin, z0, freq, Title, Width2d, Height2d, Fontsize2d);
	}

	// (2) Zin

	id = 1;
	if (IFreq[id] && nfeed && nfreq) {
		plot2dZin(nfeed, nfreq, Zin, FreqScale[id], Freqdiv, freq, Title, Width2d, Height2d, Fontsize2d);
	}

	// (3) Yin

	id = 2;
	if (IFreq[id] && nfeed && nfreq) {
		plot2dYin(nfeed, nfreq, Zin, FreqScale[id], Freqdiv, freq, Title, Width2d, Height2d, Fontsize2d);
	}

	// (4) reflection

	id = 3;
	if (IFreq[id] && nfeed && nfreq) {
		plot2dRef(nfeed, nfreq, Zin, z0, FreqScale[id], Freqdiv, freq, Title, Width2d, Height2d, Fontsize2d);
	}

	// (5) S-parameter

	id = 4;
	if (IFreq[id] && npoint && nfreq) {
		plot2dSpara(npoint, nfreq, Spara, FreqScale[id], Freqdiv, freq, Title, Width2d, Height2d, Fontsize2d);
	}

	// (6) coupling

	id = 5;
	if (IFreq[id] && nfeed && npoint && nfreq) {
		double ***couple = (double ***)malloc(nfeed * npoint * nfreq * sizeof(double **));
		for (int ifeed = 0; ifeed < nfeed; ifeed++) {
			couple[ifeed] = (double **)malloc(npoint * nfreq * sizeof(double *));
			for (int ipoint = 0; ipoint < npoint; ipoint++) {
				couple[ifeed][ipoint] = (double *)malloc(nfreq * sizeof(double));
				for (int ifreq = 0; ifreq < nfreq; ifreq++) {
					couple[ifeed][ipoint][ifreq] = d_abs(coupling(ifeed, ipoint, ifreq));
				}
			}
		}

		plot2dCoupling(nfeed, npoint, nfreq, couple, FreqScale[id], Freqdiv, freq, Title, Width2d, Height2d, Fontsize2d);

		for (int ifeed = 0; ifeed < nfeed; ifeed++) {
			for (int ipoint = 0; ipoint < npoint; ipoint++) {
				free(couple[ifeed][ipoint]);
			}
			free(couple[ifeed]);
		}
		free(couple);
	}

	// free
	if (nfeed) {
		free(z0);
	}
}
