/*
calcNear2d.c

calculate near2d field
*/

#include "orcwa.h"
#include "orcwa_prototype.h"


void calcNear2d(void)
{
	// alloc
	int *div1 = (int *)malloc(NNear2d * sizeof(int));
	int *div2 = (int *)malloc(NNear2d * sizeof(int));

	// setup
	for (int n = 0; n < NNear2d; n++) {
		if      (Near2d[n].dir == 'X') {
			// Y-Z
			div1[n] = Ny;
			div2[n] = Nz;
		}
		else if (Near2d[n].dir == 'Y') {
			// X-Z
			div1[n] = Nx;
			div2[n] = Nz;
		}
		else if (Near2d[n].dir == 'Z') {
			// X-Y
			div1[n] = Nx;
			div2[n] = Ny;
		}
	}

	// alloc
	int64_t num = 0;
	for (int n = 0; n < NNear2d; n++) {
		num += (div1[n] + 1) * (div2[n] + 1);
	}
	const size_t size = num * NFreq2 * sizeof(d_complex_t);
	Near2dEx = (d_complex_t *)malloc(size);
	Near2dEy = (d_complex_t *)malloc(size);
	Near2dEz = (d_complex_t *)malloc(size);
	Near2dHx = (d_complex_t *)malloc(size);
	Near2dHy = (d_complex_t *)malloc(size);
	Near2dHz = (d_complex_t *)malloc(size);

	// set
	int64_t adr = 0;
	for (int n = 0; n < NNear2d; n++) {
		for (int ifreq = 0; ifreq < NFreq2; ifreq++) {
			for (int l1 = 0; l1 <= div1[n]; l1++) {
			for (int l2 = 0; l2 <= div2[n]; l2++) {
				int i = 0, j = 0, k = 0;
				if      (Near2d[n].dir == 'X') {
					// Y-Z
					i = Near2d[n].id0;
					j = l1;
					k = l2;
				}
				else if (Near2d[n].dir == 'Y') {
					// X-Z
					j = Near2d[n].id0;
					k = l2;
					i = l1;
				}
				else if (Near2d[n].dir == 'Z') {
					// X-Y
					k = Near2d[n].id0;
					i = l1;
					j = l2;
				}
				NodeE_c(ifreq, i, j, k, &Near2dEx[adr], &Near2dEy[adr], &Near2dEz[adr]);
				NodeH_c(ifreq, i, j, k, &Near2dHx[adr], &Near2dHy[adr], &Near2dHz[adr]);
				adr++;
			}
			}
		}
	}
	//assert(adr * sizeof(complex_t) == size);

	// free
	free(div1);
	free(div2);
}
