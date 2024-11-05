/*
calcNear1d.c

calculate near1d field
*/

#include "orcwa.h"
#include "orcwa_prototype.h"


void calcNear1d(void)
{
	// alloc
	int *div = (int *)malloc(NNear1d * sizeof(int));
	for (int n = 0; n < NNear1d; n++) {
		if      (Near1d[n].dir == 'X') {
			div[n] = Nx;
		}
		else if (Near1d[n].dir == 'Y') {
			div[n] = Ny;
		}
		else if (Near1d[n].dir == 'Z') {
			div[n] = Nz;
		}
	}

	int num = 0;
	for (int n = 0; n < NNear1d; n++) {
		num += div[n] + 1;
	}
	const size_t size = num * NFreq2 * sizeof(d_complex_t);
	Near1dEx = (d_complex_t *)malloc(size);
	Near1dEy = (d_complex_t *)malloc(size);
	Near1dEz = (d_complex_t *)malloc(size);
	Near1dHx = (d_complex_t *)malloc(size);
	Near1dHy = (d_complex_t *)malloc(size);
	Near1dHz = (d_complex_t *)malloc(size);

	int64_t adr = 0;
	for (int n = 0; n < NNear1d; n++) {
		for (int ifreq = 0; ifreq < NFreq2; ifreq++) {
			for (int l = 0; l <= div[n]; l++) {
				int i = 0, j = 0, k = 0;
				if      (Near1d[n].dir == 'X') {
					i = l;
					j = Near1d[n].id1;
					k = Near1d[n].id2;
				}
				else if (Near1d[n].dir == 'Y') {
					j = l;
					k = Near1d[n].id1;
					i = Near1d[n].id2;
				}
				else if (Near1d[n].dir == 'Z') {
					k = l;
					i = Near1d[n].id1;
					j = Near1d[n].id2;
				}
				NodeE_c(ifreq, i, j, k, &Near1dEx[adr], &Near1dEy[adr], &Near1dEz[adr]);
				NodeH_c(ifreq, i, j, k, &Near1dHx[adr], &Near1dHy[adr], &Near1dHz[adr]);
				//printf("%d %d %d %d %d %d %zd\n", n, ifreq, l, i, j, k, adr);
				//printf("%d %f %f\n", l, Near1dEz[adr].r, Near1dEz[adr].i);
				adr++;
			}
		}
	}
	//assert(adr * sizeof(d_complex_t) == size);

	// free
	free(div);
}
