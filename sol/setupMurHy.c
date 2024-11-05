/*
setupMurHy.c

setup Mur for Hy
*/

#include "orcwa.h"
#include "orcwa_prototype.h"

void setupMurHy(int mode)
{
	int64_t num = 0;
	for (int i = iMin - 1; i <= iMax; i++) {
	for (int j = jMin - 0; j <= jMax; j++) {
	for (int k = kMin - 1; k <= kMax; k++) {
		if ((((k < 0) || (k >= Nz)) && (i >= 0) && (i < Nx)) ||
		    (((i < 0) || (i >= Nx)) && (k >= 0) && (k < Nz))) {
			if      (mode == 1) {
				fMurHy[num].i = i;
				fMurHy[num].j = j;
				fMurHy[num].k = k;
				id_t   m = 0;
				double d = 0;
				int    i1 = 0, j1 = 0, k1 = 0;
				if      (k <   0) {
					m = IEX(i,     j,     k + 1);
					d = Zn[1] - Zn[0];
					i1 = i;
					j1 = j;
					k1 = k + 1;
				}
				else if (k >= Nz) {
					m = IEX(i,     j,     k    );
					d = Zn[Nz] - Zn[Nz - 1];
					i1 = i;
					j1 = j;
					k1 = k - 1;
				}
				else if (i <   0) {
					m = IEZ(i + 1, j,     k    );
					d = Xn[1] - Xn[0];
					i1 = i + 1;
					j1 = j;
					k1 = k;
				}
				else if (i >= Nx) {
					m = IEZ(i,     j,     k    );
					d = Xn[Nx] - Xn[Nx - 1];
					i1 = i - 1;
					j1 = j;
					k1 = k;
				}
				fMurHy[num].g = (real_t)factorMur(d, m);
				fMurHy[num].i1 = i1;
				fMurHy[num].j1 = j1;
				fMurHy[num].k1 = k1;
			}
			num++;
		}
	}
	}
	}

	// array size
	if (mode == 0) {
		numMurHy = num;
	}
}
