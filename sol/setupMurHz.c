/*
setupMurHz.c

setup Mur for Hz
*/

#include "orcwa.h"
#include "orcwa_prototype.h"

void setupMurHz(int mode)
{
	int64_t num = 0;
	for (int i = iMin - 1; i <= iMax; i++) {
	for (int j = jMin - 1; j <= jMax; j++) {
	for (int k = kMin - 0; k <= kMax; k++) {
		if ((((i < 0) || (i >= Nx)) && (j >= 0) && (j < Ny)) ||
		    (((j < 0) || (j >= Ny)) && (i >= 0) && (i < Nx))) {
			if      (mode == 1) {
				fMurHz[num].i = i;
				fMurHz[num].j = j;
				fMurHz[num].k = k;
				id_t   m = 0;
				double d = 0;
				int    i1 = 0, j1 = 0, k1 = 0;
				if      (i <   0) {
					m = IEY(i + 1, j,     k    );
					d = Xn[1] - Xn[0];
					i1 = i + 1;
					j1 = j;
					k1 = k;
				}
				else if (i >= Nx) {
					m = IEY(i,     j,     k    );
					d = Xn[Nx] - Xn[Nx - 1];
					i1 = i - 1;
					j1 = j;
					k1 = k;
				}
				else if (j <   0) {
					m = IEX(i,     j + 1, k    );
					d = Yn[1] - Yn[0];
					i1 = i;
					j1 = j + 1;
					k1 = k;
				}
				else if (j >= Ny) {
					m = IEX(i,     j,     k    );
					d = Yn[Ny] - Yn[Ny - 1];
					i1 = i;
					j1 = j - 1;
					k1 = k;
				}
				fMurHz[num].g = (real_t)factorMur(d, m);
				fMurHz[num].i1 = i1;
				fMurHz[num].j1 = j1;
				fMurHz[num].k1 = k1;
			}
			num++;
		}
	}
	}
	}

	// array size
	if (mode == 0) {
		numMurHz = num;
	}
}
