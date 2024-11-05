/*
setupMurHx.c

setup Mur for Hx
*/

#include "orcwa.h"
#include "orcwa_prototype.h"

void setupMurHx(int mode)
{
	int64_t num = 0;
	for (int i = iMin - 0; i <= iMax; i++) {
	for (int j = jMin - 1; j <= jMax; j++) {
	for (int k = kMin - 1; k <= kMax; k++) {
		if ((((j < 0) || (j >= Ny)) && (k >= 0) && (k < Nz)) ||
		    (((k < 0) || (k >= Nz)) && (j >= 0) && (j < Ny))) {
			if      (mode == 1) {
				fMurHx[num].i = i;
				fMurHx[num].j = j;
				fMurHx[num].k = k;
				id_t   m = 0;
				double d = 0;
				int    i1 = 0, j1 = 0, k1 = 0;
				if      (j <   0) {
					m = IEZ(i,     j + 1, k    );
					d = Yn[1] - Yn[0];
					i1 = i;
					j1 = j + 1;
					k1 = k;
				}
				else if (j >= Ny) {
					m = IEZ(i,     j,     k    );
					d = Yn[Ny] - Yn[Ny - 1];
					i1 = i;
					j1 = j - 1;
					k1 = k;
				}
				else if (k <   0) {
					m = IEY(i,     j,     k + 1);
					d = Zn[1] - Zn[0];
					i1 = i;
					j1 = j;
					k1 = k + 1;
				}
				else if (k >= Nz) {
					m = IEY(i,     j,     k    );
					d = Zn[Nz] - Zn[Nz - 1];
					i1 = i;
					j1 = j;
					k1 = k - 1;
				}
				fMurHx[num].g = (real_t)factorMur(d, m);
				fMurHx[num].i1 = i1;
				fMurHx[num].j1 = j1;
				fMurHx[num].k1 = k1;
			}
			num++;
		}
	}
	}
	}

	// array size
	if (mode == 0) {
		numMurHx = num;
	}
}
