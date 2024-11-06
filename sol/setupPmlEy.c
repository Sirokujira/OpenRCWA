/*
setupPmlEy.c

setup PML for Ey
*/

#include "orcwa.h"

void setupPmlEy(int mode)
{
	const int lx = cPML.l;
	const int ly = cPML.l;
	const int lz = cPML.l;

	int64_t num = 0;
	for (int i = iMin - lx + 1; i < iMax + lx; i++) {
	for (int j = jMin - ly + 0; j < jMax + ly; j++) {
	for (int k = kMin - lz + 1; k < kMax + lz; k++) {
		if ((i < 0) || (Nx <  i) ||
		    (j < 0) || (Ny <= j) ||
		    (k < 0) || (Nz <  k)) {
			if      (mode == 1) {
				fPmlEy[num].i = i;
				fPmlEy[num].j = j;
				fPmlEy[num].k = k;
				const int i_ = MAX(0, MIN(Nx,     i));
				const int j_ = MAX(0, MIN(Ny - 1, j));
				const int k_ = MAX(0, MIN(Nz,     k));
				id_t m = 0;
				if      (j  <  0) m = IEY(i_,     0,      k_    );
				else if (Ny <= j) m = IEY(i_,     Ny - 1, k_    );
				else if (k  <  0) m = IEY(i_,     j_,     0     );
				else if (Nz <  k) m = IEY(i_,     j_,     Nz    );
				else if (i  <  0) m = IEY(0,      j_,     k_    );
				else if (Nx <  i) m = IEY(Nx,     j_,     k_    );
				fPmlEy[num].m = m;
			}
			num++;
		}
	}
	}
	}

	// array size
	if (mode == 0) {
		numPmlEy = num;
	}
}
