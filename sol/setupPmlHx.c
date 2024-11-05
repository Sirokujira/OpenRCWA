/*
setupPmlHx.c

setup PML for Hx
*/

#include "orcwa.h"

void setupPmlHx(int mode)
{
	const int lx = cPML.l;
	const int ly = cPML.l;
	const int lz = cPML.l;

	int64_t num = 0;
	for (int i = iMin - lx + 1; i < iMax + lx; i++) {
	for (int j = jMin - ly + 0; j < jMax + ly; j++) {
	for (int k = kMin - lz + 0; k < kMax + lz; k++) {
		if ((i < 0) || (Nx <  i) ||
		    (j < 0) || (Ny <= j) ||
		    (k < 0) || (Nz <= k)) {
			if      (mode == 1) {
				fPmlHx[num].i = i;
				fPmlHx[num].j = j;
				fPmlHx[num].k = k;
				const int i_ = MAX(0, MIN(Nx,     i));
				const int j_ = MAX(0, MIN(Ny - 1, j));
				const int k_ = MAX(0, MIN(Nz - 1, k));
				id_t m = 0;
				if      (i  <  0) m = IHX(0,      j_,     k_    );
				else if (Nx <  i) m = IHX(Nx,     j_,     k_    );
				else if (j  <  0) m = IHX(i_,     0,      k_    );
				else if (Ny <= j) m = IHX(i_,     Ny - 1, k_    );
				else if (k  <  0) m = IHX(i_,     j_,     0     );
				else if (Nz <= k) m = IHX(i_,     j_,     Nz - 1);
				fPmlHx[num].m = m;
			}
			num++;
		}
	}
	}
	}

	// array size
	if (mode == 0) {
		numPmlHx = num;
	}
}
