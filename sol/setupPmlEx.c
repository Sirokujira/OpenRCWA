/*
setupPmlEx.c

setup PML for Ex
*/

#include "orcwa.h"

void setupPmlEx(int mode)
{
	const int lx = cPML.l;
	const int ly = cPML.l;
	const int lz = cPML.l;

	int64_t num = 0;
	for (int i = iMin - lx + 0; i < iMax + lx; i++) {
	for (int j = jMin - ly + 1; j < jMax + ly; j++) {
	for (int k = kMin - lz + 1; k < kMax + lz; k++) {
		if ((i < 0) || (Nx <= i) ||
		    (j < 0) || (Ny <  j) ||
		    (k < 0) || (Nz <  k)) {
			if      (mode == 1) {
				fPmlEx[num].i = i;
				fPmlEx[num].j = j;
				fPmlEx[num].k = k;
				const int i_ = MAX(0, MIN(Nx - 1, i));
				const int j_ = MAX(0, MIN(Ny,     j));
				const int k_ = MAX(0, MIN(Nz,     k));
				id_t m = 0;
				if      (i  <  0) m = IEX(0,      j_,     k_    );
				else if (Nx <= i) m = IEX(Nx - 1, j_,     k_    );
				else if (j  <  0) m = IEX(i_,     0,      k_    );
				else if (Ny <  j) m = IEX(i_,     Ny,     k_    );
				else if (k  <  0) m = IEX(i_,     j_,     0     );
				else if (Nz <  k) m = IEX(i_,     j_,     Nz    );
				fPmlEx[num].m = m;
			}
			num++;
		}
	}
	}
	}

	// array size
	if (mode == 0) {
		numPmlEx = num;
	}
}
