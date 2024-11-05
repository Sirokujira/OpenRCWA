/*
initfield.c

initialize E and H
*/

#include "orcwa.h"

void initfield(void)
{
	size_t size, xsize, ysize, zsize;

	size = NN * sizeof(real_t);
	memset(Ex, 0, size);
	memset(Ey, 0, size);
	memset(Ez, 0, size);
	memset(Hx, 0, size);
	memset(Hy, 0, size);
	memset(Hz, 0, size);

	if      (iABC == 0) {
		for (int64_t n = 0; n < numMurHx; n++) {
			fMurHx[n].f = 0;
		}
		for (int64_t n = 0; n < numMurHy; n++) {
			fMurHy[n].f = 0;
		}
		for (int64_t n = 0; n < numMurHz; n++) {
			fMurHz[n].f = 0;
		}
	}
	else if (iABC == 1) {
		xsize = numPmlEx * sizeof(real_t);
		ysize = numPmlEy * sizeof(real_t);
		zsize = numPmlEz * sizeof(real_t);
		memset(Exy, 0, xsize);
		memset(Exz, 0, xsize);
		memset(Eyz, 0, ysize);
		memset(Eyx, 0, ysize);
		memset(Ezx, 0, zsize);
		memset(Ezy, 0, zsize);

		xsize = numPmlHx * sizeof(real_t);
		ysize = numPmlHy * sizeof(real_t);
		zsize = numPmlHz * sizeof(real_t);
		memset(Hxy, 0, xsize);
		memset(Hxz, 0, xsize);
		memset(Hyz, 0, ysize);
		memset(Hyx, 0, ysize);
		memset(Hzx, 0, zsize);
		memset(Hzy, 0, zsize);
	}

	memset(Eiter, 0, Iter_size);
	memset(Hiter, 0, Iter_size);
	Niter = 0;

	if (NFeed > 0) {
		memset(VFeed, 0, Feed_size);
		memset(IFeed, 0, Feed_size);
	}

	if (NPoint > 0) {
		memset(VPoint, 0, Point_size);
	}

	if ((NN > 0) && (NFreq2 > 0)) {
		size = NFreq2 * NN * sizeof(float);
		memset(cEx_r, 0, size);
		memset(cEx_i, 0, size);
		memset(cEy_r, 0, size);
		memset(cEy_i, 0, size);
		memset(cEz_r, 0, size);
		memset(cEz_i, 0, size);
		memset(cHx_r, 0, size);
		memset(cHx_i, 0, size);
		memset(cHy_r, 0, size);
		memset(cHy_i, 0, size);
		memset(cHz_r, 0, size);
		memset(cHz_i, 0, size);
	}
}
