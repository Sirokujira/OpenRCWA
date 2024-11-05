/*
initfield.cu (CUDA)

initialize E and H
*/

#include "orcwa.h"
#include "orcwa_cuda.h"
#include "orcwa_prototype.h"

void initfield(void)
{
	size_t size, xsize, ysize, zsize;

	size = NN * sizeof(real_t);
	cuda_memset(GPU, Ex, 0, size);
	cuda_memset(GPU, Ey, 0, size);
	cuda_memset(GPU, Ez, 0, size);
	cuda_memset(GPU, Hx, 0, size);
	cuda_memset(GPU, Hy, 0, size);
	cuda_memset(GPU, Hz, 0, size);

	if (GPU) {
		size = sumGrid.x * sumGrid.y * sumGrid.z * sizeof(real_t);
		cuda_memset(GPU, d_sumE, 0, size);
		cuda_memset(GPU, d_sumH, 0, size);
	}

	if      (iABC == 1) {
		xsize = numPmlEx * sizeof(real_t);
		ysize = numPmlEy * sizeof(real_t);
		zsize = numPmlEz * sizeof(real_t);
		cuda_memset(GPU, Exy, 0, xsize);
		cuda_memset(GPU, Exz, 0, xsize);
		cuda_memset(GPU, Eyz, 0, ysize);
		cuda_memset(GPU, Eyx, 0, ysize);
		cuda_memset(GPU, Ezx, 0, zsize);
		cuda_memset(GPU, Ezy, 0, zsize);

		xsize = numPmlHx * sizeof(real_t);
		ysize = numPmlHy * sizeof(real_t);
		zsize = numPmlHz * sizeof(real_t);
		cuda_memset(GPU, Hxy, 0, xsize);
		cuda_memset(GPU, Hxz, 0, xsize);
		cuda_memset(GPU, Hyz, 0, ysize);
		cuda_memset(GPU, Hyx, 0, ysize);
		cuda_memset(GPU, Hzx, 0, zsize);
		cuda_memset(GPU, Hzy, 0, zsize);
	}

	if (NFeed > 0) {
		if (GPU) {
			cuda_memset(GPU, d_VFeed, 0, Feed_size);
			cuda_memset(GPU, d_IFeed, 0, Feed_size);
		}
		else {
			memset(VFeed, 0, Feed_size);
			memset(IFeed, 0, Feed_size);
		}
	}

	if (NPoint > 0) {
		if (GPU) {
			cuda_memset(GPU, d_VPoint, 0, Point_size);
		}
		else {
			memset(VPoint, 0, Point_size);
		}
	}
}
