/*
memallocfree2_gpu.cu (CUDA)

alloc and free
(2) iteration variables
*/

#include "orcwa.h"
#include "orcwa_cuda.h"
#include "orcwa_prototype.h"


void memalloc2_gpu()
{
	size_t size, xsize, ysize, zsize;

	// E/H
	size = NN * sizeof(real_t);
	cuda_malloc(GPU, UM, (void **)&Ex, size);
	cuda_malloc(GPU, UM, (void **)&Ey, size);
	cuda_malloc(GPU, UM, (void **)&Ez, size);
	cuda_malloc(GPU, UM, (void **)&Hx, size);
	cuda_malloc(GPU, UM, (void **)&Hy, size);
	cuda_malloc(GPU, UM, (void **)&Hz, size);

	// ABC
	if      (iABC == 1) {
		xsize = numPmlEx * sizeof(real_t);
		ysize = numPmlEy * sizeof(real_t);
		zsize = numPmlEz * sizeof(real_t);
		cuda_malloc(GPU, UM, (void **)&Exy, xsize);
		cuda_malloc(GPU, UM, (void **)&Exz, xsize);
		cuda_malloc(GPU, UM, (void **)&Eyz, ysize);
		cuda_malloc(GPU, UM, (void **)&Eyx, ysize);
		cuda_malloc(GPU, UM, (void **)&Ezx, zsize);
		cuda_malloc(GPU, UM, (void **)&Ezy, zsize);

		xsize = numPmlHx * sizeof(real_t);
		ysize = numPmlHy * sizeof(real_t);
		zsize = numPmlHz * sizeof(real_t);
		cuda_malloc(GPU, UM, (void **)&Hxy, xsize);
		cuda_malloc(GPU, UM, (void **)&Hxz, xsize);
		cuda_malloc(GPU, UM, (void **)&Hyz, ysize);
		cuda_malloc(GPU, UM, (void **)&Hyx, ysize);
		cuda_malloc(GPU, UM, (void **)&Hzx, zsize);
		cuda_malloc(GPU, UM, (void **)&Hzy, zsize);
	}

	// feed
	if (NFeed > 0) {
		size = NFeed * sizeof(feed_t);
		cuda_malloc(GPU, UM, (void **)&d_Feed, size);
		cuda_memcpy(GPU, d_Feed, Feed, size, cudaMemcpyHostToDevice);

		cuda_malloc(GPU, UM, (void **)&d_VFeed, Feed_size);
		cuda_malloc(GPU, UM, (void **)&d_IFeed, Feed_size);
	}

	// inductor
	if (NInductor > 0) {
		size = NInductor * sizeof(inductor_t);
		cuda_malloc(GPU, UM, (void **)&d_Inductor, size);
		cuda_memcpy(GPU, d_Inductor, Inductor, size, cudaMemcpyHostToDevice);
	}

	// point
	if (NPoint > 0) {
		size = (NPoint + 2) * sizeof(point_t);
		cuda_malloc(GPU, UM, (void **)&d_Point, size);
		cuda_memcpy(GPU, d_Point, Point, size, cudaMemcpyHostToDevice);

		cuda_malloc(GPU, UM, (void **)&d_VPoint, Point_size);
	}
}


void memfree2_gpu()
{
	cuda_free(GPU, Ex);
	cuda_free(GPU, Ey);
	cuda_free(GPU, Ez);
	cuda_free(GPU, Hx);
	cuda_free(GPU, Hy);
	cuda_free(GPU, Hz);

	cuda_free(GPU, d_iEx);
	cuda_free(GPU, d_iEy);
	cuda_free(GPU, d_iEz);
	cuda_free(GPU, d_iHx);
	cuda_free(GPU, d_iHy);
	cuda_free(GPU, d_iHz);

	cuda_free(GPU, d_C1);
	cuda_free(GPU, d_C2);
	cuda_free(GPU, d_D1);
	cuda_free(GPU, d_D2);

	if (VECTOR) {
		cuda_free(GPU, d_K1Ex);
		cuda_free(GPU, d_K2Ex);
		cuda_free(GPU, d_K1Ey);
		cuda_free(GPU, d_K2Ey);
		cuda_free(GPU, d_K1Ez);
		cuda_free(GPU, d_K2Ez);
		cuda_free(GPU, d_K1Hx);
		cuda_free(GPU, d_K2Hx);
		cuda_free(GPU, d_K1Hy);
		cuda_free(GPU, d_K2Hy);
		cuda_free(GPU, d_K1Hz);
		cuda_free(GPU, d_K2Hz);
	}

	if      (iABC == 0) {
		cuda_free(GPU, d_fMurHx);
		cuda_free(GPU, d_fMurHy);
		cuda_free(GPU, d_fMurHz);
	}
	else if (iABC == 1) {
		cuda_free(GPU, Exy);
		cuda_free(GPU, Exz);
		cuda_free(GPU, Eyz);
		cuda_free(GPU, Eyx);
		cuda_free(GPU, Ezx);
		cuda_free(GPU, Ezy);

		cuda_free(GPU, Hxy);
		cuda_free(GPU, Hxz);
		cuda_free(GPU, Hyz);
		cuda_free(GPU, Hyx);
		cuda_free(GPU, Hzx);
		cuda_free(GPU, Hzy);

		cuda_free(GPU, d_fPmlEx);
		cuda_free(GPU, d_fPmlEy);
		cuda_free(GPU, d_fPmlEz);

		cuda_free(GPU, d_fPmlHx);
		cuda_free(GPU, d_fPmlHy);
		cuda_free(GPU, d_fPmlHz);

		cuda_free(GPU, d_gPmlXn);
		cuda_free(GPU, d_gPmlYn);
		cuda_free(GPU, d_gPmlZn);

		cuda_free(GPU, d_gPmlXc);
		cuda_free(GPU, d_gPmlYc);
		cuda_free(GPU, d_gPmlZc);

		cuda_free(GPU, d_rPmlE);
		cuda_free(GPU, d_rPmlH);
		cuda_free(GPU, d_rPml);
	}

	cuda_free(GPU, d_Xn);
	cuda_free(GPU, d_Yn);
	cuda_free(GPU, d_Zn);

	cuda_free(GPU, d_Xc);
	cuda_free(GPU, d_Yc);
	cuda_free(GPU, d_Zc);

	cuda_free(GPU, d_RXn);
	cuda_free(GPU, d_RYn);
	cuda_free(GPU, d_RZn);

	cuda_free(GPU, d_RXc);
	cuda_free(GPU, d_RYc);
	cuda_free(GPU, d_RZc);

	if (NFeed > 0) {
		cuda_free(GPU, d_Feed);
		cuda_free(GPU, d_VFeed);
		cuda_free(GPU, d_IFeed);
	}

	if (NPoint > 0) {
		cuda_free(GPU, d_Point);
		cuda_free(GPU, d_VPoint);
	}

	// host memory (cf. setup_host())

	free(h_Xn);
	free(h_Yn);
	free(h_Zn);

	free(h_Xc);
	free(h_Yc);
	free(h_Zc);
}
