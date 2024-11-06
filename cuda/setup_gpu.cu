/*
setup_gpu.cu (CUDA)
*/

#include "orcwa.h"
#include "orcwa_cuda.h"
#include "orcwa_prototype.h"


void setup_gpu()
{
	size_t size, xsize, ysize, zsize;

	// execution configuration

	updateBlock = dim3(32, 4, 1);

	sumBlock = dim3(32, 4, 1);
	sumGrid = dim3(CEIL(kMax - kMin, sumBlock.x),
	               CEIL(jMax - jMin, sumBlock.y),
	               CEIL(iMax - iMin, sumBlock.z));

	dispersionBlock = 256;
	murBlock = 256;
	pmlBlock = 256;
	pbcBlock = 16;

	near1dBlock = 256;
	near2dBlock = dim3(32, 8);

	// average array

	size = sumGrid.x * sumGrid.y * sumGrid.z * sizeof(real_t);
	//printf("%d %d %d %lld\n", sumGrid.x, sumGrid.y, sumGrid.z, size);
	h_sumE = (real_t *)malloc(size);
	h_sumH = (real_t *)malloc(size);
	cuda_malloc(GPU, UM, (void **)&d_sumE, size);
	cuda_malloc(GPU, UM, (void **)&d_sumH, size);

	// === host/device memory ===

	// mesh

	xsize = (Nx + 1) * sizeof(real_t);
	ysize = (Ny + 1) * sizeof(real_t);
	zsize = (Nz + 1) * sizeof(real_t);

	cuda_malloc(GPU, UM, (void **)&d_Xn,  xsize);
	cuda_malloc(GPU, UM, (void **)&d_Yn,  ysize);
	cuda_malloc(GPU, UM, (void **)&d_Zn,  zsize);
	cuda_malloc(GPU, UM, (void **)&d_RXn, xsize);
	cuda_malloc(GPU, UM, (void **)&d_RYn, ysize);
	cuda_malloc(GPU, UM, (void **)&d_RZn, zsize);

	cuda_memcpy(GPU, d_Xn, h_Xn, xsize, cudaMemcpyHostToDevice);
	cuda_memcpy(GPU, d_Yn, h_Yn, ysize, cudaMemcpyHostToDevice);
	cuda_memcpy(GPU, d_Zn, h_Zn, zsize, cudaMemcpyHostToDevice);
	cuda_memcpy(GPU, d_RXn, RXn, xsize, cudaMemcpyHostToDevice);
	cuda_memcpy(GPU, d_RYn, RYn, ysize, cudaMemcpyHostToDevice);
	cuda_memcpy(GPU, d_RZn, RZn, zsize, cudaMemcpyHostToDevice);

	xsize = (Nx + 0) * sizeof(real_t);
	ysize = (Ny + 0) * sizeof(real_t);
	zsize = (Nz + 0) * sizeof(real_t);

	cuda_malloc(GPU, UM, (void **)&d_Xc,  xsize);
	cuda_malloc(GPU, UM, (void **)&d_Yc,  ysize);
	cuda_malloc(GPU, UM, (void **)&d_Zc,  zsize);
	cuda_malloc(GPU, UM, (void **)&d_RXc, xsize);
	cuda_malloc(GPU, UM, (void **)&d_RYc, ysize);
	cuda_malloc(GPU, UM, (void **)&d_RZc, zsize);

	cuda_memcpy(GPU, d_Xc, h_Xc, xsize, cudaMemcpyHostToDevice);
	cuda_memcpy(GPU, d_Yc, h_Yc, ysize, cudaMemcpyHostToDevice);
	cuda_memcpy(GPU, d_Zc, h_Zc, zsize, cudaMemcpyHostToDevice);
	cuda_memcpy(GPU, d_RXc, RXc, xsize, cudaMemcpyHostToDevice);
	cuda_memcpy(GPU, d_RYc, RYc, ysize, cudaMemcpyHostToDevice);
	cuda_memcpy(GPU, d_RZc, RZc, zsize, cudaMemcpyHostToDevice);

	// material ID

	size = NN * sizeof(id_t);
	cuda_malloc(GPU, UM, (void **)&d_iEx, size);
	cuda_malloc(GPU, UM, (void **)&d_iEy, size);
	cuda_malloc(GPU, UM, (void **)&d_iEz, size);
	cuda_malloc(GPU, UM, (void **)&d_iHx, size);
	cuda_malloc(GPU, UM, (void **)&d_iHy, size);
	cuda_malloc(GPU, UM, (void **)&d_iHz, size);
	cuda_memcpy(GPU, d_iEx, iEx, size, cudaMemcpyHostToDevice);
	cuda_memcpy(GPU, d_iEy, iEy, size, cudaMemcpyHostToDevice);
	cuda_memcpy(GPU, d_iEz, iEz, size, cudaMemcpyHostToDevice);
	cuda_memcpy(GPU, d_iHx, iHx, size, cudaMemcpyHostToDevice);
	cuda_memcpy(GPU, d_iHy, iHy, size, cudaMemcpyHostToDevice);
	cuda_memcpy(GPU, d_iHz, iHz, size, cudaMemcpyHostToDevice);

	// material factor

	size = NMaterial * sizeof(real_t);
	cuda_malloc(GPU, UM, (void **)&d_C1, size);
	cuda_malloc(GPU, UM, (void **)&d_C2, size);
	cuda_malloc(GPU, UM, (void **)&d_D1, size);
	cuda_malloc(GPU, UM, (void **)&d_D2, size);
	cuda_memcpy(GPU, d_C1, C1, size, cudaMemcpyHostToDevice);
	cuda_memcpy(GPU, d_C2, C2, size, cudaMemcpyHostToDevice);
	cuda_memcpy(GPU, d_D1, D1, size, cudaMemcpyHostToDevice);
	cuda_memcpy(GPU, d_D2, D2, size, cudaMemcpyHostToDevice);

	if (VECTOR) {
		size = NN * sizeof(real_t);
		cuda_malloc(GPU, UM, (void **)&d_K1Ex, size);
		cuda_malloc(GPU, UM, (void **)&d_K2Ex, size);
		cuda_malloc(GPU, UM, (void **)&d_K1Ey, size);
		cuda_malloc(GPU, UM, (void **)&d_K2Ey, size);
		cuda_malloc(GPU, UM, (void **)&d_K1Ez, size);
		cuda_malloc(GPU, UM, (void **)&d_K2Ez, size);
		cuda_malloc(GPU, UM, (void **)&d_K1Hx, size);
		cuda_malloc(GPU, UM, (void **)&d_K2Hx, size);
		cuda_malloc(GPU, UM, (void **)&d_K1Hy, size);
		cuda_malloc(GPU, UM, (void **)&d_K2Hy, size);
		cuda_malloc(GPU, UM, (void **)&d_K1Hz, size);
		cuda_malloc(GPU, UM, (void **)&d_K2Hz, size);
		cuda_memcpy(GPU, d_K1Ex, K1Ex, size, cudaMemcpyHostToDevice);
		cuda_memcpy(GPU, d_K2Ex, K2Ex, size, cudaMemcpyHostToDevice);
		cuda_memcpy(GPU, d_K1Ey, K1Ey, size, cudaMemcpyHostToDevice);
		cuda_memcpy(GPU, d_K2Ey, K2Ey, size, cudaMemcpyHostToDevice);
		cuda_memcpy(GPU, d_K1Ez, K1Ez, size, cudaMemcpyHostToDevice);
		cuda_memcpy(GPU, d_K2Ez, K2Ez, size, cudaMemcpyHostToDevice);
		cuda_memcpy(GPU, d_K1Hx, K1Hx, size, cudaMemcpyHostToDevice);
		cuda_memcpy(GPU, d_K2Hx, K2Hx, size, cudaMemcpyHostToDevice);
		cuda_memcpy(GPU, d_K1Hy, K1Hy, size, cudaMemcpyHostToDevice);
		cuda_memcpy(GPU, d_K2Hy, K2Hy, size, cudaMemcpyHostToDevice);
		cuda_memcpy(GPU, d_K1Hz, K1Hz, size, cudaMemcpyHostToDevice);
		cuda_memcpy(GPU, d_K2Hz, K2Hz, size, cudaMemcpyHostToDevice);
	}

	// disperson
	//printf("%zd %zd %zd\n", numDispersionEx, numDispersionEy, numDispersionEz);
	if (numDispersionEx > 0) {
		size = numDispersionEx * sizeof(dispersion_t);
		cuda_malloc(GPU, UM, (void **)&d_mDispersionEx, size);
		cuda_memcpy(GPU, d_mDispersionEx, mDispersionEx, size, cudaMemcpyHostToDevice);

		size = numDispersionEx * sizeof(real_t);
		cuda_malloc(GPU, UM, (void **)&d_DispersionEx, size);
		cuda_memcpy(GPU, d_DispersionEx, DispersionEx, size, cudaMemcpyHostToDevice);
	}
	if (numDispersionEy > 0) {
		size = numDispersionEy * sizeof(dispersion_t);
		cuda_malloc(GPU, UM, (void **)&d_mDispersionEy, size);
		cuda_memcpy(GPU, d_mDispersionEy, mDispersionEy, size, cudaMemcpyHostToDevice);

		size = numDispersionEy * sizeof(real_t);
		cuda_malloc(GPU, UM, (void **)&d_DispersionEy, size);
		cuda_memcpy(GPU, d_DispersionEy, DispersionEy, size, cudaMemcpyHostToDevice);
	}
	if (numDispersionEz > 0) {
		size = numDispersionEz * sizeof(dispersion_t);
		cuda_malloc(GPU, UM, (void **)&d_mDispersionEz, size);
		cuda_memcpy(GPU, d_mDispersionEz, mDispersionEz, size, cudaMemcpyHostToDevice);

		size = numDispersionEz * sizeof(real_t);
		cuda_malloc(GPU, UM, (void **)&d_DispersionEz, size);
		cuda_memcpy(GPU, d_DispersionEz, DispersionEz, size, cudaMemcpyHostToDevice);
	}

	// ABC

	if      (iABC == 0) {
		xsize = numMurHx * sizeof(mur_t);
		ysize = numMurHy * sizeof(mur_t);
		zsize = numMurHz * sizeof(mur_t);
		cuda_malloc(GPU, UM, (void **)&d_fMurHx, xsize);
		cuda_malloc(GPU, UM, (void **)&d_fMurHy, ysize);
		cuda_malloc(GPU, UM, (void **)&d_fMurHz, zsize);
		cuda_memcpy(GPU, d_fMurHx, fMurHx, xsize, cudaMemcpyHostToDevice);
		cuda_memcpy(GPU, d_fMurHy, fMurHy, ysize, cudaMemcpyHostToDevice);
		cuda_memcpy(GPU, d_fMurHz, fMurHz, zsize, cudaMemcpyHostToDevice);
	}
	else if (iABC == 1) {
		xsize = numPmlEx * sizeof(pml_t);
		ysize = numPmlEy * sizeof(pml_t);
		zsize = numPmlEz * sizeof(pml_t);
		cuda_malloc(GPU, UM, (void **)&d_fPmlEx, xsize);
		cuda_malloc(GPU, UM, (void **)&d_fPmlEy, ysize);
		cuda_malloc(GPU, UM, (void **)&d_fPmlEz, zsize);
		cuda_memcpy(GPU, d_fPmlEx, fPmlEx, xsize, cudaMemcpyHostToDevice);
		cuda_memcpy(GPU, d_fPmlEy, fPmlEy, ysize, cudaMemcpyHostToDevice);
		cuda_memcpy(GPU, d_fPmlEz, fPmlEz, zsize, cudaMemcpyHostToDevice);

		xsize = numPmlHx * sizeof(pml_t);
		ysize = numPmlHy * sizeof(pml_t);
		zsize = numPmlHz * sizeof(pml_t);
		cuda_malloc(GPU, UM, (void **)&d_fPmlHx, xsize);
		cuda_malloc(GPU, UM, (void **)&d_fPmlHy, ysize);
		cuda_malloc(GPU, UM, (void **)&d_fPmlHz, zsize);
		cuda_memcpy(GPU, d_fPmlHx, fPmlHx, xsize, cudaMemcpyHostToDevice);
		cuda_memcpy(GPU, d_fPmlHy, fPmlHy, ysize, cudaMemcpyHostToDevice);
		cuda_memcpy(GPU, d_fPmlHz, fPmlHz, zsize, cudaMemcpyHostToDevice);

		xsize = (Nx + (2 * cPML.l)) * sizeof(real_t);
		ysize = (Ny + (2 * cPML.l)) * sizeof(real_t);
		zsize = (Nz + (2 * cPML.l)) * sizeof(real_t);
		cuda_malloc(GPU, UM, (void **)&d_gPmlXn, xsize);
		cuda_malloc(GPU, UM, (void **)&d_gPmlXc, xsize);
		cuda_malloc(GPU, UM, (void **)&d_gPmlYn, ysize);
		cuda_malloc(GPU, UM, (void **)&d_gPmlYc, ysize);
		cuda_malloc(GPU, UM, (void **)&d_gPmlZn, zsize);
		cuda_malloc(GPU, UM, (void **)&d_gPmlZc, zsize);
		cuda_memcpy(GPU, d_gPmlXn, gPmlXn, xsize, cudaMemcpyHostToDevice);
		cuda_memcpy(GPU, d_gPmlXc, gPmlXc, xsize, cudaMemcpyHostToDevice);
		cuda_memcpy(GPU, d_gPmlYn, gPmlYn, ysize, cudaMemcpyHostToDevice);
		cuda_memcpy(GPU, d_gPmlYc, gPmlYc, ysize, cudaMemcpyHostToDevice);
		cuda_memcpy(GPU, d_gPmlZn, gPmlZn, zsize, cudaMemcpyHostToDevice);
		cuda_memcpy(GPU, d_gPmlZc, gPmlZc, zsize, cudaMemcpyHostToDevice);

		size = NMaterial * sizeof(real_t);
		cuda_malloc(GPU, UM, (void **)&d_rPmlE, size);
		cuda_malloc(GPU, UM, (void **)&d_rPmlH, size);
		cuda_malloc(GPU, UM, (void **)&d_rPml,  size);
		cuda_memcpy(GPU, d_rPmlE, rPmlE, size, cudaMemcpyHostToDevice);
		cuda_memcpy(GPU, d_rPmlH, rPmlH, size, cudaMemcpyHostToDevice);
		cuda_memcpy(GPU, d_rPml,  rPml,  size, cudaMemcpyHostToDevice);
	}
}


void setup_host()
{
	// constant memory

	h_Param.Ni = Ni;
	h_Param.Nj = Nj;
	h_Param.Nk = Nk;
	h_Param.N0 = N0;
	h_Param.Nx = Nx;
	h_Param.Ny = Ny;
	h_Param.Nz = Nz;
	h_Param.iMin = iMin;
	h_Param.iMax = iMax;
	h_Param.jMin = jMin;
	h_Param.jMax = jMax;
	h_Param.kMin = kMin;
	h_Param.kMax = kMax;
	h_Param.NFeed = NFeed;
	h_Param.IPlanewave = IPlanewave;
	for (int m = 0; m < 3; m++) {
		h_Param.ei[m] = (real_t)Planewave.ei[m];
		h_Param.hi[m] = (real_t)Planewave.hi[m];
		h_Param.r0[m] = (real_t)Planewave.r0[m];
		h_Param.ri[m] = (real_t)Planewave.ri[m];
	}
	h_Param.ai = (real_t)Planewave.ai;
	h_Param.dt = (real_t)Dt;

	// mesh (real_t)

	h_Xn = (real_t *)malloc((Nx + 1) * sizeof(real_t));
	h_Yn = (real_t *)malloc((Ny + 1) * sizeof(real_t));
	h_Zn = (real_t *)malloc((Nz + 1) * sizeof(real_t));
	for (int i = 0; i <= Nx; i++) {
		h_Xn[i] = (real_t)Xn[i];
	}
	for (int j = 0; j <= Ny; j++) {
		h_Yn[j] = (real_t)Yn[j];
	}
	for (int k = 0; k <= Nz; k++) {
		h_Zn[k] = (real_t)Zn[k];
	}

	h_Xc = (real_t *)malloc((Nx + 0) * sizeof(real_t));
	h_Yc = (real_t *)malloc((Ny + 0) * sizeof(real_t));
	h_Zc = (real_t *)malloc((Nz + 0) * sizeof(real_t));
	for (int i = 0; i < Nx; i++) {
		h_Xc[i] = (real_t)Xc[i];
	}
	for (int j = 0; j < Ny; j++) {
		h_Yc[j] = (real_t)Yc[j];
	}
	for (int k = 0; k < Nz; k++) {
		h_Zc[k] = (real_t)Zc[k];
	}
}
