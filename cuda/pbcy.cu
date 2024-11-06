/*
pbcy.cu (CUDA)

PBC on +/- Y boundary
*/

#include "orcwa.h"
#include "orcwa_cuda.h"


__host__ __device__
static void _pbcyhz(int k, int i, real_t *hz, param_t *p)
{
	hz[LA(p, i,    -1, k)] = hz[LA(p, i, p->Ny - 1, k)];
	hz[LA(p, i, p->Ny, k)] = hz[LA(p, i,         0, k)];
}


__host__ __device__
static void _pbcyhx(int k, int i, real_t *hx, param_t *p)
{
	hx[LA(p, i,    -1, k)] = hx[LA(p, i, p->Ny - 1, k)];
	hx[LA(p, i, p->Ny, k)] = hx[LA(p, i,         0, k)];
}


__global__
static void pbcyhz_gpu(real_t *hz)
{
	int k = d_Param.kMin - 0 + (blockDim.x * blockIdx.x) + threadIdx.x;
	int i = d_Param.iMin - 1 + (blockDim.y * blockIdx.y) + threadIdx.y;

	if (k <= d_Param.kMax) {
	if (i <= d_Param.iMax) {
		_pbcyhz(k, i, hz, &d_Param);
	}
	}
}


__global__
static void pbcyhx_gpu(real_t *hx)
{
	const int k = d_Param.kMin - 1 + (blockDim.x * blockIdx.x) + threadIdx.x;
	const int i = d_Param.iMin - 0 + (blockDim.y * blockIdx.y) + threadIdx.y;

	if (k <= d_Param.kMax) {
	if (i <= d_Param.iMax) {
		_pbcyhx(k, i, hx, &d_Param);
	}
	}
}


static void pbcyhz_cpu(real_t *hz)
{
	for (int k = h_Param.kMin - 0; k <= h_Param.kMax; k++) {
	for (int i = h_Param.iMin - 1; i <= h_Param.iMax; i++) {
		_pbcyhz(k, i, hz, &h_Param);
	}
	}
}


static void pbcyhx_cpu(real_t *hx)
{
	for (int k = h_Param.kMin - 1; k <= h_Param.kMax; k++) {
	for (int i = h_Param.iMin - 0; i <= h_Param.iMax; i++) {
		_pbcyhx(k, i, hx, &h_Param);
	}
	}
}


void pbcy()
{
	if (GPU) {
		cudaMemcpyToSymbol(d_Param, &h_Param, sizeof(param_t));
		dim3 block(pbcBlock, pbcBlock);
		dim3 grid_hz(CEIL(kMax - kMin + 1, block.x),
		             CEIL(iMax - iMin + 2, block.y));
		dim3 grid_hx(CEIL(kMax - kMin + 2, block.x),
		             CEIL(iMax - iMin + 1, block.y));
		pbcyhz_gpu<<<grid_hz, block>>>(Hz);
		pbcyhx_gpu<<<grid_hx, block>>>(Hx);
		if (UM) cudaDeviceSynchronize();
	}
	else {
		pbcyhz_cpu(Hz);
		pbcyhx_cpu(Hx);
	}
}
