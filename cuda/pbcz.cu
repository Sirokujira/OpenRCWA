/*
pbcz.cu (CUDA)

PBC on +/- Z boundary
*/

#include "orcwa.h"
#include "orcwa_cuda.h"


__host__ __device__
static void _pbczhx(int i, int j, real_t *hx, param_t *p)
{
	hx[LA(p, i, j,    -1)] = hx[LA(p, i, j, p->Nz - 1)];
	hx[LA(p, i, j, p->Nz)] = hx[LA(p, i, j,         0)];
}


__host__ __device__
static void _pbczhy(int i, int j, real_t *hy, param_t *p)
{
	hy[LA(p, i, j,    -1)] = hy[LA(p, i, j, p->Nz - 1)];
	hy[LA(p, i, j, p->Nz)] = hy[LA(p, i, j,         0)];
}


__global__
static void pbczhx_gpu(real_t *hx)
{
	int i = d_Param.iMin - 0 + (blockDim.x * blockIdx.x) + threadIdx.x;
	int j = d_Param.jMin - 1 + (blockDim.y * blockIdx.y) + threadIdx.y;

	if (i <= d_Param.iMax) {
	if (j <= d_Param.jMax) {
		_pbczhx(i, j, hx, &d_Param);
	}
	}
}


__global__
static void pbczhy_gpu(real_t *hy)
{
	const int i = d_Param.iMin - 1 + (blockDim.x * blockIdx.x) + threadIdx.x;
	const int j = d_Param.jMin - 0 + (blockDim.y * blockIdx.y) + threadIdx.y;

	if (i <= d_Param.iMax) {
	if (j <= d_Param.jMax) {
		_pbczhy(i, j, hy, &d_Param);
	}
	}
}


static void pbczhx_cpu(real_t *hx)
{
	for (int i = h_Param.iMin - 0; i <= h_Param.iMax; i++) {
	for (int j = h_Param.jMin - 1; j <= h_Param.jMax; j++) {
		_pbczhx(i, j, hx, &h_Param);
	}
	}
}


static void pbczhy_cpu(real_t *hy)
{
	for (int i = h_Param.iMin - 1; i <= h_Param.iMax; i++) {
	for (int j = h_Param.jMin - 0; j <= h_Param.jMax; j++) {
		_pbczhy(i, j, hy, &h_Param);
	}
	}
}


void pbcz()
{
	if (GPU) {
		cudaMemcpyToSymbol(d_Param, &h_Param, sizeof(param_t));
		dim3 block(pbcBlock, pbcBlock);
		dim3 grid_hx(CEIL(iMax - iMin + 1, block.x),
		             CEIL(jMax - jMin + 2, block.y));
		dim3 grid_hy(CEIL(iMax - iMin + 2, block.x),
		             CEIL(jMax - jMin + 1, block.y));
		pbczhx_gpu<<<grid_hx, block>>>(Hx);
		pbczhy_gpu<<<grid_hy, block>>>(Hy);
		if (UM) cudaDeviceSynchronize();
	}
	else {
		pbczhx_cpu(Hx);
		pbczhy_cpu(Hy);
	}
}
