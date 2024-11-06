/*
memallocfree3_gpu.cu (CUDA)

alloc and free
(3) near3d
*/

#include "orcwa.h"
#include "orcwa_cuda.h"
#include "orcwa_prototype.h"

void memalloc3_gpu()
{
	if ((NN > 0) && (NFreq2 > 0)) {
		const size_t size = NN * NFreq2 * sizeof(float);

		cuda_malloc(GPU, UM, (void **)&d_cEx_r, size);
		cuda_malloc(GPU, UM, (void **)&d_cEy_r, size);
		cuda_malloc(GPU, UM, (void **)&d_cEz_r, size);
		cuda_malloc(GPU, UM, (void **)&d_cHx_r, size);
		cuda_malloc(GPU, UM, (void **)&d_cHy_r, size);
		cuda_malloc(GPU, UM, (void **)&d_cHz_r, size);
		cuda_malloc(GPU, UM, (void **)&d_cEx_i, size);
		cuda_malloc(GPU, UM, (void **)&d_cEy_i, size);
		cuda_malloc(GPU, UM, (void **)&d_cEz_i, size);
		cuda_malloc(GPU, UM, (void **)&d_cHx_i, size);
		cuda_malloc(GPU, UM, (void **)&d_cHy_i, size);
		cuda_malloc(GPU, UM, (void **)&d_cHz_i, size);
	}
}


void memcopy3_gpu()
{
	if ((NN > 0) && (NFreq2 > 0)) {
		const size_t size = NN * NFreq2 * sizeof(float);

		cEx_r = (float *)malloc(size);
		cEx_i = (float *)malloc(size);
		cEy_r = (float *)malloc(size);
		cEy_i = (float *)malloc(size);
		cEz_r = (float *)malloc(size);
		cEz_i = (float *)malloc(size);
		cHx_r = (float *)malloc(size);
		cHx_i = (float *)malloc(size);
		cHy_r = (float *)malloc(size);
		cHy_i = (float *)malloc(size);
		cHz_r = (float *)malloc(size);
		cHz_i = (float *)malloc(size);

		cuda_memcpy(GPU, cEx_r, d_cEx_r, size, cudaMemcpyDeviceToHost);
		cuda_memcpy(GPU, cEy_r, d_cEy_r, size, cudaMemcpyDeviceToHost);
		cuda_memcpy(GPU, cEz_r, d_cEz_r, size, cudaMemcpyDeviceToHost);
		cuda_memcpy(GPU, cHx_r, d_cHx_r, size, cudaMemcpyDeviceToHost);
		cuda_memcpy(GPU, cHy_r, d_cHy_r, size, cudaMemcpyDeviceToHost);
		cuda_memcpy(GPU, cHz_r, d_cHz_r, size, cudaMemcpyDeviceToHost);
		cuda_memcpy(GPU, cEx_i, d_cEx_i, size, cudaMemcpyDeviceToHost);
		cuda_memcpy(GPU, cEy_i, d_cEy_i, size, cudaMemcpyDeviceToHost);
		cuda_memcpy(GPU, cEz_i, d_cEz_i, size, cudaMemcpyDeviceToHost);
		cuda_memcpy(GPU, cHx_i, d_cHx_i, size, cudaMemcpyDeviceToHost);
		cuda_memcpy(GPU, cHy_i, d_cHy_i, size, cudaMemcpyDeviceToHost);
		cuda_memcpy(GPU, cHz_i, d_cHz_i, size, cudaMemcpyDeviceToHost);
	}
}


void memfree3_gpu()
{
	if ((NN > 0) && (NFreq2 > 0)) {
		cuda_free(GPU, d_cEx_r);
		cuda_free(GPU, d_cEy_r);
		cuda_free(GPU, d_cEz_r);
		cuda_free(GPU, d_cHx_r);
		cuda_free(GPU, d_cHy_r);
		cuda_free(GPU, d_cHz_r);
		cuda_free(GPU, d_cEx_i);
		cuda_free(GPU, d_cEy_i);
		cuda_free(GPU, d_cEz_i);
		cuda_free(GPU, d_cHx_i);
		cuda_free(GPU, d_cHy_i);
		cuda_free(GPU, d_cHz_i);
	}
}

