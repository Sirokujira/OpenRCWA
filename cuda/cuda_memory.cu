/*
cuda_memory.cu

CUDA memory utilities
*/

#include <stdlib.h>

// malloc and clear
void cuda_malloc(int gpu, int um, void **ptr, size_t size)
{
	if (gpu) {
		if (um) {
			cudaMallocManaged(ptr, size);
		}
		else {
			cudaMalloc(ptr, size);
		}
		cudaMemset(*ptr, 0, size);
	}
	else {
		*ptr = malloc(size);
		memset(*ptr, 0, size);
	}
}

// free
void cuda_free(int gpu, void *ptr)
{
	if (gpu) {
		cudaFree(ptr);
	}
	else {
		free(ptr);
	}
}

// memset
void cuda_memset(int gpu, void *ptr, int c, size_t size)
{
	if (gpu) {
		cudaMemset(ptr, c, size);
	}
	else {
		memset(ptr, c, size);
	}
}

// memcpy
void cuda_memcpy(int gpu, void *dst, const void *src, size_t size, cudaMemcpyKind kind)
{
	if (gpu) {
		cudaMemcpy(dst, src, size, kind);
	}
	else {
		memcpy(dst, src, size);
	}
}
