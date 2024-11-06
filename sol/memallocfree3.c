/*
memallocfree3.c

alloc and free
(3) near3d
*/

#include "orcwa.h"
#include "orcwa_prototype.h"

void memalloc3(void)
{
	if ((NN > 0) && (NFreq2 > 0)) {
		const size_t size = NFreq2 * NN * sizeof(float);

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
	}
}


void memfree3(void)
{
	if ((NN > 0) && (NFreq2 > 0)) {
		free(cEx_r);
		free(cEx_i);
		free(cEy_r);
		free(cEy_i);
		free(cEz_r);
		free(cEz_i);
		free(cHx_r);
		free(cHx_i);
		free(cHy_r);
		free(cHy_i);
		free(cHz_r);
		free(cHz_i);
	}
}
