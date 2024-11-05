/*
check_gpu.cu

return = 0/1 : OK/NG
device : I : device number (=0,1,...)
msg    : O : GPU properties or error message
*/

#include <stdio.h>
#include <string.h>

int check_gpu(int idevice, char msg[])
{
	cudaError_t ierr;

	// check CUDA support
	int deviceCount;
	cudaGetDeviceCount(&deviceCount);
	if (deviceCount <= 0) {
		strcpy(msg, "*** There is no device supporting CUDA");
		return 1;
	}
	//printf("deviceCount=%d\n", deviceCount);

	// check device number
	if ((idevice < 0) || (idevice >= deviceCount)) {
		sprintf(msg, "*** Invalid device number = %d", idevice);
		return 1;
	}

	// set device
	ierr = cudaSetDevice(idevice);
	if (ierr != cudaSuccess) {
		strcpy(msg, cudaGetErrorString(ierr));
		return 1;
	}

	// properties
	cudaDeviceProp prop;
	ierr = cudaGetDeviceProperties(&prop, idevice);
	if (ierr != cudaSuccess) {
		sprintf(msg, cudaGetErrorString(ierr));
		return 1;
	}
	if (prop.major < 3) {
		sprintf(msg, "*** Compute Capability < 3.0");
		return 1;
	}

	// GPU info
	sprintf(msg, "%s, %dMB, %dMP, C.C.%d.%d",
		prop.name,
		(int)(prop.totalGlobalMem / 1024 / 1024),
		prop.multiProcessorCount,
		prop.major,
		prop.minor);

	return 0;
}
