/*
average.cu (CUDA)
*/

#include "orcwa.h"
#include "orcwa_cuda.h"
#include "reduction_sum.cu"

__host__ __device__
static real_t average_e(real_t *ex, real_t *ey, real_t *ez, int i, int j, int k, param_t *p)
{
	return
		fabs(
			+ ex[LA(p, i,     j,     k    )]
			+ ex[LA(p, i,     j + 1, k    )]
			+ ex[LA(p, i,     j,     k + 1)]
			+ ex[LA(p, i,     j + 1, k + 1)]) +
		fabs(
			+ ey[LA(p, i,     j,     k    )]
			+ ey[LA(p, i,     j,     k + 1)]
			+ ey[LA(p, i + 1, j,     k    )]
			+ ey[LA(p, i + 1, j,     k + 1)]) +
		fabs(
			+ ez[LA(p, i,     j,     k    )]
			+ ez[LA(p, i + 1, j,     k    )]
			+ ez[LA(p, i,     j + 1, k    )]
			+ ez[LA(p, i + 1, j + 1, k    )]);
}

__host__ __device__
static real_t average_h(real_t *hx, real_t *hy, real_t *hz, int i, int j, int k, param_t *p)
{
	return
		fabs(
			+ hx[LA(p, i,     j,     k    )]
			+ hx[LA(p, i + 1, j,     k    )]) +
		fabs(
			+ hy[LA(p, i,     j,     k    )]
			+ hy[LA(p, i,     j + 1, k    )]) +
		fabs(
			+ hz[LA(p, i,     j,     k    )]
			+ hz[LA(p, i,     j,     k + 1)]);
}

// GPU
__global__
static void average_gpu(
	int imin, int imax, int jmin, int jmax, int kmin, int kmax,
	real_t *ex, real_t *ey, real_t *ez, real_t *hx, real_t *hy, real_t *hz,
	real_t *sume, real_t *sumh)
{
	extern __shared__ real_t sm[];

	const int i = imin + threadIdx.z + (blockIdx.z * blockDim.z);
	const int j = jmin + threadIdx.y + (blockIdx.y * blockDim.y);
	const int k = kmin + threadIdx.x + (blockIdx.x * blockDim.x);

	const int nthread = blockDim.x * blockDim.y * blockDim.z;
	const int tid = threadIdx.x + (threadIdx.y * blockDim.x) + (threadIdx.z * blockDim.x * blockDim.y);
	const int bid = blockIdx.x + (blockIdx.y * gridDim.x) + (blockIdx.z * gridDim.x * gridDim.y);

	if ((i < imax) && (j < jmax) && (k < kmax)) {
		sm[tid] = average_e(ex, ey, ez, i, j, k, &d_Param);
	}
	else {
		sm[tid] = 0;
	}
	reduction_sum(tid, nthread, sm, &sume[bid]);

	if ((i < imax) && (j < jmax) && (k < kmax)) {
		sm[tid] = average_h(hx, hy, hz, i, j, k, &d_Param);
	}
	else {
		sm[tid] = 0;
	}
	reduction_sum(tid, nthread, sm, &sumh[bid]);
}

// CPU
static void average_cpu(real_t *sum)
{
	sum[0] = 0;
	sum[1] = 0;
	for (int i = iMin; i < iMax; i++) {
	for (int j = jMin; j < jMax; j++) {
	for (int k = kMin; k < kMax; k++) {
		sum[0] += average_e(Ex, Ey, Ez, i, j, k, &h_Param);
		sum[1] += average_h(Hx, Hy, Hz, i, j, k, &h_Param);
	}
	}
	}
}

void average(double fsum[])
{
	real_t sum[2];

	// sum
	if (GPU) {
		cudaMemcpyToSymbol(d_Param, &h_Param, sizeof(param_t));

		const int sm_size = sumBlock.x * sumBlock.y * sumBlock.z * sizeof(real_t);
		average_gpu<<<sumGrid, sumBlock, sm_size>>>(
			iMin, iMax, jMin, jMax, kMin, kMax, Ex, Ey, Ez, Hx, Hy, Hz, d_sumE, d_sumH);

		// device to host
		const int ngrid = sumGrid.x * sumGrid.y * sumGrid.z;
		const size_t size = ngrid * sizeof(real_t);
		cudaMemcpy(h_sumE, d_sumE, size, cudaMemcpyDeviceToHost);
		cudaMemcpy(h_sumH, d_sumH, size, cudaMemcpyDeviceToHost);

		// sum of partial sum
		sum[0] = 0;
		sum[1] = 0;
		for (int n = 0; n < ngrid; n++) {
			sum[0] += h_sumE[n];
			sum[1] += h_sumH[n];
		}
	}
	else {
		average_cpu(sum);
	}

	// average
	fsum[0] = sum[0] / (4.0 * Nx * Ny * Nz);
	fsum[1] = sum[1] / (2.0 * Nx * Ny * Nz);
}
