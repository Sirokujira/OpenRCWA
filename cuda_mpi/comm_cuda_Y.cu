/*
comm_cuda_Y.cu (CUDA + MPI)
share Y plane Hz and Hx
mode = 0/1 : boundary/PBC
*/

#ifdef _MPI
#include <mpi.h>
#endif

#include "orcwa.h"
#include "orcwa_cuda.h"
#include "orcwa_prototype.h"

// device to host
__host__ __device__ __inline__
static void _d2h(int i, int j, int k, real_t *h, real_t *buf, int k0, int k1, int i0, int i1, param_t *p)
{
	const int64_t m = (k - k0) * (i1 - i0 + 1) + (i - i0);
	const int64_t n = (p->Ni * i) + (p->Nj * j) + (p->Nk * k) + p->N0;
	buf[m] = h[n];
}
__global__
static void d2h_gpu(int j, real_t *h, real_t *buf, int k0, int k1, int i0, int i1)
{
	const int k = k0 + threadIdx.x + (blockIdx.x * blockDim.x);
	const int i = i0 + threadIdx.y + (blockIdx.y * blockDim.y);
	if ((k <= k1) &&
	    (i <= i1)) {
		_d2h(i, j, k, h, buf, k0, k1, i0, i1, &d_Param);
	}
}
static void d2h_cpu(int j, real_t *h, real_t *buf, int k0, int k1, int i0, int i1)
{
	for (int k = k0; k <= k1; k++) {
	for (int i = i0; i <= i1; i++) {
		_d2h(i, j, k, h, buf, k0, k1, i0, i1, &h_Param);
	}
	}
}
static void d2h(int j, real_t h[], real_t buf[], real_t d_buf[], size_t size, int k0, int k1, int i0, int i1)
{
	if (GPU) {
		// parameter
		cudaMemcpyToSymbol(d_Param, &h_Param, sizeof(param_t));

		// grid
		dim3 grid(CEIL(k1 - k0 + 1, bufBlock.x),
		          CEIL(i1 - i0 + 1, bufBlock.y));

		// device
		d2h_gpu<<<grid, bufBlock>>>(j, h, d_buf, k0, k1, i0, i1);

		// device to host
		cuda_memcpy(GPU, buf, d_buf, size, cudaMemcpyDeviceToHost);

		if (UM) cudaDeviceSynchronize();
	}
	else {
		d2h_cpu(j, h, buf, k0, k1, i0, i1);
	}
}


// host to device
__host__ __device__ __inline__
static void _h2d(int i, int j, int k, real_t *h, real_t *buf, int k0, int k1, int i0, int i1, param_t *p)
{
	const int64_t m = (k - k0) * (i1 - i0 + 1) + (i - i0);
	const int64_t n = (p->Ni * i) + (p->Nj * j) + (p->Nk * k) + p->N0;
	h[n] = buf[m];
}
__global__
static void h2d_gpu(int j, real_t *h, real_t *buf, int k0, int k1, int i0, int i1)
{
	const int k = k0 + threadIdx.x + (blockIdx.x * blockDim.x);
	const int i = i0 + threadIdx.y + (blockIdx.y * blockDim.y);
	if ((k <= k1) &&
	    (i <= i1)) {
		_h2d(i, j, k, h, buf, k0, k1, i0, i1, &d_Param);
	}
}
static void h2d_cpu(int j, real_t *h, real_t *buf, int k0, int k1, int i0, int i1)
{
	for (int k = k0; k <= k1; k++) {
	for (int i = i0; i <= i1; i++) {
		_h2d(i, j, k, h, buf, k0, k1, i0, i1, &h_Param);
	}
	}
}
static void h2d(int j, real_t h[], real_t buf[], real_t d_buf[], size_t size, int k0, int k1, int i0, int i1)
{
	if (GPU) {
		// parameter
		cudaMemcpyToSymbol(d_Param, &h_Param, sizeof(param_t));

		// grid
		dim3 grid(CEIL(k1 - k0 + 1, bufBlock.x),
		          CEIL(i1 - i0 + 1, bufBlock.y));

		// host to device
		cuda_memcpy(GPU, d_buf, buf, size, cudaMemcpyHostToDevice);

		// device
		h2d_gpu<<<grid, bufBlock>>>(j, h, d_buf, k0, k1, i0, i1);

		if (UM) cudaDeviceSynchronize();
	}
	else {
		h2d_cpu(j, h, buf, k0, k1, i0, i1);
	}
}


// -Y/+Y boundary
void comm_cuda_Y(int mode)
{
#ifdef _MPI
	MPI_Status status;
	const int tag = 0;
	int by[][2] = {{Ipy > 0, Ipy < Npy - 1}, {Ipy == 0, Ipy == Npy - 1}};
	int py[][2] = {{Ipy - 1, Ipy + 1}, {Npy - 1, 0}};
	int jrecv[] = {jMin - 1, jMax + 0};
	int jsend[] = {jMin + 0, jMax - 1};
	int khz[] = {Bid.khz_y[0], Bid.khz_y[1]};
	int ihz[] = {Bid.ihz_y[0], Bid.ihz_y[1]};
	int khx[] = {Bid.khx_y[0], Bid.khx_y[1]};
	int ihx[] = {Bid.ihx_y[0], Bid.ihx_y[1]};
	const size_t size_hz = Bid.numhz_y * sizeof(real_t);
	const size_t size_hx = Bid.numhx_y * sizeof(real_t);
	int j;

	for (int side = 0; side < 2; side++) {
		if (by[mode][side]) {
			// H to buffer
			j = jsend[side];
			d2h(j, Hz, Sendbuf_hz_y, d_Sendbuf_hz_y, size_hz, khz[0], khz[1], ihz[0], ihz[1]);
			d2h(j, Hx, Sendbuf_hx_y, d_Sendbuf_hx_y, size_hx, khx[0], khx[1], ihx[0], ihx[1]);

			// MPI
			const int ipy = py[mode][side];
			const int dst = (Ipx * Npy * Npz) + (ipy * Npz) + Ipz;
			MPI_Sendrecv(Sendbuf_hz_y, Bid.numhz_y, MPI_REAL_T, dst, tag,
			             Recvbuf_hz_y, Bid.numhz_y, MPI_REAL_T, dst, tag, MPI_COMM_WORLD, &status);
			MPI_Sendrecv(Sendbuf_hx_y, Bid.numhx_y, MPI_REAL_T, dst, tag,
			             Recvbuf_hx_y, Bid.numhx_y, MPI_REAL_T, dst, tag, MPI_COMM_WORLD, &status);

			// buffer to H
			j = jrecv[side];
			h2d(j, Hz, Recvbuf_hz_y, d_Recvbuf_hz_y, size_hz, khz[0], khz[1], ihz[0], ihz[1]);
			h2d(j, Hx, Recvbuf_hx_y, d_Recvbuf_hx_y, size_hx, khx[0], khx[1], ihx[0], ihx[1]);
		}
	}
#else
	mode = 0;  // dummy
#endif
}
