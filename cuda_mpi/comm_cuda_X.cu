/*
comm_cuda_X.cu (CUDA + MPI)
share X plane Hy and Hz
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
static void _d2h(int i, int j, int k, real_t *h, real_t *buf, int j0, int j1, int k0, int k1, param_t *p)
{
	const int64_t m = (j - j0) * (k1 - k0 + 1) + (k - k0);
	const int64_t n = (p->Ni * i) + (p->Nj * j) + (p->Nk * k) + p->N0;
	buf[m] = h[n];
}
__global__
static void d2h_gpu(int i, real_t *h, real_t *buf, int j0, int j1, int k0, int k1)
{
	const int j = j0 + threadIdx.x + (blockIdx.x * blockDim.x);
	const int k = k0 + threadIdx.y + (blockIdx.y * blockDim.y);
	if ((j <= j1) &&
	    (k <= k1)) {
		_d2h(i, j, k, h, buf, j0, j1, k0, k1, &d_Param);
	}
}
static void d2h_cpu(int i, real_t *h, real_t *buf, int j0, int j1, int k0, int k1)
{
	for (int j = j0; j <= j1; j++) {
	for (int k = k0; k <= k1; k++) {
		_d2h(i, j, k, h, buf, j0, j1, k0, k1, &h_Param);
	}
	}
}
static void d2h(int i, real_t h[], real_t buf[], real_t d_buf[], size_t size, int j0, int j1, int k0, int k1)
{
	if (GPU) {
		// parameter
		cudaMemcpyToSymbol(d_Param, &h_Param, sizeof(param_t));

		// grid
		dim3 grid(CEIL(j1 - j0 + 1, bufBlock.x),
		          CEIL(k1 - k0 + 1, bufBlock.y));

		// device
		d2h_gpu<<<grid, bufBlock>>>(i, h, d_buf, j0, j1, k0, k1);

		// device to host
		cuda_memcpy(GPU, buf, d_buf, size, cudaMemcpyDeviceToHost);

		if (UM) cudaDeviceSynchronize();
	}
	else {
		d2h_cpu(i, h, buf, j0, j1, k0, k1);
	}
}


// host to device
__host__ __device__ __inline__
static void _h2d(int i, int j, int k, real_t *h, real_t *buf, int j0, int j1, int k0, int k1, param_t *p)
{
	const int64_t m = (j - j0) * (k1 - k0 + 1) + (k - k0);
	const int64_t n = (p->Ni * i) + (p->Nj * j) + (p->Nk * k) + p->N0;
	h[n] = buf[m];
}
__global__
static void h2d_gpu(int i, real_t *h, real_t *buf, int j0, int j1, int k0, int k1)
{
	const int j = j0 + threadIdx.x + (blockIdx.x * blockDim.x);
	const int k = k0 + threadIdx.y + (blockIdx.y * blockDim.y);
	if ((j <= j1) &&
	    (k <= k1)) {
		_h2d(i, j, k, h, buf, j0, j1, k0, k1, &d_Param);
	}
}
static void h2d_cpu(int i, real_t *h, real_t *buf, int j0, int j1, int k0, int k1)
{
	for (int j = j0; j <= j1; j++) {
	for (int k = k0; k <= k1; k++) {
		_h2d(i, j, k, h, buf, j0, j1, k0, k1, &h_Param);
	}
	}
}
static void h2d(int i, real_t h[], real_t buf[], real_t d_buf[], size_t size, int j0, int j1, int k0, int k1)
{
	if (GPU) {
		// parameter
		cudaMemcpyToSymbol(d_Param, &h_Param, sizeof(param_t));

		// grid
		dim3 grid(CEIL(j1 - j0 + 1, bufBlock.x),
		          CEIL(k1 - k0 + 1, bufBlock.y));

		// host to device
		cuda_memcpy(GPU, d_buf, buf, size, cudaMemcpyHostToDevice);

		// device
		h2d_gpu<<<grid, bufBlock>>>(i, h, d_buf, j0, j1, k0, k1);

		if (UM) cudaDeviceSynchronize();
	}
	else {
		h2d_cpu(i, h, buf, j0, j1, k0, k1);
	}
}


// -X/+X boundary
void comm_cuda_X(int mode)
{
#ifdef _MPI
	MPI_Status status;
	const int tag = 0;
	int bx[][2] = {{Ipx > 0, Ipx < Npx - 1}, {Ipx == 0, Ipx == Npx - 1}};
	int px[][2] = {{Ipx - 1, Ipx + 1}, {Npx - 1, 0}};
	int irecv[] = {iMin - 1, iMax + 0};
	int isend[] = {iMin + 0, iMax - 1};
	int jhy[] = {Bid.jhy_x[0], Bid.jhy_x[1]};
	int khy[] = {Bid.khy_x[0], Bid.khy_x[1]};
	int jhz[] = {Bid.jhz_x[0], Bid.jhz_x[1]};
	int khz[] = {Bid.khz_x[0], Bid.khz_x[1]};
	const size_t size_hy = Bid.numhy_x * sizeof(real_t);
	const size_t size_hz = Bid.numhz_x * sizeof(real_t);
	int i;

	for (int side = 0; side < 2; side++) {
		if (bx[mode][side]) {
			// H to buffer
			i = isend[side];
			d2h(i, Hy, Sendbuf_hy_x, d_Sendbuf_hy_x, size_hy, jhy[0], jhy[1], khy[0], khy[1]);
			d2h(i, Hz, Sendbuf_hz_x, d_Sendbuf_hz_x, size_hz, jhz[0], jhz[1], khz[0], khz[1]);

			// MPI
			const int ipx = px[mode][side];
			const int dst = (ipx * Npy * Npz) + (Ipy * Npz) + Ipz;
			MPI_Sendrecv(Sendbuf_hy_x, Bid.numhy_x, MPI_REAL_T, dst, tag,
			             Recvbuf_hy_x, Bid.numhy_x, MPI_REAL_T, dst, tag, MPI_COMM_WORLD, &status);
			MPI_Sendrecv(Sendbuf_hz_x, Bid.numhz_x, MPI_REAL_T, dst, tag,
			             Recvbuf_hz_x, Bid.numhz_x, MPI_REAL_T, dst, tag, MPI_COMM_WORLD, &status);

			// buffer to H
			i = irecv[side];
			h2d(i, Hy, Recvbuf_hy_x, d_Recvbuf_hy_x, size_hy, jhy[0], jhy[1], khy[0], khy[1]);
			h2d(i, Hz, Recvbuf_hz_x, d_Recvbuf_hz_x, size_hz, jhz[0], jhz[1], khz[0], khz[1]);
		}
	}
#else
	mode = 0;  // dummy
#endif
}
