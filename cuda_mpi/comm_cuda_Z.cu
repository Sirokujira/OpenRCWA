/*
comm_cuda_Z.cu (CUDA + MPI)
share Z plane Hx and Hy
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
static void _d2h(int i, int j, int k, real_t *h, real_t *buf, int i0, int i1, int j0, int j1, param_t *p)
{
	const int64_t m = (i - i0) * (j1 - j0 + 1) + (j - j0);
	const int64_t n = (p->Ni * i) + (p->Nj * j) + (p->Nk * k) + p->N0;
	buf[m] = h[n];
}
__global__
static void d2h_gpu(int k, real_t *h, real_t *buf, int i0, int i1, int j0, int j1)
{
	const int i = i0 + threadIdx.x + (blockIdx.x * blockDim.x);
	const int j = j0 + threadIdx.y + (blockIdx.y * blockDim.y);
	if ((i <= i1) &&
	    (j <= j1)) {
		_d2h(i, j, k, h, buf, i0, i1, j0, j1, &d_Param);
	}
}
static void d2h_cpu(int k, real_t *h, real_t *buf, int i0, int i1, int j0, int j1)
{
	for (int i = i0; i <= i1; i++) {
	for (int j = j0; j <= j1; j++) {
		_d2h(i, j, k, h, buf, i0, i1, j0, j1, &h_Param);
	}
	}
}
static void d2h(int k, real_t h[], real_t buf[], real_t d_buf[], size_t size, int i0, int i1, int j0, int j1)
{
	if (GPU) {
		// parameter
		cudaMemcpyToSymbol(d_Param, &h_Param, sizeof(param_t));

		// grid
		dim3 grid(CEIL(i1 - i0 + 1, bufBlock.x),
		          CEIL(j1 - j0 + 1, bufBlock.y));

		// device
		d2h_gpu<<<grid, bufBlock>>>(k, h, d_buf, i0, i1, j0, j1);

		// device to host
		cuda_memcpy(GPU, buf, d_buf, size, cudaMemcpyDeviceToHost);

		if (UM) cudaDeviceSynchronize();
	}
	else {
		d2h_cpu(k, h, buf, i0, i1, j0, j1);
	}
}


// host to device
__host__ __device__ __inline__
static void _h2d(int i, int j, int k, real_t *h, real_t *buf, int i0, int i1, int j0, int j1, param_t *p)
{
	const int64_t m = (i - i0) * (j1 - j0 + 1) + (j - j0);
	const int64_t n = (p->Ni * i) + (p->Nj * j) + (p->Nk * k) + p->N0;
	h[n] = buf[m];
}
__global__
static void h2d_gpu(int k, real_t *h, real_t *buf, int i0, int i1, int j0, int j1)
{
	const int i = i0 + threadIdx.x + (blockIdx.x * blockDim.x);
	const int j = j0 + threadIdx.y + (blockIdx.y * blockDim.y);
	if ((i <= i1) &&
	    (j <= j1)) {
		_h2d(i, j, k, h, buf, i0, i1, j0, j1, &d_Param);
	}
}
static void h2d_cpu(int k, real_t *h, real_t *buf, int i0, int i1, int j0, int j1)
{
	for (int i = i0; i <= i1; i++) {
	for (int j = j0; j <= j1; j++) {
		_h2d(i, j, k, h, buf, i0, i1, j0, j1, &h_Param);
	}
	}
}
static void h2d(int k, real_t h[], real_t buf[], real_t d_buf[], size_t size, int i0, int i1, int j0, int j1)
{
	if (GPU) {
		// parameter
		cudaMemcpyToSymbol(d_Param, &h_Param, sizeof(param_t));

		// grid
		dim3 grid(CEIL(i1 - i0 + 1, bufBlock.x),
		          CEIL(j1 - j0 + 1, bufBlock.y));

		// host to device
		cuda_memcpy(GPU, d_buf, buf, size, cudaMemcpyHostToDevice);

		// device
		h2d_gpu<<<grid, bufBlock>>>(k, h, d_buf, i0, i1, j0, j1);

		if (UM) cudaDeviceSynchronize();
	}
	else {
		h2d_cpu(k, h, buf, i0, i1, j0, j1);
	}
}


// -Z/+Z boundary
void comm_cuda_Z(int mode)
{
#ifdef _MPI
	MPI_Status status;
	const int tag = 0;
	int bz[][2] = {{Ipz > 0, Ipz < Npz - 1}, {Ipz == 0, Ipz == Npz - 1}};
	int pz[][2] = {{Ipz - 1, Ipz + 1}, {Npz - 1, 0}};
	int krecv[] = {kMin - 1, kMax + 0};
	int ksend[] = {kMin + 0, kMax - 1};
	int ihx[] = {Bid.ihx_z[0], Bid.ihx_z[1]};
	int jhx[] = {Bid.jhx_z[0], Bid.jhx_z[1]};
	int ihy[] = {Bid.ihy_z[0], Bid.ihy_z[1]};
	int jhy[] = {Bid.jhy_z[0], Bid.jhy_z[1]};
	const size_t size_hx = Bid.numhx_z * sizeof(real_t);
	const size_t size_hy = Bid.numhy_z * sizeof(real_t);
	int k;

	for (int side = 0; side < 2; side++) {
		if (bz[mode][side]) {
			// H to buffer
			k = ksend[side];
			d2h(k, Hx, Sendbuf_hx_z, d_Sendbuf_hx_z, size_hx, ihx[0], ihx[1], jhx[0], jhx[1]);
			d2h(k, Hy, Sendbuf_hy_z, d_Sendbuf_hy_z, size_hy, ihy[0], ihy[1], jhy[0], jhy[1]);

			// MPI
			const int ipz = pz[mode][side];
			const int dst = (Ipx * Npy * Npz) + (Ipy * Npz) + ipz;
			MPI_Sendrecv(Sendbuf_hx_z, Bid.numhx_z, MPI_REAL_T, dst, tag,
			             Recvbuf_hx_z, Bid.numhx_z, MPI_REAL_T, dst, tag, MPI_COMM_WORLD, &status);
			MPI_Sendrecv(Sendbuf_hy_z, Bid.numhy_z, MPI_REAL_T, dst, tag,
			             Recvbuf_hy_z, Bid.numhy_z, MPI_REAL_T, dst, tag, MPI_COMM_WORLD, &status);

			// buffer to H
			k = krecv[side];
			h2d(k, Hx, Recvbuf_hx_z, d_Recvbuf_hx_z, size_hx, ihx[0], ihx[1], jhx[0], jhx[1]);
			h2d(k, Hy, Recvbuf_hy_z, d_Recvbuf_hy_z, size_hy, ihy[0], ihy[1], jhy[0], jhy[1]);
		}
	}
#else
	mode = 0;  // dummy
#endif
}
