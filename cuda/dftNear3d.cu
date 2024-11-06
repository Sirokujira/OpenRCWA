/*
dftNear3d.cu (CUDA)
*/

#include "orcwa.h"
#include "orcwa_cuda.h"

__host__ __device__ __forceinline__
static void dftadd(float *f_r, float *f_i, real_t f, float fctr_r, float fctr_i)
{
	*f_r += (float)f * fctr_r;
	*f_i += (float)f * fctr_i;
}


// Ex
__global__
static void dft_near3dEx_gpu(
	real_t *ex, float *cex_r, float *cex_i, float f_r, float f_i,
	int imin, int imax, int jmin, int jmax, int kmin, int kmax,
	int64_t adr0, int64_t ni, int64_t nj, int64_t nk, int64_t n0)
{
	const int i = imin + threadIdx.z + (blockIdx.z * blockDim.z);
	const int j = jmin + threadIdx.y + (blockIdx.y * blockDim.y);
	const int k = kmin + threadIdx.x + (blockIdx.x * blockDim.x);
	if (i < imax + 0) {
	if (j < jmax + 1) {
	if (k < kmax + 1) {
		const int64_t m = (ni * i) + (nj * j) + (nk * k) + n0;
		const int64_t n = adr0 + m;
		dftadd(&cex_r[n], &cex_i[n], ex[m], f_r, f_i);
	}
	}
	}
}


// Ey
__global__
static void dft_near3dEy_gpu(
	real_t *ey, float *cey_r, float *cey_i, float f_r, float f_i,
	int imin, int imax, int jmin, int jmax, int kmin, int kmax,
	int64_t adr0, int64_t ni, int64_t nj, int64_t nk, int64_t n0)
{
	const int i = imin + threadIdx.z + (blockIdx.z * blockDim.z);
	const int j = jmin + threadIdx.y + (blockIdx.y * blockDim.y);
	const int k = kmin + threadIdx.x + (blockIdx.x * blockDim.x);
	if (i < imax + 1) {
	if (j < jmax + 0) {
	if (k < kmax + 1) {
		const int64_t m = (ni * i) + (nj * j) + (nk * k) + n0;
		const int64_t n = adr0 + m;
		dftadd(&cey_r[n], &cey_i[n], ey[m], f_r, f_i);
	}
	}
	}
}


// Ez
__global__
static void dft_near3dEz_gpu(
	real_t *ez, float *cez_r, float *cez_i, float f_r, float f_i,
	int imin, int imax, int jmin, int jmax, int kmin, int kmax,
	int64_t adr0, int64_t ni, int64_t nj, int64_t nk, int64_t n0)
{
	const int i = imin + threadIdx.z + (blockIdx.z * blockDim.z);
	const int j = jmin + threadIdx.y + (blockIdx.y * blockDim.y);
	const int k = kmin + threadIdx.x + (blockIdx.x * blockDim.x);
	if (i < imax + 1) {
	if (j < jmax + 1) {
	if (k < kmax + 0) {
		const int64_t m = (ni * i) + (nj * j) + (nk * k) + n0;
		const int64_t n = adr0 + m;
		dftadd(&cez_r[n], &cez_i[n], ez[m], f_r, f_i);
	}
	}
	}
}


// Hx
__global__
static void dft_near3dHx_gpu(
	real_t *hx, float *chx_r, float *chx_i, float f_r, float f_i,
	int imin, int imax, int jmin, int jmax, int kmin, int kmax,
	int64_t adr0, int64_t ni, int64_t nj, int64_t nk, int64_t n0)
{
	const int i = imin - 0 + threadIdx.z + (blockIdx.z * blockDim.z);
	const int j = jmin - 1 + threadIdx.y + (blockIdx.y * blockDim.y);
	const int k = kmin - 1 + threadIdx.x + (blockIdx.x * blockDim.x);
	if (i < imax + 1) {
	if (j < jmax + 1) {
	if (k < kmax + 1) {
		const int64_t m = (ni * i) + (nj * j) + (nk * k) + n0;
		const int64_t n = adr0 + m;
		dftadd(&chx_r[n], &chx_i[n], hx[m], f_r, f_i);
	}
	}
	}
}


// Hy
__global__
static void dft_near3dHy_gpu(
	real_t *hy, float *chy_r, float *chy_i, float f_r, float f_i,
	int imin, int imax, int jmin, int jmax, int kmin, int kmax,
	int64_t adr0, int64_t ni, int64_t nj, int64_t nk, int64_t n0)
{
	const int i = imin - 1 + threadIdx.z + (blockIdx.z * blockDim.z);
	const int j = jmin - 0 + threadIdx.y + (blockIdx.y * blockDim.y);
	const int k = kmin - 1 + threadIdx.x + (blockIdx.x * blockDim.x);
	if (i < imax + 1) {
	if (j < jmax + 1) {
	if (k < kmax + 1) {
		const int64_t m = (ni * i) + (nj * j) + (nk * k) + n0;
		const int64_t n = adr0 + m;
		dftadd(&chy_r[n], &chy_i[n], hy[m], f_r, f_i);
	}
	}
	}
}


// Hz
__global__
static void dft_near3dHz_gpu(
	real_t *hz, float *chz_r, float *chz_i, float f_r, float f_i,
	int imin, int imax, int jmin, int jmax, int kmin, int kmax,
	int64_t adr0, int64_t ni, int64_t nj, int64_t nk, int64_t n0)
{
	const int i = imin - 1 + threadIdx.z + (blockIdx.z * blockDim.z);
	const int j = jmin - 1 + threadIdx.y + (blockIdx.y * blockDim.y);
	const int k = kmin - 0 + threadIdx.x + (blockIdx.x * blockDim.x);
	if (i < imax + 1) {
	if (j < jmax + 1) {
	if (k < kmax + 1) {
		const int64_t m = (ni * i) + (nj * j) + (nk * k) + n0;
		const int64_t n = adr0 + m;
		dftadd(&chz_r[n], &chz_i[n], hz[m], f_r, f_i);
	}
	}
	}
}

// CPU

// Ex
static void dft_near3dEx_cpu(
	real_t *ex, float *cex_r, float *cex_i, float f_r, float f_i,
	int imin, int imax, int jmin, int jmax, int kmin, int kmax,
	int64_t adr0, param_t *p)
{
	for (int i = imin; i < imax + 0; i++) {
	for (int j = jmin; j < jmax + 1; j++) {
	for (int k = kmin; k < kmax + 1; k++) {
		const int64_t m = LA(p, i, j, k);
		assert((m >= 0) && (m < NN));
		const int64_t n = adr0 + m;
		dftadd(&cex_r[n], &cex_i[n], ex[m], f_r, f_i);
	}
	}
	}
}


// Ey
static void dft_near3dEy_cpu(
	real_t *ey, float *cey_r, float *cey_i, float f_r, float f_i,
	int imin, int imax, int jmin, int jmax, int kmin, int kmax,
	int64_t adr0, param_t *p)
{
	for (int i = imin; i < imax + 1; i++) {
	for (int j = jmin; j < jmax + 0; j++) {
	for (int k = kmin; k < kmax + 1; k++) {
		const int64_t m = LA(p, i, j, k);
		assert((m >= 0) && (m < NN));
		const int64_t n = adr0 + m;
		dftadd(&cey_r[n], &cey_i[n], ey[m], f_r, f_i);
	}
	}
	}
}


// Ez
static void dft_near3dEz_cpu(
	real_t *ez, float *cez_r, float *cez_i, float f_r, float f_i,
	int imin, int imax, int jmin, int jmax, int kmin, int kmax,
	int64_t adr0, param_t *p)
{
	for (int i = imin; i < imax + 1; i++) {
	for (int j = jmin; j < jmax + 1; j++) {
	for (int k = kmin; k < kmax + 0; k++) {
		const int64_t m = LA(p, i, j, k);
		assert((m >= 0) && (m < NN));
		const int64_t n = adr0 + m;
		dftadd(&cez_r[n], &cez_i[n], ez[m], f_r, f_i);
	}
	}
	}
}


// Hx
static void dft_near3dHx_cpu(
	real_t *hx, float *chx_r, float *chx_i, float f_r, float f_i,
	int imin, int imax, int jmin, int jmax, int kmin, int kmax,
	int64_t adr0, param_t *p)
{
	for (int i = imin - 0; i < imax + 1; i++) {
	for (int j = jmin - 1; j < jmax + 1; j++) {
	for (int k = kmin - 1; k < kmax + 1; k++) {
		const int64_t m = LA(p, i, j, k);
		assert((m >= 0) && (m < NN));
		const int64_t n = adr0 + m;
		dftadd(&chx_r[n], &chx_i[n], hx[m], f_r, f_i);
	}
	}
	}
}


// Hy
static void dft_near3dHy_cpu(
	real_t *hy, float *chy_r, float *chy_i, float f_r, float f_i,
	int imin, int imax, int jmin, int jmax, int kmin, int kmax,
	int64_t adr0, param_t *p)
{
	for (int i = imin - 1; i < imax + 1; i++) {
	for (int j = jmin - 0; j < jmax + 1; j++) {
	for (int k = kmin - 1; k < kmax + 1; k++) {
		const int64_t m = LA(p, i, j, k);
		assert((m >= 0) && (m < NN));
		const int64_t n = adr0 + m;
		dftadd(&chy_r[n], &chy_i[n], hy[m], f_r, f_i);
	}
	}
	}
}


// Hz
static void dft_near3dHz_cpu(
	real_t *hz, float *chz_r, float *chz_i, float f_r, float f_i,
	int imin, int imax, int jmin, int jmax, int kmin, int kmax,
	int64_t adr0, param_t *p)
{
	for (int i = imin - 1; i < imax + 1; i++) {
	for (int j = jmin - 1; j < jmax + 1; j++) {
	for (int k = kmin - 0; k < kmax + 1; k++) {
		const int64_t m = LA(p, i, j, k);
		assert((m >= 0) && (m < NN));
		const int64_t n = adr0 + m;
		dftadd(&chz_r[n], &chz_i[n], hz[m], f_r, f_i);
	}
	}
	}
}

void dftNear3d(int itime)
{
	for (int ifreq = 0; ifreq < NFreq2; ifreq++) {
		const int64_t adr0 = ifreq * NN;
		const int id = (itime * NFreq2) + ifreq;

		const float fe_r = (float)cEdft[id].r;
		const float fe_i = (float)cEdft[id].i;
		const float fh_r = (float)cHdft[id].r;
		const float fh_i = (float)cHdft[id].i;

		if (GPU) {
			// Ex
			dim3 gridEx(
				CEIL(kMax - kMin + 1, updateBlock.x),
				CEIL(jMax - jMin + 1, updateBlock.y),
				CEIL(iMax - iMin + 0, updateBlock.z));
			dft_near3dEx_gpu<<<gridEx, updateBlock>>>(
				Ex, d_cEx_r, d_cEx_i, fe_r, fe_i, iMin, iMax, jMin, jMax, kMin, kMax, adr0, Ni, Nj, Nk, N0);

			// Ey
			dim3 gridEy(
				CEIL(kMax - kMin + 1, updateBlock.x),
				CEIL(jMax - jMin + 0, updateBlock.y),
				CEIL(iMax - iMin + 1, updateBlock.z));
			dft_near3dEy_gpu<<<gridEy, updateBlock>>>(
				Ey, d_cEy_r, d_cEy_i, fe_r, fe_i, iMin, iMax, jMin, jMax, kMin, kMax, adr0, Ni, Nj, Nk, N0);

			// Ez
			dim3 gridEz(
				CEIL(kMax - kMin + 0, updateBlock.x),
				CEIL(jMax - jMin + 1, updateBlock.y),
				CEIL(iMax - iMin + 1, updateBlock.z));
			dft_near3dEz_gpu<<<gridEz, updateBlock>>>(
				Ez, d_cEz_r, d_cEz_i, fe_r, fe_i, iMin, iMax, jMin, jMax, kMin, kMax, adr0, Ni, Nj, Nk, N0);

			// Hx
			dim3 gridHx(
				CEIL(kMax - kMin + 2, updateBlock.x),
				CEIL(jMax - jMin + 2, updateBlock.y),
				CEIL(iMax - iMin + 1, updateBlock.z));
			dft_near3dHx_gpu<<<gridHx, updateBlock>>>(
				Hx, d_cHx_r, d_cHx_i, fh_r, fh_i, iMin, iMax, jMin, jMax, kMin, kMax, adr0, Ni, Nj, Nk, N0);

			// Hy
			dim3 gridHy(
				CEIL(kMax - kMin + 2, updateBlock.x),
				CEIL(jMax - jMin + 1, updateBlock.y),
				CEIL(iMax - iMin + 2, updateBlock.z));
			dft_near3dHy_gpu<<<gridHy, updateBlock>>>(
				Hy, d_cHy_r, d_cHy_i, fh_r, fh_i, iMin, iMax, jMin, jMax, kMin, kMax, adr0, Ni, Nj, Nk, N0);

			// Hz
			dim3 gridHz(
				CEIL(kMax - kMin + 1, updateBlock.x),
				CEIL(jMax - jMin + 2, updateBlock.y),
				CEIL(iMax - iMin + 2, updateBlock.z));
			dft_near3dHz_gpu<<<gridHz, updateBlock>>>(
				Hz, d_cHz_r, d_cHz_i, fh_r, fh_i, iMin, iMax, jMin, jMax, kMin, kMax, adr0, Ni, Nj, Nk, N0);

			if (UM) cudaDeviceSynchronize();
		}
		else {
			// CPU TODO
			dft_near3dEx_cpu(Ex, d_cEx_r, d_cEx_i, fe_r, fe_i, iMin, iMax, jMin, jMax, kMin, kMax, adr0, &h_Param);
			dft_near3dEy_cpu(Ey, d_cEy_r, d_cEy_i, fe_r, fe_i, iMin, iMax, jMin, jMax, kMin, kMax, adr0, &h_Param);
			dft_near3dEz_cpu(Ez, d_cEz_r, d_cEz_i, fe_r, fe_i, iMin, iMax, jMin, jMax, kMin, kMax, adr0, &h_Param);
			dft_near3dHx_cpu(Hx, d_cHx_r, d_cHx_i, fh_r, fh_i, iMin, iMax, jMin, jMax, kMin, kMax, adr0, &h_Param);
			dft_near3dHy_cpu(Hy, d_cHy_r, d_cHy_i, fh_r, fh_i, iMin, iMax, jMin, jMax, kMin, kMax, adr0, &h_Param);
			dft_near3dHz_cpu(Hz, d_cHz_r, d_cHz_i, fh_r, fh_i, iMin, iMax, jMin, jMax, kMin, kMax, adr0, &h_Param);
		}
	}
}
