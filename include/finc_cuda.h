/*
finc_cuda.h (CUDA)

incidence function (gauss derivative)
*/

__host__ __device__ __forceinline__
static void finc_cuda(
	real_t x, real_t y, real_t z, real_t t,
	const real_t r0[], const real_t ri[], real_t fc, real_t ai, real_t dt,
	real_t *fi, real_t *dfi)
{
	const real_t c = (real_t)2.99792458e8;

	t -= ((x - r0[0]) * ri[0]
	    + (y - r0[1]) * ri[1]
	    + (z - r0[2]) * ri[2]) / c;

	const real_t at = ai * t;
	const real_t ex = (at * at < 16) ? (real_t)exp(-at * at) : 0;
	//const real_t ex = (real_t)exp(-at * at);
	*fi = at * ex * fc;
	*dfi = dt * ai * (1 - 2 * at * at) * ex * fc;
}
