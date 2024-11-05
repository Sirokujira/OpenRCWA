/*
dispersionEz.cu (CUDA)

update Ez (dispersion)
*/

#include "orcwa.h"
#include "orcwa_cuda.h"
#include "finc_cuda.h"


__host__ __device__
static void dispersion(
	real_t e[], real_t *de, dispersion_t *me, param_t *p,
	real_t x, real_t y, real_t z, real_t t)
{
	const int64_t n = LA(p, me->i, me->j, me->k);

	real_t fi = 0;
	if (p->IPlanewave) {
		real_t dfi;
		finc_cuda(x, y, z, t, p->r0, p->ri, p->ei[2], p->ai, p->dt, &fi, &dfi);
	}

	e[n] += me->f1 * (*de);

	*de = me->f2 * (e[n] + fi)
	    + me->f3 * (*de);
}


__global__
static void dispersionEz_gpu(
	int64_t num, real_t e[], real_t de[], dispersion_t me[],
	const real_t xn[], const real_t yn[], const real_t zc[], real_t t)
	
{
	const int64_t n = threadIdx.x + (blockIdx.x * blockDim.x);
	if (n < num) {
		const real_t x = xn[me[n].i];
		const real_t y = yn[me[n].j];
		const real_t z = zc[me[n].k];
		dispersion(
			e, &de[n], &me[n], &d_Param,
			x, y, z, t);
	}
}


static void dispersionEz_cpu(
	int64_t num, real_t e[], real_t de[], dispersion_t me[],
	const real_t xn[], const real_t yn[], const real_t zc[], real_t t)
{
	for (int64_t n = 0; n < num; n++) {
		const real_t x = xn[me[n].i];
		const real_t y = yn[me[n].j];
		const real_t z = zc[me[n].k];
		dispersion(
			e, &de[n], &me[n], &h_Param,
			x, y, z, t);
	}
}


void dispersionEz(double t)
{
	if (GPU) {
		cudaMemcpyToSymbol(d_Param, &h_Param, sizeof(param_t));
		dispersionEz_gpu<<<(int)CEIL(numDispersionEz, dispersionBlock), dispersionBlock>>>(
			numDispersionEz, Ez, d_DispersionEz, d_mDispersionEz,
			d_Xn, d_Yn, d_Zc, (real_t)t);
		if (UM) cudaDeviceSynchronize();
	}
	else {
		dispersionEz_cpu(
			numDispersionEz, Ez, DispersionEz, mDispersionEz,
			h_Xn, h_Yn, h_Zc, (real_t)t);
	}
}
