/*
efeed.cu

E on feeds
*/

#include "orcwa.h"
#include "orcwa_cuda.h"

// feed voltage : gauss derivative
__host__ __device__
static double vfeed_(double t, double tw, double td)
{
	double arg = (t - tw - td) / (tw / 4.0);
	double v = sqrt(2.0) * exp(0.5) * arg * exp(-arg * arg);

	return v;
}

__host__ __device__
static void efeed_(
	feed_t feed, double *v, double *c,
	double t, double tw, double rfeed, double eta0,
	real_t *ex, real_t *ey, real_t *ez, real_t *hx, real_t *hy, real_t *hz, id_t *iex, id_t *iey, id_t *iez,
	param_t *p)
{
	const double eps= 1e-6;

	// V
	const double v0 = vfeed_(t, tw, feed.delay);
	*v = v0 * feed.volt;

	const int i = feed.i;
	const int j = feed.j;
	const int k = feed.k;
	const double dx = feed.dx;
	const double dy = feed.dy;
	const double dz = feed.dz;

	// E and I
	if      ((feed.dir == 'X') &&
	         (p->iMin <= i) && (i <  p->iMax) &&
	         (p->jMin <= j) && (j <= p->jMax) &&
	         (p->kMin <= k) && (k <= p->kMax)) {  // MPI
		*c = dz * (hz[LA(p, i, j, k)] - hz[LA(p, i,     j - 1, k    )])
		   - dy * (hy[LA(p, i, j, k)] - hy[LA(p, i,     j,     k - 1)]);
		*c /= eta0;
		*v -= rfeed * (*c);
		if ((iex[LA(p, i, j, k)] == PEC) || (fabs(v0) > eps)) {
			ex[LA(p, i, j, k)] = -(real_t)(*v / dx);
		}
	}
	else if ((feed.dir == 'Y') &&
	         (p->iMin <= i) && (i <= p->iMax) &&
	         (p->jMin <= j) && (j <  p->jMax) &&
	         (p->kMin <= k) && (k <= p->kMax)) {  // MPI
		*c = dx * (hx[LA(p, i, j, k)] - hx[LA(p, i,     j,     k - 1)])
		   - dz * (hz[LA(p, i, j, k)] - hz[LA(p, i - 1, j,     k    )]);
		*c /= eta0;
		*v -= rfeed * (*c);
		if ((iey[LA(p, i, j, k)] == PEC) || (fabs(v0) > eps)) {
			ey[LA(p, i, j, k)] = -(real_t)(*v / dy);
		}
	}
	else if ((feed.dir == 'Z') &&
	         (p->iMin <= i) && (i <= p->iMax) &&
	         (p->jMin <= j) && (j <= p->jMax) &&
	         (p->kMin <= k) && (k <  p->kMax)) {  // MPI
		*c = dy * (hy[LA(p, i, j, k)] - hy[LA(p, i - 1, j,     k    )])
		   - dx * (hx[LA(p, i, j, k)] - hx[LA(p, i,     j - 1, k    )]);
		*c /= eta0;
		*v -= rfeed * (*c);
		if ((iez[LA(p, i, j, k)] == PEC) || (fabs(v0) > eps)) {
			ez[LA(p, i, j, k)] = -(real_t)(*v / dz);
		}
	}
}

// gpu
__global__
static void efeed_gpu(
	int nfeed, feed_t *feed, double *vin, double *iin,
	double t, double tw, int itime, int maxiter, double rfeed, double eta0,
	real_t *ex, real_t *ey, real_t *ez, real_t *hx, real_t *hy, real_t *hz, id_t *iex, id_t *iey, id_t *iez)
{
	const int n = threadIdx.x + (blockIdx.x * blockDim.x);
	if (n < nfeed) {
		const int adr = n * (maxiter + 1) + itime;
		efeed_(
			feed[n], &vin[adr], &iin[adr],
			t, tw, rfeed, eta0,
			ex, ey, ez, hx, hy, hz, iex, iey, iez,
			&d_Param);
	}
}

// cpu
static void efeed_cpu(
	int nfeed, feed_t *feed, double *vin, double *iin,
	double t, double tw, int itime, int maxiter, double rfeed, double eta0,
	real_t *ex, real_t *ey, real_t *ez, real_t *hx, real_t *hy, real_t *hz, id_t *iex, id_t *iey, id_t *iez)
{
	for (int n = 0; n < nfeed; n++) {
		const int adr = n * (maxiter + 1) + itime;
		efeed_(
			feed[n], &vin[adr], &iin[adr],
			t, tw, rfeed, eta0,
			ex, ey, ez, hx, hy, hz, iex, iey, iez,
			&h_Param);
	}
}

void efeed(int itime)
{
	if (NFeed <= 0) return;

	const double t = (itime + 1) * Dt;

	if (GPU) {
		cudaMemcpyToSymbol(d_Param, &h_Param, sizeof(param_t));
		const int block = MIN(256, NFeed);
		const int grid = CEIL(NFeed, block);
		efeed_gpu<<<grid, block>>>(
			NFeed, d_Feed, d_VFeed, d_IFeed,
			t, Tw, itime, Solver.maxiter, rFeed, ETA0,
			Ex, Ey, Ez, Hx, Hy, Hz, d_iEx, d_iEy, d_iEz);
		if (UM) cudaDeviceSynchronize();
	}
	else {
		efeed_cpu(
			NFeed, Feed, VFeed, IFeed,
			t, Tw, itime, Solver.maxiter, rFeed, ETA0,
			Ex, Ey, Ez, Hx, Hy, Hz, iEx, iEy, iEz);
	}
}
