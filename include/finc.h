/*
finc.h

incidence function (gauss derivative)
*/

#include <math.h>

static inline void finc(
	double x, double y, double z, double t,
	const double r0[], const double ri[], double fc, double ai, double dt,
	real_t *fi, real_t *dfi)
{
	const double c = 2.99792458e8;

	t -= ((x - r0[0]) * ri[0]
	    + (y - r0[1]) * ri[1]
	    + (z - r0[2]) * ri[2]) / c;

	const double at = ai * t;
	const double ex = (at * at < 16) ? exp(-at * at) : 0;
	//const double ex = exp(-at * at);
	*fi = (real_t)(at * ex * fc);
	*dfi = (real_t)(dt * ai * (1 - 2 * at * at) * ex * fc);
}
