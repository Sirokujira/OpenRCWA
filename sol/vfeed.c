/*
vfeed.c

feed voltage : gauss derivative
*/

#include <math.h>

double vfeed(double t, double tw, double td)
{
	const double arg = (t - tw - td) / (tw / 4);

	return sqrt(2.0) * exp(0.5) * arg * exp(-arg * arg);
}
