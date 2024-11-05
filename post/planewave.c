/*
planewave.c
*/

#include "orcwa.h"
#include "complex.h"

void planewave(double freq, double x, double y, double z, d_complex_t e[], d_complex_t h[])
{
	if (!IPlanewave) return;

	const double x0 = (Xn[0] + Xn[Nx]) / 2;
	const double y0 = (Yn[0] + Yn[Ny]) / 2;
	const double z0 = (Zn[0] + Zn[Nz]) / 2;

	const double rri = (x - x0) * Planewave.ri[0]
	                 + (y - y0) * Planewave.ri[1]
	                 + (z - z0) * Planewave.ri[2];

	const double k = (2 * PI * freq) / C;

	const d_complex_t phs = d_exp(-k * rri);

	for (int m = 0; m < 3; m++) {
		e[m] = d_rmul(Planewave.ei[m], phs);
		h[m] = d_rmul(Planewave.hi[m], phs);
	}
}
