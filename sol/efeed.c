/*
efeed.c

E on feeds
*/

#include "orcwa.h"
#include "orcwa_prototype.h"

void efeed(int itime)
{
	if (NFeed <= 0) return;

	const double eps = 1e-6;

	const double t = (itime + 1) * Dt;

	for (int ifeed = 0; ifeed < NFeed; ifeed++) {
		const int     i = Feed[ifeed].i;
		const int     j = Feed[ifeed].j;
		const int     k = Feed[ifeed].k;
		const double dx = Feed[ifeed].dx;
		const double dy = Feed[ifeed].dy;
		const double dz = Feed[ifeed].dz;

		// V
		const double v0 = vfeed(t, Tw, Feed[ifeed].delay);
		double v = v0 * Feed[ifeed].volt;

		// E, V, I
		double c = 0;
		if      ((Feed[ifeed].dir == 'X') &&
		         (iMin <= i) && (i <  iMax) &&
		         (jMin <= j) && (j <= jMax) &&
		         (kMin <= k) && (k <= kMax)) {  // MPI
			c = dz * (HZ(i, j, k) - HZ(i,     j - 1, k    ))
			  - dy * (HY(i, j, k) - HY(i,     j,     k - 1));
			c /= ETA0;
			v -= rFeed * c;
			if ((IEX(i, j, k) == PEC) || (fabs(v0) > eps)) {
				EX(i, j, k) = -(real_t)(v / dx);
			}
		}
		else if ((Feed[ifeed].dir == 'Y') &&
		         (iMin <= i) && (i <= iMax) &&
		         (jMin <= j) && (j <  jMax) &&
		         (kMin <= k) && (k <= kMax)) {  // MPI
			c = dx * (HX(i, j, k) - HX(i,     j,     k - 1))
			  - dz * (HZ(i, j, k) - HZ(i - 1, j,     k    ));
			c /= ETA0;
			v -= rFeed * c;
			if ((IEY(i, j, k) == PEC) || (fabs(v0) > eps)) {
				EY(i, j, k) = -(real_t)(v / dy);
			}
		}
		else if ((Feed[ifeed].dir == 'Z') &&
		         (iMin <= i) && (i <= iMax) &&
		         (jMin <= j) && (j <= jMax) &&
		         (kMin <= k) && (k <  kMax)) {  // MPI
			c = dy * (HY(i, j, k) - HY(i - 1, j,     k    ))
			  - dx * (HX(i, j, k) - HX(i,     j - 1, k    ));
			c /= ETA0;
			v -= rFeed * c;
			if ((IEZ(i, j, k) == PEC) || (fabs(v0) > eps)) {
				EZ(i, j, k) = -(real_t)(v / dz);
			}
		}

		// V/I waveform
		const int id = ifeed * (Solver.maxiter + 1) + itime;
		VFeed[id] = v;
		IFeed[id] = c;
	}
}
