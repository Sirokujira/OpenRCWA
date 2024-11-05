/*
vpoint.c

V waveform on points
*/

#include "orcwa.h"
#include "finc.h"

void vpoint(int itime)
{
	if (NPoint <= 0) return;
	real_t fi, dfi;

	for (int n = 0; n < NPoint + 2; n++) {
		const int i = Point[n].i;
		const int j = Point[n].j;
		const int k = Point[n].k;

		double e = 0;
		double d = 0;
		if      ((Point[n].dir == 'X') &&
		         (iMin <= i) && (i <  iMax) &&
		         (jMin <= j) && (j <= jMax) &&
		         (kMin <= k) && (k <= kMax)) {  // MPI
			e = EX(i, j, k);
			d = Point[n].dx;
			if (IPlanewave) {
				const double t = (itime + 1) * Dt;
				finc(Xc[i], Yn[j], Zn[k], t, Planewave.r0, Planewave.ri, Planewave.ei[0], Planewave.ai, Dt, &fi, &dfi);
				e += fi;
			}
		}
		else if ((Point[n].dir == 'Y') &&
		         (iMin <= i) && (i <= iMax) &&
		         (jMin <= j) && (j <  jMax) &&
		         (kMin <= k) && (k <= kMax)) {  // MPI
			e = EY(i, j, k);
			d = Point[n].dy;
			if (IPlanewave) {
				const double t = (itime + 1) * Dt;
				finc(Xn[i], Yc[j], Zn[k], t, Planewave.r0, Planewave.ri, Planewave.ei[1], Planewave.ai, Dt, &fi, &dfi);
				e += fi;
			}
		}
		else if ((Point[n].dir == 'Z') &&
		         (iMin <= i) && (i <= iMax) &&
		         (jMin <= j) && (j <= jMax) &&
		         (kMin <= k) && (k <  kMax)) {  // MPI
			e = EZ(i, j, k);
			d = Point[n].dz;
			if (IPlanewave) {
				const double t = (itime + 1) * Dt;
				finc(Xn[i], Yn[j], Zc[k], t, Planewave.r0, Planewave.ri, Planewave.ei[2], Planewave.ai, Dt, &fi, &dfi);
				e += fi;
			}
		}
		const int id = n * (Solver.maxiter + 1) + itime;
		VPoint[id] = e * (-d);
	}
}
