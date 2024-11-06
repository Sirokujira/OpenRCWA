/*
eload.c

E on loads (inductors)
*/

#include "orcwa.h"
#include "orcwa_prototype.h"

void eload(void)
{
	if (NInductor <= 0) return;

	const double cdt = C * Dt;

	for (int n = 0; n < NInductor; n++) {
		inductor_t *ptr = &Inductor[n];

		int     i = ptr->i;
		int     j = ptr->j;
		int     k = ptr->k;
		double dx = ptr->dx;
		double dy = ptr->dy;
		double dz = ptr->dz;

		if      ((ptr->dir == 'X') &&
		         (iMin <= i) && (i <  iMax) &&
		         (jMin <= j) && (j <= jMax) &&
		         (kMin <= k) && (k <= kMax)) {  // MPI
			const double roth = (HZ(i, j, k) - HZ(i,     j - 1, k    )) / dy
			                  - (HY(i, j, k) - HY(i,     j,     k - 1)) / dz;
			EX(i, j, k) = (real_t)(ptr->e + (cdt * roth) - (ptr->fctr * cdt * cdt * ptr->esum));
			ptr->e = EX(i, j, k);
			ptr->esum += ptr->e;
		}
		else if ((ptr->dir == 'Y') &&
		         (iMin <= i) && (i <= iMax) &&
		         (jMin <= j) && (j <  jMax) &&
		         (kMin <= k) && (k <= kMax)) {  // MPI
			const double roth = (HX(i, j, k) - HX(i,     j,     k - 1)) / dz
			                  - (HZ(i, j, k) - HZ(i - 1, j,     k    )) / dx;
			EY(i, j, k) = (real_t)(ptr->e + (cdt * roth) - (ptr->fctr * cdt * cdt * ptr->esum));
			ptr->e = EY(i, j, k);
			ptr->esum += ptr->e;
		}
		else if ((ptr->dir == 'Z') &&
		         (iMin <= i) && (i <= iMax) &&
		         (jMin <= j) && (j <= jMax) &&
		         (kMin <= k) && (k <  kMax)) {  // MPI
			const double roth = (HY(i, j, k) - HY(i - 1, j,     k    )) / dx
			                  - (HX(i, j, k) - HX(i,     j - 1, k    )) / dy;
			EZ(i, j, k) = (real_t)(ptr->e + (cdt * roth) - (ptr->fctr * cdt * cdt * ptr->esum));
			ptr->e = EZ(i, j, k);
			ptr->esum += ptr->e;
		}
	}
}
