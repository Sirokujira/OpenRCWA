/*
farfield.c

far field
*/

#include "orcwa.h"
#include "complex.h"
#include "orcwa_prototype.h"

void farfield(int ifreq, double theta, double phi, double ffctr, d_complex_t *etheta, d_complex_t *ephi)
{
	// wave number

	const double kwave = (2 * PI * Freq2[ifreq]) / C;

	// unit vector in r, theta, phi

	double r1[3], t1[3], p1[3];

	const double cos_t = cos(theta * DTOR);
	const double sin_t = sin(theta * DTOR);
	const double cos_p = cos(phi   * DTOR);
	const double sin_p = sin(phi   * DTOR);

	r1[0] = +sin_t * cos_p;
	r1[1] = +sin_t * sin_p;
	r1[2] = +cos_t;

	t1[0] = +cos_t * cos_p;
	t1[1] = +cos_t * sin_p;
	t1[2] = -sin_t;

	p1[0] = -sin_p;
	p1[1] = +cos_p;
	p1[2] = 0;

	d_complex_t plx = d_complex(0, 0);
	d_complex_t ply = d_complex(0, 0);
	d_complex_t plz = d_complex(0, 0);
	d_complex_t pnx = d_complex(0, 0);
	d_complex_t pny = d_complex(0, 0);
	d_complex_t pnz = d_complex(0, 0);

	d_complex_t *ex = SurfaceEx[ifreq];
	d_complex_t *ey = SurfaceEy[ifreq];
	d_complex_t *ez = SurfaceEz[ifreq];
	d_complex_t *hx = SurfaceHx[ifreq];
	d_complex_t *hy = SurfaceHy[ifreq];
	d_complex_t *hz = SurfaceHz[ifreq];

	for (int n = 0; n < NSurface; n++) {
		surface_t *p = &Surface[n];

		// Z0 * J = n X (Z0 * H)
		const d_complex_t cjx = d_sub(d_rmul(+p->ny, hz[n]), d_rmul(+p->nz, hy[n]));
		const d_complex_t cjy = d_sub(d_rmul(+p->nz, hx[n]), d_rmul(+p->nx, hz[n]));
		const d_complex_t cjz = d_sub(d_rmul(+p->nx, hy[n]), d_rmul(+p->ny, hx[n]));

		// M = -n X E
		const d_complex_t cmx = d_sub(d_rmul(-p->ny, ez[n]), d_rmul(-p->nz, ey[n]));
		const d_complex_t cmy = d_sub(d_rmul(-p->nz, ex[n]), d_rmul(-p->nx, ez[n]));
		const d_complex_t cmz = d_sub(d_rmul(-p->nx, ey[n]), d_rmul(-p->ny, ex[n]));

		// exp(jkr * r) * dS
		const double rr = (r1[0] * p->x) + (r1[1] * p->y) + (r1[2] * p->z);
		const d_complex_t expds = d_rmul(p->ds, d_exp(kwave * rr));

		// L += M * exp(jkr * r) * dS
		plx = d_add(plx, d_mul(cmx, expds));
		ply = d_add(ply, d_mul(cmy, expds));
		plz = d_add(plz, d_mul(cmz, expds));

		// Z0 * N += (Z0 * J) * exp(jkr * r) * dS
		pnx = d_add(pnx, d_mul(cjx, expds));
		pny = d_add(pny, d_mul(cjy, expds));
		pnz = d_add(pnz, d_mul(cjz, expds));
	}

	// Z0 * N-theta, Z0 * N-phi
	const d_complex_t pnt = d_add3(d_rmul(t1[0], pnx), d_rmul(t1[1], pny), d_rmul(t1[2], pnz));
	const d_complex_t pnp = d_add3(d_rmul(p1[0], pnx), d_rmul(p1[1], pny), d_rmul(p1[2], pnz));

	// L-theta, L-phi
	const d_complex_t plt = d_add3(d_rmul(t1[0], plx), d_rmul(t1[1], ply), d_rmul(t1[2], plz));
	const d_complex_t plp = d_add3(d_rmul(p1[0], plx), d_rmul(p1[1], ply), d_rmul(p1[2], plz));

	// F-theta, F-phi
	*etheta = d_rmul(ffctr, d_add(pnt, plp));
	*ephi   = d_rmul(ffctr, d_sub(pnp, plt));
}


// alloc far field array
void alloc_farfield(void)
{
	size_t size;

	assert((Nx > 0) && (Ny > 0) && (Nz > 0) && (NFreq2 > 0));

	NSurface = 2 * ((Nx * Ny) + (Ny * Nz) + (Nz * Nx));
	//printf("%zd %d\n", NSurface, NFreq2);

	size = NSurface * sizeof(surface_t);
	Surface = (surface_t *)malloc(size);

	size = NFreq2 * sizeof(d_complex_t *);
	SurfaceEx = (d_complex_t **)malloc(size);
	SurfaceEy = (d_complex_t **)malloc(size);
	SurfaceEz = (d_complex_t **)malloc(size);
	SurfaceHx = (d_complex_t **)malloc(size);
	SurfaceHy = (d_complex_t **)malloc(size);
	SurfaceHz = (d_complex_t **)malloc(size);

	size = NSurface * sizeof(d_complex_t);
	for (int ifreq = 0; ifreq < NFreq2; ifreq++) {
		SurfaceEx[ifreq] = (d_complex_t *)malloc(size);
		SurfaceEy[ifreq] = (d_complex_t *)malloc(size);
		SurfaceEz[ifreq] = (d_complex_t *)malloc(size);
		SurfaceHx[ifreq] = (d_complex_t *)malloc(size);
		SurfaceHy[ifreq] = (d_complex_t *)malloc(size);
		SurfaceHz[ifreq] = (d_complex_t *)malloc(size);
	}
}


// setup surface E and H
void setup_farfield(void)
{
	d_complex_t cex[2][2], cey[2][2], cez[2][2];
	d_complex_t chx[2][2], chy[2][2], chz[2][2];

	// counter
	int64_t n = 0;

	// X surface
	for (int side = 0; side < 2; side++) {
		const int i = (side == 0) ? 0 : Nx;
		for (int j = 0; j < Ny; j++) {
		for (int k = 0; k < Nz; k++) {
			for (int ifreq = 0; ifreq < NFreq2; ifreq++) {
				for (int jn = 0; jn < 2; jn++) {
				for (int kn = 0; kn < 2; kn++) {
					NodeE_c(ifreq, i, j + jn, k + kn, &cex[jn][kn], &cey[jn][kn], &cez[jn][kn]);
					NodeH_c(ifreq, i, j + jn, k + kn, &chx[jn][kn], &chy[jn][kn], &chz[jn][kn]);
				}
				}
				SurfaceEx[ifreq][n] = d_complex(0, 0);
				SurfaceEy[ifreq][n] = d_rmul(0.25, d_add4(cey[0][0],
				                                          cey[0][1],
				                                          cey[1][0],
				                                          cey[1][1]));
				SurfaceEz[ifreq][n] = d_rmul(0.25, d_add4(cez[0][0],
				                                          cez[0][1],
				                                          cez[1][0],
				                                          cez[1][1]));
				SurfaceHx[ifreq][n] = d_complex(0, 0);
				SurfaceHy[ifreq][n] = d_rmul(0.25, d_add4(chy[0][0],
				                                          chy[0][1],
				                                          chy[1][0],
				                                          chy[1][1]));
				SurfaceHz[ifreq][n] = d_rmul(0.25, d_add4(chz[0][0],
				                                          chz[0][1],
				                                          chz[1][0],
				                                          chz[1][1]));
			}
			Surface[n].nx = (side == 0) ? -1 : +1;
			Surface[n].ny = 0;
			Surface[n].nz = 0;
			Surface[n].x = Xn[i];
			Surface[n].y = Yc[j];
			Surface[n].z = Zc[k];
			Surface[n].ds = (Yn[j + 1] - Yn[j]) * (Zn[k + 1] - Zn[k]);
			if (PBCx) Surface[n].ds = 0;  // X PBC -> skip X boundaries
			n++;
		}
		}
	}

	// Y surface
	for (int side = 0; side < 2; side++) {
		const int j = (side == 0) ? 0 : Ny;
		for (int k = 0; k < Nz; k++) {
		for (int i = 0; i < Nx; i++) {
			for (int ifreq = 0; ifreq < NFreq2; ifreq++) {
				for (int kn = 0; kn < 2; kn++) {
				for (int in = 0; in < 2; in++) {
					NodeE_c(ifreq, i + in, j, k + kn, &cex[kn][in], &cey[kn][in], &cez[kn][in]);
					NodeH_c(ifreq, i + in, j, k + kn, &chx[kn][in], &chy[kn][in], &chz[kn][in]);
				}
				}
				SurfaceEx[ifreq][n] = d_rmul(0.25, d_add4(cex[0][0],
				                                          cex[0][1],
				                                          cex[1][0],
				                                          cex[1][1]));
				SurfaceEy[ifreq][n] = d_complex(0, 0);
				SurfaceEz[ifreq][n] = d_rmul(0.25, d_add4(cez[0][0],
				                                          cez[0][1],
				                                          cez[1][0],
				                                          cez[1][1]));
				SurfaceHx[ifreq][n] = d_rmul(0.25, d_add4(chx[0][0],
				                                          chx[0][1],
				                                          chx[1][0],
				                                          chx[1][1]));
				SurfaceHy[ifreq][n] = d_complex(0, 0);
				SurfaceHz[ifreq][n] = d_rmul(0.25, d_add4(chz[0][0],
				                                          chz[0][1],
				                                          chz[1][0],
				                                          chz[1][1]));
			}
			Surface[n].nx = 0;
			Surface[n].ny = (side == 0) ? -1 : +1;
			Surface[n].nz = 0;
			Surface[n].x = Xc[i];
			Surface[n].y = Yn[j];
			Surface[n].z = Zc[k];
			Surface[n].ds = (Zn[k + 1] - Zn[k]) * (Xn[i + 1] - Xn[i]);
			if (PBCy) Surface[n].ds = 0;  // Y PBC -> skip Y boundaries
			n++;
		}
		}
	}

	// Z surface
	for (int side = 0; side < 2; side++) {
		const int k = (side == 0) ? 0 : Nz;
		for (int i = 0; i < Nx; i++) {
		for (int j = 0; j < Ny; j++) {
			for (int ifreq = 0; ifreq < NFreq2; ifreq++) {
				for (int in = 0; in < 2; in++) {
				for (int jn = 0; jn < 2; jn++) {
					NodeE_c(ifreq, i + in, j + jn, k, &cex[in][jn], &cey[in][jn], &cez[in][jn]);
					NodeH_c(ifreq, i + in, j + jn, k, &chx[in][jn], &chy[in][jn], &chz[in][jn]);
				}
				}
				SurfaceEx[ifreq][n] = d_rmul(0.25, d_add4(cex[0][0],
				                                          cex[0][1],
				                                          cex[1][0],
				                                          cex[1][1]));
				SurfaceEy[ifreq][n] = d_rmul(0.25, d_add4(cey[0][0],
				                                          cey[0][1],
				                                          cey[1][0],
				                                          cey[1][1]));
				SurfaceEz[ifreq][n] = d_complex(0, 0);
				SurfaceHx[ifreq][n] = d_rmul(0.25, d_add4(chx[0][0],
				                                          chx[0][1],
				                                          chx[1][0],
				                                          chx[1][1]));
				SurfaceHy[ifreq][n] = d_rmul(0.25, d_add4(chy[0][0],
				                                          chy[0][1],
				                                          chy[1][0],
				                                          chy[1][1]));
				SurfaceHz[ifreq][n] = d_complex(0, 0);
			}
			Surface[n].nx = 0;
			Surface[n].ny = 0;
			Surface[n].nz = (side == 0) ? -1 : +1;
			Surface[n].x = Xc[i];
			Surface[n].y = Yc[j];
			Surface[n].z = Zn[k];
			Surface[n].ds = (Xn[i + 1] - Xn[i]) * (Yn[j + 1] - Yn[j]);
			if (PBCz) Surface[n].ds = 0;  // Z PBC -> skip Z boundaries
			n++;
		}
		}
	}
}


// far field components
void farComponent(d_complex_t etheta, d_complex_t ephi, double e[])
{
	// abs
	e[0] = sqrt(d_norm(etheta) + d_norm(ephi));

	// theta/phi
	e[1] = d_abs(etheta);
	e[2] = d_abs(ephi);

	// major/minor
	double tmp = d_abs(d_add(d_mul(etheta, etheta), d_mul(ephi, ephi)));
	e[3] = sqrt((d_norm(etheta) + d_norm(ephi) + tmp) / 2);
	e[4] = sqrt((d_norm(etheta) + d_norm(ephi) - tmp) / 2);

	// RHCP/LHCP
	e[5] = d_abs(d_add(etheta, d_mul(d_complex(0, 1), ephi))) / sqrt(2);
	e[6] = d_abs(d_sub(etheta, d_mul(d_complex(0, 1), ephi))) / sqrt(2);
}


// far field factor
double farfactor(int ifreq)
{
	double ffctr = 0;

	const double kwave = (2 * PI * Freq2[ifreq]) / C;

	if (NFeed) {
		// feed (post only)
		double sum = 0;
		for (int ifeed = 0; ifeed < NFeed; ifeed++) {
			sum += 0.5 * Pin[MatchingLoss ? 1 : 0][(ifeed * NFreq2) + ifreq];
		}
		ffctr = kwave / sqrt(8 * PI * ETA0 * sum);
	}
	else if (IPlanewave) {
		// plane wave (solver + post)
		const double einc = 1;
		ffctr = kwave / (einc * sqrt(4 * PI));
	}

	return ffctr;
}
