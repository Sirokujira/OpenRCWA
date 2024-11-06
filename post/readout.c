/*
readout.c

read orcwa.out
*/

#include "orcwa.h"

void readout(FILE *fp)
{
	size_t num = 0;
	size_t size;

	num += fread(Title,           sizeof(char),        256, fp);
	num += fread(&Nx,             sizeof(int),         1,   fp);
	num += fread(&Ny,             sizeof(int),         1,   fp);
	num += fread(&Nz,             sizeof(int),         1,   fp);
	num += fread(&Ni,             sizeof(int),         1,                             fp);
	num += fread(&Nj,             sizeof(int),         1,                             fp);
	num += fread(&Nk,             sizeof(int),         1,                             fp);
	num += fread(&N0,             sizeof(int),         1,                             fp);
	num += fread(&NN,             sizeof(int64_t),     1,   fp);
	num += fread(&NFreq1,         sizeof(int),         1,   fp);
	num += fread(&NFreq2,         sizeof(int),         1,   fp);
	num += fread(&NFeed,          sizeof(int),         1,   fp);
	num += fread(&NPoint,         sizeof(int),         1,   fp);
	num += fread(&Niter,          sizeof(int),         1,   fp);
	num += fread(&Ntime,          sizeof(int),         1,   fp);
	num += fread(&Solver.maxiter, sizeof(int),         1,   fp);
	num += fread(&Solver.nout,    sizeof(int),         1,   fp);
	num += fread(&Dt,             sizeof(double),      1,   fp);
	num += fread(&NGline,         sizeof(int),         1,   fp);
	num += fread(&IPlanewave,     sizeof(int),         1,   fp);
	num += fread(&Planewave,      sizeof(planewave_t), 1,   fp);

	Xn     =         (double *)malloc(sizeof(double)      * (Nx + 1));
	Yn     =         (double *)malloc(sizeof(double)      * (Ny + 1));
	Zn     =         (double *)malloc(sizeof(double)      * (Nz + 1));
	Xc     =         (double *)malloc(sizeof(double)      * Nx);
	Yc     =         (double *)malloc(sizeof(double)      * Ny);
	Zc     =         (double *)malloc(sizeof(double)      * Nz);
	Eiter  =         (double *)malloc(sizeof(double)      * Niter);
	Hiter  =         (double *)malloc(sizeof(double)      * Niter);
	VFeed  =         (double *)malloc(sizeof(double)      * NFeed  * (Solver.maxiter + 1));
	IFeed  =         (double *)malloc(sizeof(double)      * NFeed  * (Solver.maxiter + 1));
	VPoint =         (double *)malloc(sizeof(double)      * NPoint * (Solver.maxiter + 1));
	Freq1  =         (double *)malloc(sizeof(double)      * NFreq1);
	Freq2  =         (double *)malloc(sizeof(double)      * NFreq2);
	Feed   =         (feed_t *)malloc(sizeof(feed_t)      * NFeed);
	Zin    =    (d_complex_t *)malloc(sizeof(d_complex_t) * NFeed * NFreq1);
	Ref    =         (double *)malloc(sizeof(double)      * NFeed * NFreq1);
	Pin[0] =         (double *)malloc(sizeof(double)      * NFeed * NFreq2);
	Pin[1] =         (double *)malloc(sizeof(double)      * NFeed * NFreq2);
	Spara  =    (d_complex_t *)malloc(sizeof(d_complex_t) * NPoint * NFreq1);
	Gline  = (double (*)[2][3])malloc(sizeof(double)      * NGline * 2 * 3);

	num += fread(Xn,     sizeof(double),      Nx + 1,                        fp);
	num += fread(Yn,     sizeof(double),      Ny + 1,                        fp);
	num += fread(Zn,     sizeof(double),      Nz + 1,                        fp);
	num += fread(Xc,     sizeof(double),      Nx,                            fp);
	num += fread(Yc,     sizeof(double),      Ny,                            fp);
	num += fread(Zc,     sizeof(double),      Nz,                            fp);
	num += fread(Eiter,  sizeof(double),      Niter,                         fp);
	num += fread(Hiter,  sizeof(double),      Niter,                         fp);
	num += fread(VFeed,  sizeof(double),      NFeed  * (Solver.maxiter + 1), fp);
	num += fread(IFeed,  sizeof(double),      NFeed  * (Solver.maxiter + 1), fp);
	num += fread(VPoint, sizeof(double),      NPoint * (Solver.maxiter + 1), fp);
	num += fread(Freq1,  sizeof(double),      NFreq1,                        fp);
	num += fread(Freq2,  sizeof(double),      NFreq2,                        fp);
	num += fread(Feed,   sizeof(feed_t),      NFeed,                         fp);
	num += fread(Zin,    sizeof(d_complex_t), NFeed * NFreq1,                fp);
	num += fread(Ref,    sizeof(double),      NFeed * NFreq1,                fp);
	num += fread(Pin[0], sizeof(double),      NFeed * NFreq2,                fp);
	num += fread(Pin[1], sizeof(double),      NFeed * NFreq2,                fp);
	num += fread(Spara,  sizeof(d_complex_t), NPoint * NFreq1,               fp);
	num += fread(Gline,  sizeof(double),      NGline * 2 * 3,                fp);

	num += fread(&NSurface, sizeof(int),         1,        fp);

	size = NSurface * sizeof(surface_t);
	Surface = (surface_t *)malloc(size);
	num += fread(Surface,   sizeof(surface_t),   NSurface, fp);

	size = NFreq2 * sizeof(d_complex_t *);
	SurfaceEx = (d_complex_t **)malloc(size);
	SurfaceEy = (d_complex_t **)malloc(size);
	SurfaceEz = (d_complex_t **)malloc(size);
	SurfaceHx = (d_complex_t **)malloc(size);
	SurfaceHy = (d_complex_t **)malloc(size);
	SurfaceHz = (d_complex_t **)malloc(size);
	for (int ifreq = 0; ifreq < NFreq2; ifreq++) {
		size = NSurface * sizeof(d_complex_t);
		SurfaceEx[ifreq] = (d_complex_t *)malloc(size);
		SurfaceEy[ifreq] = (d_complex_t *)malloc(size);
		SurfaceEz[ifreq] = (d_complex_t *)malloc(size);
		SurfaceHx[ifreq] = (d_complex_t *)malloc(size);
		SurfaceHy[ifreq] = (d_complex_t *)malloc(size);
		SurfaceHz[ifreq] = (d_complex_t *)malloc(size);

		num += fread(SurfaceEx[ifreq], sizeof(d_complex_t), NSurface, fp);
		num += fread(SurfaceEy[ifreq], sizeof(d_complex_t), NSurface, fp);
		num += fread(SurfaceEz[ifreq], sizeof(d_complex_t), NSurface, fp);
		num += fread(SurfaceHx[ifreq], sizeof(d_complex_t), NSurface, fp);
		num += fread(SurfaceHy[ifreq], sizeof(d_complex_t), NSurface, fp);
		num += fread(SurfaceHz[ifreq], sizeof(d_complex_t), NSurface, fp);
	}

	size = NN * NFreq2 * sizeof(float *);
	cEx_r = (float *)malloc(size);
	cEx_i = (float *)malloc(size);
	cEy_r = (float *)malloc(size);
	cEy_i = (float *)malloc(size);
	cEz_r = (float *)malloc(size);
	cEz_i = (float *)malloc(size);
	cHx_r = (float *)malloc(size);
	cHx_i = (float *)malloc(size);
	cHy_r = (float *)malloc(size);
	cHy_i = (float *)malloc(size);
	cHz_r = (float *)malloc(size);
	cHz_i = (float *)malloc(size);
	for (int ifreq = 0; ifreq < NFreq2; ifreq++) {
		int64_t n0 = ifreq * NN;
		num += fread(&cEx_r[n0], sizeof(float), NN, fp);
		num += fread(&cEx_i[n0], sizeof(float), NN, fp);
		num += fread(&cEy_r[n0], sizeof(float), NN, fp);
		num += fread(&cEy_i[n0], sizeof(float), NN, fp);
		num += fread(&cEz_r[n0], sizeof(float), NN, fp);
		num += fread(&cEz_i[n0], sizeof(float), NN, fp);
		num += fread(&cHx_r[n0], sizeof(float), NN, fp);
		num += fread(&cHx_i[n0], sizeof(float), NN, fp);
		num += fread(&cHy_r[n0], sizeof(float), NN, fp);
		num += fread(&cHy_i[n0], sizeof(float), NN, fp);
		num += fread(&cHz_r[n0], sizeof(float), NN, fp);
		num += fread(&cHz_i[n0], sizeof(float), NN, fp);
	}

	size_t num0;
	size = fread(&num0, sizeof(size_t), 1, fp);
	size = size;  // suppress gcc warning
	//printf("%zd %zd\n", num, num0);

	if (num != num0) {
		fprintf(stderr, "*** invalid file length : (%zd, %zd)\n", num0, num);
	}
}
