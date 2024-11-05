/*
writeout.c

write orcwa.out
*/

#include "orcwa.h"

void writeout(FILE *fp)
{
	size_t num = 0;

	num += fwrite(Title,           sizeof(char),        256,                           fp);
	num += fwrite(&Nx,             sizeof(int),         1,                             fp);
	num += fwrite(&Ny,             sizeof(int),         1,                             fp);
	num += fwrite(&Nz,             sizeof(int),         1,                             fp);
	num += fwrite(&Ni,             sizeof(int),         1,                             fp);
	num += fwrite(&Nj,             sizeof(int),         1,                             fp);
	num += fwrite(&Nk,             sizeof(int),         1,                             fp);
	num += fwrite(&N0,             sizeof(int),         1,                             fp);
	num += fwrite(&NN,             sizeof(int64_t),     1,                             fp);
	num += fwrite(&NFreq1,         sizeof(int),         1,                             fp);
	num += fwrite(&NFreq2,         sizeof(int),         1,                             fp);
	num += fwrite(&NFeed,          sizeof(int),         1,                             fp);
	num += fwrite(&NPoint,         sizeof(int),         1,                             fp);
	num += fwrite(&Niter,          sizeof(int),         1,                             fp);
	num += fwrite(&Ntime,          sizeof(int),         1,                             fp);
	num += fwrite(&Solver.maxiter, sizeof(int),         1,                             fp);
	num += fwrite(&Solver.nout,    sizeof(int),         1,                             fp);
	num += fwrite(&Dt,             sizeof(double),      1,                             fp);
	num += fwrite(&NGline,         sizeof(int),         1,                             fp);
	num += fwrite(&IPlanewave,     sizeof(int),         1,                             fp);
	num += fwrite(&Planewave,      sizeof(planewave_t), 1,                             fp);

	num += fwrite(Xn,              sizeof(double),      Nx + 1,                        fp);
	num += fwrite(Yn,              sizeof(double),      Ny + 1,                        fp);
	num += fwrite(Zn,              sizeof(double),      Nz + 1,                        fp);
	num += fwrite(Xc,              sizeof(double),      Nx,                            fp);
	num += fwrite(Yc,              sizeof(double),      Ny,                            fp);
	num += fwrite(Zc,              sizeof(double),      Nz,                            fp);
	num += fwrite(Eiter,           sizeof(double),      Niter,                         fp);
	num += fwrite(Hiter,           sizeof(double),      Niter,                         fp);
	num += fwrite(VFeed,           sizeof(double),      NFeed  * (Solver.maxiter + 1), fp);
	num += fwrite(IFeed,           sizeof(double),      NFeed  * (Solver.maxiter + 1), fp);
	num += fwrite(VPoint,          sizeof(double),      NPoint * (Solver.maxiter + 1), fp);
	num += fwrite(Freq1,           sizeof(double),      NFreq1,                        fp);
	num += fwrite(Freq2,           sizeof(double),      NFreq2,                        fp);
	num += fwrite(Feed,            sizeof(feed_t),      NFeed,                         fp);
	num += fwrite(Zin,             sizeof(d_complex_t), NFeed * NFreq1,                fp);
	num += fwrite(Ref,             sizeof(double),      NFeed * NFreq1,                fp);
	num += fwrite(Pin[0],          sizeof(double),      NFeed * NFreq2,                fp);
	num += fwrite(Pin[1],          sizeof(double),      NFeed * NFreq2,                fp);
	num += fwrite(Spara,           sizeof(d_complex_t), NPoint * NFreq1,               fp);
	num += fwrite(Gline,           sizeof(double),      NGline * 2 * 3,                fp);

	num += fwrite(&NSurface,       sizeof(int),         1,        fp);
	num += fwrite(Surface,         sizeof(surface_t),   NSurface, fp);
	for (int ifreq = 0; ifreq < NFreq2; ifreq++) {
		num += fwrite(SurfaceEx[ifreq], sizeof(d_complex_t), NSurface, fp);
		num += fwrite(SurfaceEy[ifreq], sizeof(d_complex_t), NSurface, fp);
		num += fwrite(SurfaceEz[ifreq], sizeof(d_complex_t), NSurface, fp);
		num += fwrite(SurfaceHx[ifreq], sizeof(d_complex_t), NSurface, fp);
		num += fwrite(SurfaceHy[ifreq], sizeof(d_complex_t), NSurface, fp);
		num += fwrite(SurfaceHz[ifreq], sizeof(d_complex_t), NSurface, fp);
	}

	for (int ifreq = 0; ifreq < NFreq2; ifreq++) {
		int64_t n0 = ifreq * NN;
		num += fwrite(&cEx_r[n0], sizeof(float), NN, fp);
		num += fwrite(&cEx_i[n0], sizeof(float), NN, fp);
		num += fwrite(&cEy_r[n0], sizeof(float), NN, fp);
		num += fwrite(&cEy_i[n0], sizeof(float), NN, fp);
		num += fwrite(&cEz_r[n0], sizeof(float), NN, fp);
		num += fwrite(&cEz_i[n0], sizeof(float), NN, fp);
		num += fwrite(&cHx_r[n0], sizeof(float), NN, fp);
		num += fwrite(&cHx_i[n0], sizeof(float), NN, fp);
		num += fwrite(&cHy_r[n0], sizeof(float), NN, fp);
		num += fwrite(&cHy_i[n0], sizeof(float), NN, fp);
		num += fwrite(&cHz_r[n0], sizeof(float), NN, fp);
		num += fwrite(&cHz_i[n0], sizeof(float), NN, fp);
	}

	fwrite(&num, sizeof(size_t), 1, fp);
}
