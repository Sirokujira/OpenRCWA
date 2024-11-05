/*
solve.cu (CUDA)
*/

#include "orcwa.h"
#include "orcwa_cuda.h"
#include "orcwa_prototype.h"

static void copy_to_host();

void solve(int io, double *tdft, FILE *fp)
{
	double fmax[] = {0, 0};
	char   str[BUFSIZ];
	int    converged = 0;

	// setup host memory
	setup_host();

	// setup (GPU)
	if (GPU) {
		setup_gpu();
	}

	// initial field
	initfield();

	// time step iteration
	int itime;
	double t = 0;
	for (itime = 0; itime <= Solver.maxiter; itime++) {

		// update H
		t += 0.5 * Dt;
		updateHx(t);
		updateHy(t);
		updateHz(t);

		// ABC H
		if      (iABC == 0) {
			murH(numMurHx, (GPU ? d_fMurHx : fMurHx), Hx);
			murH(numMurHy, (GPU ? d_fMurHy : fMurHy), Hy);
			murH(numMurHz, (GPU ? d_fMurHz : fMurHz), Hz);
		}
		else if (iABC == 1) {
			pmlHx();
			pmlHy();
			pmlHz();
		}

		// PBC H
		if (PBCx) {
			pbcx();
		}
		if (PBCy) {
			pbcy();
		}
		if (PBCz) {
			pbcz();
		}

		// update E
		t += 0.5 * Dt;
		updateEx(t);
		updateEy(t);
		updateEz(t);

		// dispersion E
		if (numDispersionEx) {
			dispersionEx(t);
		}
		if (numDispersionEy) {
			dispersionEy(t);
		}
		if (numDispersionEz) {
			dispersionEz(t);
		}

		// ABC E
		if      (iABC == 1) {
			pmlEx();
			pmlEy();
			pmlEz();
		}

		// feed
		if (NFeed) {
			efeed(itime);
		}

		// inductor
		if (NInductor) {
			eload();
		}

		// point
		if (NPoint) {
			vpoint(itime);
		}

		// DFT
		if (GPU) cudaDeviceSynchronize();
		const double t0 = cputime();
		dftNear3d(itime);
		if (GPU) cudaDeviceSynchronize();
		*tdft += cputime() - t0;

		// average and convergence
		if ((itime % Solver.nout == 0) || (itime == Solver.maxiter)) {
			// average
			double fsum[2];
			average(fsum);

			// average (plot)
			Eiter[Niter] = fsum[0];
			Hiter[Niter] = fsum[1];
			Niter++;

			// monitor
			if (io) {
				sprintf(str, "%7d %.6f %.6f", itime, fsum[0], fsum[1]);
				fprintf(fp,     "%s\n", str);
				fprintf(stdout, "%s\n", str);
				fflush(fp);
				fflush(stdout);
			}

			// check convergence
			fmax[0] = MAX(fmax[0], fsum[0]);
			fmax[1] = MAX(fmax[1], fsum[1]);
			if ((fsum[0] < fmax[0] * Solver.converg) &&
			    (fsum[1] < fmax[1] * Solver.converg)) {
				converged = 1;
				break;
			}
		}
	}

	// result
	if (io) {
		sprintf(str, "    --- %s ---", (converged ? "converged" : "max steps"));
		fprintf(fp,     "%s\n", str);
		fprintf(stdout, "%s\n", str);
		fflush(fp);
		fflush(stdout);
	}

	// time steps
	Ntime = itime + converged;

	// copy point from device to host
	if (GPU) {
		copy_to_host();
	}

	// free
	memfree2_gpu();

	// copy near3d from device to host
	memcopy3_gpu();

	// free
	memfree3_gpu();
}


// copy from device to host
static void copy_to_host()
{
	if (NFeed) {
		cuda_memcpy(GPU, VFeed, d_VFeed, Feed_size, cudaMemcpyDeviceToHost);
		cuda_memcpy(GPU, IFeed, d_IFeed, Feed_size, cudaMemcpyDeviceToHost);
	}

	if (NPoint) {
		cuda_memcpy(GPU, VPoint, d_VPoint, Point_size, cudaMemcpyDeviceToHost);
	}
}
