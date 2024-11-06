/*
monitor.c
*/

#include "orcwa.h"

static void memory_size(int *, int *);
static int output_size(void);


// title or message
static void monitor1_(FILE *fp, const char msg[])
{
	fprintf(fp, "%s\n", msg);
	fflush(fp);
}


// condition
static void monitor2_(FILE *fp, int gpu, int commsize)
{
	int cpu_mem, gpu_mem;
	memory_size(&cpu_mem, &gpu_mem);

	time_t now;
	time(&now);

	fprintf(fp, "%s", ctime(&now));
	fprintf(fp, "Title = %s\n", Title);
	fprintf(fp, "Source = %s\n", (NFeed ? "feed" : "plane wave"));
	fprintf(fp, "Cells = %d x %d x %d = %lld\n", Nx, Ny, Nz, ((long long int)Nx * Ny * Nz));
	fprintf(fp, "No. of Materials   = %zd\n", NMaterial);
	fprintf(fp, "No. of Geometries  = %zd\n", NGeometry);
	if (NFeed) fprintf(fp, "No. of Feeds       = %d\n", NFeed);
	fprintf(fp, "No. of Points      = %d\n", NPoint);
	fprintf(fp, "No. of Freq.s (1)  = %d\n", NFreq1);
	fprintf(fp, "No. of Freq.s (2)  = %d\n", NFreq2);
	fprintf(fp, "CPU Memory size    = %d [MB]\n", cpu_mem);
	if (gpu) fprintf(fp, "GPU Memory %s = %d [MB]\n", (commsize > 1 ? "/ proc." : "size   "), gpu_mem);
	fprintf(fp, "Output filesize    = %d [MB]\n", output_size());
	if (iABC == 0) fprintf(fp, "ABC = Mur-1st\n");
	if (iABC == 1) fprintf(fp, "ABC = PML (L=%d, M=%.2f, R0=%.2e)\n", cPML.l, cPML.m, cPML.r0);
	if (PBCx || PBCy || PBCz) fprintf(fp, "PBC : %s%s%s\n", (PBCx ? "X" : ""), (PBCy ? "Y" : ""), (PBCz ? "Z" : ""));
	fprintf(fp, "Dt[sec] = %.4e, Tw[sec] = %.4e, Tw/Dt = %.3f\n", Dt, Tw, Tw / Dt);
	fprintf(fp, "Iterations = %d, Convergence = %.3e\n", Solver.maxiter, Solver.converg);
	fprintf(fp, "=== iteration start ===\n");
	fprintf(fp, "   step   <E>      <H>\n");

	fflush(fp);
}


// output files
static void monitor3_(FILE *fp, const char fn_log[], const char fn_out[])
{
	fprintf(fp, "=== output files ===\n");
	fprintf(fp, "%s, %s\n", fn_log, fn_out);
	fflush(fp);
}


// cpu time
static void monitor4_(FILE *fp, const double cpu[])
{
	time_t now;
	time(&now);

	fprintf(fp, "%s\n", "=== cpu time [sec] ===");
	fprintf(fp, "  part-1 : %11.3f\n", cpu[1] - cpu[0]);
	fprintf(fp, "  part-2 : %11.3f\n", cpu[2] - cpu[1]);
	fprintf(fp, "  part-3 : %11.3f\n", cpu[3] - cpu[2]);
	fprintf(fp, "  part-4 : %11.3f\n", cpu[4] - cpu[3]);
	fprintf(fp, "  %s\n", "--------------------");
	fprintf(fp, "  total  : %11.3f\n", cpu[4] - cpu[0]);
	fprintf(fp, "%s\n", "=== normal end ===");
	fprintf(fp, "%s", ctime(&now));

	fflush(fp);
}


void monitor1(FILE *fp, const char msg[])
{
	monitor1_(fp,     msg);
	monitor1_(stdout, msg);
}


void monitor2(FILE *fp, int gpu, int commsize)
{
	monitor2_(fp,     gpu, commsize);
	monitor2_(stdout, gpu, commsize);
}


void monitor3(FILE *fp, const char fn_log[], const char fn_out[])
{
	monitor3_(fp,     fn_log, fn_out);
	monitor3_(stdout, fn_log, fn_out);
}


void monitor4(FILE *fp, const double cpu[])
{
	monitor4_(fp,     cpu);
	monitor4_(stdout, cpu);
}


// memory size [MB]
static void memory_size(int *cpu, int *gpu)
{
	int64_t mem = 0;

	mem += 6 * NN * sizeof(real_t);        // Ex, Ey, Ez, Hx, Hy, Hz
	mem += 6 * NN * sizeof(id_t);          // iEx, iEy, iEz, iHx, iHy, iHz

	if (VECTOR) {
		mem += 12 * NN * sizeof(real_t);   // K1Ex, K2Ex, ... K1Hx, K2Hx, ...
	}

	// PML
	if (iABC == 1) {
		mem += (numPmlEx + numPmlEy + numPmlEz
		      + numPmlHx + numPmlHy + numPmlHz)
		     * (sizeof(pml_t) + 2 * sizeof(real_t));
	}

	// cEx, cEy, cEz, cHx, cHy, cHz
	int64_t cpu_mem = mem + NFreq2 * Nx * Ny * Nz * 12 * sizeof(float);
	int64_t gpu_mem = mem + NFreq2 * NN * 12 * sizeof(float);

	*cpu = (int)(cpu_mem / 1024 / 1024) + 1;
	*gpu = (int)(gpu_mem / 1024 / 1024) + 1;
}


// output filesize [MB]
static int output_size(void)
{
	return (int)(NFreq2 * Nx * Ny * Nz * 12 * sizeof(float) / 1024 / 1024) + 1;
}
