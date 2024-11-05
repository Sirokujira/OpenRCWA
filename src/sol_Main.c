/*
OpenRCWA Version 4.2.1 (CPU + OpenMP)

solver
*/

#define MAIN
#include "orcwa.h"
#undef MAIN

#include "orcwa_prototype.h"

static void args(int, char *[], int *, int *, char [], char []);
static void error_check(int, int);

int main(int argc, char *argv[])
{
	const char prog[] = "(CPU+OpenMP)";
	const char errfmt[] = "*** file %s open error.\n";
	char str[BUFSIZ];
	int ierr = 0;
	double cpu[] = {0, 0, 0, 0, 0};
	FILE *fp_in = NULL, *fp_out = NULL, *fp_log = NULL;

	// process (MPI)
	Npx = Npy = Npz = 1;
	const int io = 1;

	// arguments
	VECTOR = 0;
	int nthread = 1;
	int prompt = 0;
	char fn_in[BUFSIZ] = "";
	char fn_out[BUFSIZ] = "orcwa.out";
	args(argc, argv, &nthread, &prompt, fn_in, fn_out);
	//printf("%d %d %s %s\n", nthread, prompt, fn_in, fn_out);

	// set number of threads
#ifdef _OPENMP
	omp_set_num_threads(nthread);
#endif

	// cpu time
	cpu[0] = cputime();

	// input data
	if (io) {
		if ((fp_in = fopen(fn_in, "r")) == NULL) {
			printf(errfmt, fn_in);
			ierr = 1;
		}
		if (!ierr) {
			ierr = input_data(fp_in);
			fclose(fp_in);
		}
	}
	error_check(ierr, prompt);

	// open log file
	if (io) {
		if ((fp_log = fopen(FN_log, "w")) == NULL) {
			printf(errfmt, FN_log);
			ierr = 1;
		}
	}
	error_check(ierr, prompt);

	// monitor
	if (io) {
		// logo
		sprintf(str, "<<< %s %s Ver.%d.%d.%d >>>", PROGRAM, prog, VERSION_MAJOR, VERSION_MINOR, VERSION_BUILD);
		monitor1(fp_log, str);
		// thread and process
		sprintf(str, "CPU, thread=%d, process=1, vector=%s", nthread, (VECTOR ? "on" : "off"));
		monitor1(fp_log, str);
	}

	// plot geometry 3d and exit
	if (io && Plot3dGeom) {
		plot3dGeom();
		ierr = 1;
	}
	error_check(ierr, prompt);

	// setup
	setupSize(1, 1, 1, 0);
	setupABCsize();
	memalloc1();
	memalloc2();
	memalloc3();
	setup();

	// monitor
	if (io) {
		monitor2(fp_log, 0, 1);
	}

	// cpu time
	cpu[1] = cputime();
	double tdft = 0;

	// solve
	solve(io, &tdft, fp_log);

	// cpu time
	cpu[3] = cputime();
	cpu[2] = cpu[3] - tdft;

	// output
	if (io) {
		// calculation and output
		outputChars(fp_log);

		// output filenames
		monitor3(fp_log, FN_log, fn_out);

		// write orcwa.out
		if ((fp_out = fopen(fn_out, "wb")) == NULL) {
			printf(errfmt, fn_out);
			ierr = 1;
		}
		if (!ierr) {
			writeout(fp_out);
			fclose(fp_out);
		}
	}
	error_check(ierr, prompt);

	// free
	memfree1();
	memfree3();

	// cpu time
	cpu[4] = cputime();

	if (io) {
		// cpu time
		monitor4(fp_log, cpu);

		// close log file
		if (fp_log != NULL) {
			fclose(fp_log);
		}
	}

	// prompt
	if (io && prompt) getchar();

	return 0;
}


static void args(int argc, char *argv[],
	int *nthread, int *prompt, char fn_in[], char fn_out[])
{
	const char usage[] = "Usage : orcwa [-n <thread>] [-no-vector|-vector] [-out <outfile>] <datafile>";

	if (argc < 2) {
		printf("%s\n", usage);
		exit(0);
	}

	while (--argc) {
		++argv;
		if (!strcmp(*argv, "-n")) {
			if (--argc) {
				*nthread = atoi(*++argv);
				if (*nthread < 1) *nthread = 1;
			}
			else {
				break;
			}
		}
		else if (!strcmp(*argv, "-no-vector")) {
			VECTOR = 0;
		}
		else if (!strcmp(*argv, "-vector")) {
			VECTOR = 1;
		}
		else if (!strcmp(*argv, "-prompt")) {
			*prompt = 1;
		}
		else if (!strcmp(*argv, "-out")) {
			argc--;
			strcpy(fn_out, *++argv);
		}
		else if (!strcmp(*argv, "--help")) {
			printf("%s\n", usage);
			exit(0);
		}
		else if (!strcmp(*argv, "--version")) {
			printf("%s Ver.%d.%d.%d\n", PROGRAM, VERSION_MAJOR, VERSION_MINOR, VERSION_BUILD);
			exit(0);
		}
		else {
			strcpy(fn_in, *argv);
		}
	}
}


// error check
static void error_check(int ierr, int prompt)
{
	if (ierr) {
		if (prompt) {
			fflush(stdout);
			getchar();
		}
		exit(1);
	}
}
