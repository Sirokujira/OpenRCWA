/*
OpenRCWA Version 4.2.1

post process
*/

#define MAIN
#include "orcwa.h"
#undef MAIN

#include "orcwa_prototype.h"

static void args(int, char *[], int *, int *, char [], char []);

int main(int argc, char *argv[])
{
	const char errfmt[] = "*** file %s open error.\n";
	FILE *fp_in = NULL, *fp_out = NULL;

	// arguments
	HTML = 0;
	int nthread = 1;
	int prompt = 0;
	char fn_in[BUFSIZ] = "";
	char fn_out[BUFSIZ] = "orcwa.out";
	args(argc, argv, &nthread, &prompt, fn_in, fn_out);

	// set number of threads
#ifdef _OPENMP
	omp_set_num_threads(nthread);
#endif

	// input data
	if ((fp_in = fopen(fn_in, "r")) == NULL) {
		printf(errfmt, fn_in);
		if (prompt) getchar();
		exit(1);
	}
	if (post_data(fp_in)) {
		fclose(fp_in);
		if (prompt) getchar();
		exit(1);
	}

	// read orcwa.out
	if ((fp_out = fopen(fn_out, "rb")) == NULL) {
		printf(errfmt, fn_out);
		if (prompt) getchar();
		exit(1);
	}
	readout(fp_out);
	fclose(fp_out);

	// post process
	post();

	// prompt
	if (prompt) getchar();

	return 0;
}


static void args(int argc, char *argv[],
	int *nthread, int *prompt, char fn_in[], char fn_out[])
{
	const char usage[] = "Usage : orcwa_post [-n <thread>] [-html] [-out <outfile>] <datafile>";

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
		else if (!strcmp(*argv, "-html")) {
			HTML = 1;
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
