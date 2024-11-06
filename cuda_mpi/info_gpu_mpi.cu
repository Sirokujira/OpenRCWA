/*
info_gpu_mpi.cu (CUDA + MPI) (OpenRCWA/OpenTHFD/OpenSTF)

check GPU, set device, show info
*/

#include <stdio.h>
#include <string.h>

extern int rank2device(int, int, const int []);
extern int check_gpu(int, char []);

extern "C" {
extern void comm_check(int, int, int);
extern void comm_string(const char *, char *);
}

void info_gpu_mpi(FILE *fp, int nhost, const int ndevice[], int gpu, int um, int commsize, int commrank, int prompt)
{
	if (gpu) {
		char msg[BUFSIZ], str[BUFSIZ];
		char *lstr = (char *)malloc(commsize * BUFSIZ * sizeof(char));

		int device = rank2device(commrank, nhost, ndevice);

		int ierr = check_gpu(device, msg);
		sprintf(str, "  GPU-%d: %s, U.M.%s, device=%d", commrank, msg, (um ? "ON" : "OFF"), device);
		comm_string(str, lstr);
		if (commrank == 0) {
			fprintf(fp,     "%s\n", lstr);
			fprintf(stdout, "%s\n", lstr);
			fflush(fp);
			fflush(stdout);
		}
		comm_check(ierr, 1, prompt);

		free(lstr);
	}
}
