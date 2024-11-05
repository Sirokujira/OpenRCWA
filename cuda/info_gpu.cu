/*
info_gpu.cu (CUDA) (OpenRCWA/OpenTHFD/OpenSTF)

check GPU, set device, show info
*/

#include <stdio.h>
#include <string.h>

extern int check_gpu(int, char []);

void info_gpu(FILE *fp, int device, int gpu, int um)
{
	if (gpu) {
		char str[BUFSIZ], msg[BUFSIZ];
		if (check_gpu(device, msg)) {
			printf("%s\n", msg);
			getchar();
			exit(1);
		}
		sprintf(str, "GPU : %s, U.M.%s, device=%d", msg, (um ? "ON" : "OFF"), device);
		fprintf(fp,     "%s\n", str);
		fprintf(stdout, "%s\n", str);
		fflush(fp);
		fflush(stdout);
	}
}
