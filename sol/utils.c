/*
utils.c

utilities
*/

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>

typedef struct {double r, i;} d_complex_t;
typedef struct {float  r, i;} f_complex_t;

#include "complex.h"

// tokenize a string
int tokenize(char *str, const char *tokensep, char *token[], int maxtoken)
{
	if ((str == NULL) || !maxtoken) return 0;

	char *thistoken = strtok(str, tokensep);

	int   count;
	for (count = 0; (count < maxtoken) && (thistoken != NULL); ) {
		token[count++] = thistoken;
		thistoken = strtok(NULL, tokensep);
	}

	token[count] = NULL;

	return count;
}


// nearest id : p[n1]...p[n2]
int nearest(double x, int n1, int n2, const double *p)
{
	int    imin = n1;
	double dmin = fabs(x - p[n1]);
	for (int i = n1 + 1; i <= n2; i++) {
		double d = fabs(x - p[i]);
		if (d < dmin) {
			dmin = d;
			imin = i;
		}
	}

	return imin;
}


// p[i1] <= p[n1] <= p1 <= p2 <= p[n2] <= p[i2]
// output : n1, n2
void getspan(const double p[], int n, int i1, int i2, double p1, double p2, int *n1, int *n2, double eps)
{
	if (i1 < 0) i1 = 0;
	if (i2 < 0) i2 = 0;
	if (i1 > n - 1) i1 = n - 1;
	if (i2 > n - 1) i2 = n - 1;

	// k1 <= k2
	int k1 = (i1 < i2) ? i1 : i2;
	int k2 = (i1 > i2) ? i1 : i2;

	// q1 <= q2
	double q1 = (p1 < p2) ? p1 : p2;
	double q2 = (p1 > p2) ? p1 : p2;

	*n1 = k1;
	*n2 = k2;

	if ((*n2 - *n1) < 1) return;

	if      ((q1 < p[k1] - eps) && (q2 < p[k1] - eps)) {
		// p1, p2 < p[k1] -> n1 > n2
		*n1 = k1 - 1;
		*n2 = k1 - 2;
	}
	else if ((q1 > p[k2] + eps) && (q2 > p[k2] + eps)) {
		// p1, p2 > p[k2] -> n1 > n2
		*n1 = k2 + 2;
		*n2 = k2 + 1;
	}
	else {
		// p[n1] <= p1
		for (int k = k1; k < k2; k++) {
			if ((q1 > p[k] - eps) && (q1 < p[k + 1] - eps)) {
				*n1 = k;
				//printf("A %d\n", *n1);
				break;
			}
		}
		if (fabs(q1 - p[k2]) < eps) {
			*n1 = k2;
			//printf("A2 %d\n", *n1);
		}
		// p2 <= p[n2]
		for (int k = k2; k > k1; k--) {
			if ((q2 < p[k] + eps) && (q2 > p[k - 1] + eps)) {
				*n2 = k;
				//printf("B %d\n", *n2);
				break;
			}
		}
		if (fabs(q2 - p[k1]) < eps) {
			*n2 = k1;
			//printf("B2 %d\n", *n2);
		}

		//assert((i1 <= *n1) && (*n1 <= *n2) && (*n2 <= i2));
	}
}


// DFT
d_complex_t calcdft(int ntime, const double f[], double freq, double dt, double shift)
{
	const double pi = 4 * atan(1);
	const double omega = 2 * pi * freq;

	double sum_r = 0;
	double sum_i = 0;
	for (int n = 0; n < ntime; n++) {
		const double ot = omega * (n + shift) * dt;
		sum_r += cos(ot) * f[n];
		sum_i -= sin(ot) * f[n];
	}

	return d_complex(sum_r, sum_i);
}
