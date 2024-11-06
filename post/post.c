/*
post.c

post process
*/

#include "orcwa.h"
#include "orcwa_prototype.h"
#include "ev.h"

static void setup_near1d(void);
static void setup_near2d(void);


void post(void)
{
	ev2d_init(Width2d, Height2d);

	ev3d_init();

	if (Piter) {
		plot2dIter();
	}

	if (Pfeed) {
		plot2dFeed();
	}

	if (Ppoint) {
		plot2dPoint();
	}

	if (NFreq1) {
		plot2dFreq();
	}

	if (NFreq2) {
		if (IFar0d) {
			outputFar0d();
		}

		if (NFar1d) {
			outputFar1d();
		}

		if (NFar2d) {
			outputFar2d();
		}

		if (NNear1d) {
			setup_near1d();
			calcNear1d();
			outputNear1d();
		}

		if (NNear2d) {
			setup_near2d();
			calcNear2d();
			outputNear2d();
		}
	}

	ev2d_file(!HTML, (!HTML ? FN_ev2d_1 : FN_ev2d_0));
	ev2d_output();

	ev3d_file(!HTML, (!HTML ? FN_ev3d_1 : FN_ev3d_0), 0);
	ev3d_output();
}


// setup near1d
static void setup_near1d(void)
{
	for (int n = 0; n < NNear1d; n++) {
		if      (Near1d[n].dir == 'X') {
			Near1d[n].id1 = nearest(Near1d[n].pos1, 0, Ny, Yn);
			Near1d[n].id2 = nearest(Near1d[n].pos2, 0, Nz, Zn);
		}
		else if (Near1d[n].dir == 'Y') {
			Near1d[n].id1 = nearest(Near1d[n].pos1, 0, Nz, Zn);
			Near1d[n].id2 = nearest(Near1d[n].pos2, 0, Nx, Xn);
		}
		else if (Near1d[n].dir == 'Z') {
			Near1d[n].id1 = nearest(Near1d[n].pos1, 0, Nx, Xn);
			Near1d[n].id2 = nearest(Near1d[n].pos2, 0, Ny, Yn);
		}
		//printf("near1d : %d %c %f %d %f %d\n", n, Near1d[n].dir, Near1d[n].pos1, Near1d[n].id1, Near1d[n].pos2, Near1d[n].id2);
	}
}


// setup near2d
static void setup_near2d(void)
{
	for (int n = 0; n < NNear2d; n++) {
		if      (Near2d[n].dir == 'X') {
			Near2d[n].id0 = nearest(Near2d[n].pos0, 0, Nx, Xn);
		}
		else if (Near2d[n].dir == 'Y') {
			Near2d[n].id0 = nearest(Near2d[n].pos0, 0, Ny, Yn);
		}
		else if (Near2d[n].dir == 'Z') {
			Near2d[n].id0 = nearest(Near2d[n].pos0, 0, Nz, Zn);
		}
		//printf("near2d : %d %s %c %f %d\n", n, Near2d[n].cmp, Near2d[n].dir, Near2d[n].pos0, Near2d[n].id0);
	}
}
