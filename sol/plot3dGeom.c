/*
plot3dGeom.c

plot geometry 3D
*/

#include "orcwa.h"
#include "ev.h"
#include "orcwa_prototype.h"

// geometry
static void plot3dGeom_g(void)
{
	const unsigned char rgb[][3] = {
		{  0,   0,   0},	// PEC
		{255,   0, 255}		// dielectrics
	};

	for (int n = 0; n < NGline; n++) {
		int m = (MGline[n] == PEC) ? 0 : 1;
		ev3d_setColor(rgb[m][0], rgb[m][1], rgb[m][2]);
		ev3d_drawLine(Gline[n][0][0], Gline[n][0][1], Gline[n][0][2],
		              Gline[n][1][0], Gline[n][1][1], Gline[n][1][2]);
	}
}


// mesh
static void plot3dGeom_m(void)
{
	if ((Nx <= 0) || (Ny <= 0) || (Nz <= 0)) return;

	// gray
	ev3d_setColor(200, 200, 200);

	double x1 = Xn[0];
	double x2 = Xn[Nx];
	double y1 = Yn[0];
	double y2 = Yn[Ny];
	double z1 = Zn[0];
	double z2 = Zn[Nz];

	// X constant
	for (int i = 0; i <= Nx; i++) {
		double x = Xn[i];
		ev3d_drawLine(x, y1, z1, x, y2, z1);
		ev3d_drawLine(x, y1, z1, x, y1, z2);
	}

	// Y constant
	for (int j = 0; j <= Ny; j++) {
		double y = Yn[j];
		ev3d_drawLine(x1, y, z1, x1, y, z2);
		ev3d_drawLine(x1, y, z1, x2, y, z1);
	}

	// Z constant
	for (int k = 0; k <= Nz; k++) {
		double z = Zn[k];
		ev3d_drawLine(x1, y1, z, x2, y1, z);
		ev3d_drawLine(x1, y1, z, x1, y2, z);
	}
}


// feed
static void plot3dGeom_f(void)
{
	if (NFeed <= 0) return;

	// red
	ev3d_setColor(255, 0, 0);

	for (int n = 0; n < NFeed; n++) {
		int i = Feed[n].i;
		int j = Feed[n].j;
		int k = Feed[n].k;
		if      (Feed[n].dir == 'X') {
			ev3d_drawLine(Xn[i], Yn[j], Zn[k], Xn[i + 1], Yn[j], Zn[k]);
		}
		else if (Feed[n].dir == 'Y') {
			ev3d_drawLine(Xn[i], Yn[j], Zn[k], Xn[i], Yn[j + 1], Zn[k]);
		}
		else if (Feed[n].dir == 'Z') {
			ev3d_drawLine(Xn[i], Yn[j], Zn[k], Xn[i], Yn[j], Zn[k + 1]);
		}
	}
}


// load (inductor)
static void plot3dGeom_l(void)
{
	if (NInductor <= 0) return;

	// orange
	ev3d_setColor(255, 165, 0);

	for (int n = 0; n < NInductor; n++) {
		char dir = Inductor[n].dir;
		int i    = Inductor[n].i;
		int j    = Inductor[n].j;
		int k    = Inductor[n].k;
		if      (dir == 'X') {
			ev3d_drawLine(Xn[i], Yn[j], Zn[k], Xn[i + 1], Yn[j], Zn[k]);
		}
		else if (dir == 'Y') {
			ev3d_drawLine(Xn[i], Yn[j], Zn[k], Xn[i], Yn[j + 1], Zn[k]);
		}
		else if (dir == 'Z') {
			ev3d_drawLine(Xn[i], Yn[j], Zn[k], Xn[i], Yn[j], Zn[k + 1]);
		}
	}
}


// point
static void plot3dGeom_p(void)
{
	if (NPoint <= 0) return;

	// green
	ev3d_setColor(0, 255, 0);

	for (int n = 0; n < NPoint + 2; n++) {
		int i = Point[n].i;
		int j = Point[n].j;
		int k = Point[n].k;
		if      (Point[n].dir == 'X') {
			ev3d_drawLine(Xn[i], Yn[j], Zn[k], Xn[i + 1], Yn[j], Zn[k]);
		}
		else if (Point[n].dir == 'Y') {
			ev3d_drawLine(Xn[i], Yn[j], Zn[k], Xn[i], Yn[j + 1], Zn[k]);
		}
		else if (Point[n].dir == 'Z') {
			ev3d_drawLine(Xn[i], Yn[j], Zn[k], Xn[i], Yn[j], Zn[k + 1]);
		}
	}
}


void plot3dGeom(void)
{
	// initialize
	ev3d_init();

	// new page
	ev3d_newPage();

	// mesh
	plot3dGeom_m();

	// geometry
	plot3dGeom_g();

	// feed
	plot3dGeom_f();

	// load (inductor)
	plot3dGeom_l();

	// point
	plot3dGeom_p();

	// title
	char str[BUFSIZ];
	ev3d_setColor(0, 0, 0);
	ev3d_drawTitle(Title);
	sprintf(str, "No. of geometries = %zd", NGeometry);
	ev3d_drawTitle(str);
	sprintf(str, "Nx=%d Ny=%d Nz=%d", Nx, Ny, Nz);
	ev3d_drawTitle(str);

	// output HTML
	ev3d_html_size(500, 500);
	ev3d_file(0, FN_geom3d_0, 0);
	ev3d_output();

	// output ev3
	ev3d_file(1, FN_geom3d_1, 0);
	ev3d_output();

	// message
	printf("output : %s, %s\n", FN_geom3d_0, FN_geom3d_1);
}
