/*
ev3dlibs.c

ev3d utilities (3D)
*/

#include <math.h>
#include "ev.h"


// plot r=r(theta,phi)
void ev3dlib_func(int nth, double th1, double th2, int nph, double ph1, double ph2, double **rdata)
{
	const double dtor = atan(1) / 45;

	// error check
	if (nth <= 0) return;
	if (nph <= 0) return;

	// min and max
	double rmin = rdata[0][0];
	double rmax = rmin;
	for (int ith = 0; ith <= nth; ith++) {
	for (int iph = 0; iph <= nph; iph++) {
		double r = rdata[ith][iph];
		if (r < rmin) rmin = r;
		if (r > rmax) rmax = r;
	}
	}

	// offset
	double off = 0;
	if (rmin < 0) {
		off -= rmin;
		rmax -= rmin;
	}

	// division
	const double dth = (th2 - th1) / nth;
	const double dph = (ph2 - ph1) / nph;

	// plot
	for (int ith = 0; ith < nth; ith++) {
	for (int iph = 0; iph < nph; iph++) {
		const double t1 = (th1 + (ith + 0) * dth) * dtor;
		const double t2 = (th1 + (ith + 1) * dth) * dtor;
		const double p1 = (ph1 + (iph + 0) * dph) * dtor;
		const double p2 = (ph1 + (iph + 1) * dph) * dtor;
		const double st1 = sin(t1);
		const double st2 = sin(t2);
		const double sp1 = sin(p1);
		const double sp2 = sin(p2);
		const double ct1 = cos(t1);
		const double ct2 = cos(t2);
		const double cp1 = cos(p1);
		const double cp2 = cos(p2);

		const double r1 = rdata[ith + 0][iph + 0] + off;
		const double r2 = rdata[ith + 0][iph + 1] + off;
		const double r3 = rdata[ith + 1][iph + 1] + off;
		const double r4 = rdata[ith + 1][iph + 0] + off;
		const double rav = (r1 + r2 + r3 + r4) / 4;

		const double x1 = r1 * st1 * cp1;
		const double x2 = r2 * st1 * cp2;
		const double x3 = r3 * st2 * cp2;
		const double x4 = r4 * st2 * cp1;
		const double y1 = r1 * st1 * sp1;
		const double y2 = r2 * st1 * sp2;
		const double y3 = r3 * st2 * sp2;
		const double y4 = r4 * st2 * sp1;
		const double z1 = r1 * ct1;
		const double z2 = r2 * ct1;
		const double z3 = r3 * ct2;
		const double z4 = r4 * ct2;

		ev3d_setColorV(rav / rmax, 1);
		ev3d_fillQuadrangle(x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4);
	}
	}
}


// array argument
void ev3dlib_drawLineA(double g[2][3])
{
	ev3d_drawLine(g[0][0], g[0][1], g[0][2],
	              g[1][0], g[1][1], g[1][2]);
}
void ev3dlib_drawTriangleA(double g[3][3])
{
	ev3d_drawTriangle(g[0][0], g[0][1], g[0][2],
	                  g[1][0], g[1][1], g[1][2],
	                  g[2][0], g[2][1], g[2][2]);
}
void ev3dlib_fillTriangleA(double g[3][3])
{
	ev3d_fillTriangle(g[0][0], g[0][1], g[0][2],
	                  g[1][0], g[1][1], g[1][2],
	                  g[2][0], g[2][1], g[2][2]);
}
void ev3dlib_drawQuadrangleA(double g[4][3])
{
	ev3d_drawQuadrangle(g[0][0], g[0][1], g[0][2],
	                    g[1][0], g[1][1], g[1][2],
	                    g[2][0], g[2][1], g[2][2],
	                    g[3][0], g[3][1], g[3][2]);
}
void ev3dlib_fillQuadrangleA(double g[4][3])
{
	ev3d_fillQuadrangle(g[0][0], g[0][1], g[0][2],
	                    g[1][0], g[1][1], g[1][2],
	                    g[2][0], g[2][1], g[2][2],
	                    g[3][0], g[3][1], g[3][2]);
}


// X/Y/Z axis (not used)
void ev3dlib_axis(double x0, double y0, double z0, double r, int xyz)
{
	const double d = r / 15;

	// X
	ev3d_drawLine(x0 - r, y0, z0, x0 + r, y0, z0);
	if (xyz) {
		ev3d_drawLine(x0 + r + 0.2 * d, y0, z0 - d, x0 + r + 1.8 * d, y0, z0 + d);
		ev3d_drawLine(x0 + r + 1.8 * d, y0, z0 - d, x0 + r + 0.2 * d, y0, z0 + d);
	}

	// Y
	ev3d_drawLine(x0, y0 - r, z0, x0, y0 + r, z0);
	if (xyz) {
		ev3d_drawLine(x0, y0 + r + 0.4 * d, z0 + d, x0, y0 + r + d, z0);
		ev3d_drawLine(x0, y0 + r + 1.6 * d, z0 + d, x0, y0 + r + d, z0);
		ev3d_drawLine(x0, y0 + r + 1.0 * d, z0 - d, x0, y0 + r + d, z0);
	}

	// Z
	ev3d_drawLine(x0, y0, z0 - r, x0, y0, z0 + r);
	if (xyz) {
		ev3d_drawLine(x0 - 0.8 * d, y0, z0 + r + 2.2 * d, x0 + 0.8 * d, y0, z0 + r + 2.2 * d);
		ev3d_drawLine(x0 - 0.8 * d, y0, z0 + r + 2.2 * d, x0 + 0.8 * d, y0, z0 + r + 0.2 * d);
		ev3d_drawLine(x0 - 0.8 * d, y0, z0 + r + 0.2 * d, x0 + 0.8 * d, y0, z0 + r + 0.2 * d);
	}
}
