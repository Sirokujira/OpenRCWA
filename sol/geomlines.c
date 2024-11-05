/*
geomlines.c (OpenRCWA/OpenTHFD/OpenSTF)

make 3d line data from geometries
*/

#include <math.h>

static void line_data(int *nline, double (*gline)[2][3], int *mline, int m,
	double x1, double y1, double z1, double x2, double y2, double z2)
{
	gline[*nline][0][0] = x1;
	gline[*nline][0][1] = y1;
	gline[*nline][0][2] = z1;
	gline[*nline][1][0] = x2;
	gline[*nline][1][1] = y2;
	gline[*nline][1][2] = z2;
	mline[*nline] = m;
	(*nline)++;
}


static void rectangle_data(int *nline, double (*gline)[2][3], int *mline, int m,
	char dir, double c0, double p1, double q1, double p2, double q2)
{
	if      (dir == 'X') {
		line_data(nline, gline, mline, m, c0, p1, q1, c0, p2, q1);
		line_data(nline, gline, mline, m, c0, p2, q1, c0, p2, q2);
		line_data(nline, gline, mline, m, c0, p2, q2, c0, p1, q2);
		line_data(nline, gline, mline, m, c0, p1, q2, c0, p1, q1);
	}
	else if (dir == 'Y') {
		line_data(nline, gline, mline, m, q1, c0, p1, q1, c0, p2);
		line_data(nline, gline, mline, m, q1, c0, p2, q2, c0, p2);
		line_data(nline, gline, mline, m, q2, c0, p2, q2, c0, p1);
		line_data(nline, gline, mline, m, q2, c0, p1, q1, c0, p1);
	}
	else if (dir == 'Z') {
		line_data(nline, gline, mline, m, p1, q1, c0, p2, q1, c0);
		line_data(nline, gline, mline, m, p2, q1, c0, p2, q2, c0);
		line_data(nline, gline, mline, m, p2, q2, c0, p1, q2, c0);
		line_data(nline, gline, mline, m, p1, q2, c0, p1, q1, c0);
	}
}


static void ellipse_data(int *nline, double (*gline)[2][3], int *mline, int m,
	char dir, double c0, double p1, double q1, double p2, double q2, int rdiv)
{
	const double pi = 4 * atan(1);

	for (int i = 0; i < rdiv; i++) {
		double p0 = (p1 + p2) / 2;
		double q0 = (q1 + q2) / 2;
		double pr = fabs(p1 - p2) / 2;
		double qr = fabs(q1 - q2) / 2;
		double a1 = (i + 0) * (2 * pi) / rdiv;
		double a2 = (i + 1) * (2 * pi) / rdiv;
		double x1 = 0, y1 = 0, z1 = 0, x2 = 0, y2 = 0, z2 = 0;
		if      (dir == 'X') {
			x1 = c0;
			y1 = p0 + pr * cos(a1);
			z1 = q0 + qr * sin(a1);
			x2 = c0;
			y2 = p0 + pr * cos(a2);
			z2 = q0 + qr * sin(a2);
		}
		else if (dir == 'Y') {
			y1 = c0;
			z1 = p0 + pr * cos(a1);
			x1 = q0 + qr * sin(a1);
			y2 = c0;
			z2 = p0 + pr * cos(a2);
			x2 = q0 + qr * sin(a2);
		}
		else if (dir == 'Z') {
			z1 = c0;
			x1 = p0 + pr * cos(a1);
			y1 = q0 + qr * sin(a1);
			z2 = c0;
			x2 = p0 + pr * cos(a2);
			y2 = q0 + qr * sin(a2);
		}
		else {
			return;
		}
		line_data(nline, gline, mline, m, x1, y1, z1, x2, y2, z2);
	}
}


/*
mode = 0 : get number of lines
     = 1 : set line data
*/
int geomlines(int mode, int ngeometry, int *shape, int *id, double (*geometry)[8], double (*gline)[2][3], int *mline, double eps)
{
	int nline = 0;
	const int rdiv = 72;

	if (mode == 0) {
		// get number of lines
		for (int n = 0; n < ngeometry; n++) {
			const int s = shape[n];
			nline +=
				(s == 1) ? 12 :
				(s == 2) ? (3 * rdiv) :
				(s == 11) || (s == 12) || (s == 13) ? (2 * rdiv + 4) :
				(s == 31) || (s == 32) || (s == 33) ? 9 :
				(s == 41) || (s == 42) || (s == 43) ? 12 :
				(s == 51) || (s == 52) || (s == 53) ? (2 * rdiv + 4) : 0;
		}
	}
	else if (mode == 1) {
		// set line data
		for (int n = 0; n < ngeometry; n++) {
			const int m = id[n];
			const int s = shape[n];
			double *g = geometry[n];

			double x1 = g[0];
			double x2 = g[1];
			double y1 = g[2];
			double y2 = g[3];
			double z1 = g[4];
			double z2 = g[5];

			double x0 = (x1 + x2) / 2;
			double y0 = (y1 + y2) / 2;
			double z0 = (z1 + z2) / 2;

			if      (s == 1) {
				// rectangle
				const double dx = fabs(x2 - x1);
				const double dy = fabs(y2 - y1);
				const double dz = fabs(z2 - z1);
				// wire
				if      ((dy < eps) && (dz < eps)) {
					line_data(&nline, gline, mline, m, x1, y1, z1, x2, y1, z1);
				}
				else if ((dz < eps) && (dx < eps)) {
					line_data(&nline, gline, mline, m, x1, y1, z1, x1, y2, z1);
				}
				else if ((dx < eps) && (dy < eps)) {
					line_data(&nline, gline, mline, m, x1, y1, z1, x1, y1, z2);
				}
				// plane
				else if (dx < eps) {
					rectangle_data(&nline, gline, mline, m, 'X', x1, y1, z1, y2, z2);
				}
				else if (dy < eps) {
					rectangle_data(&nline, gline, mline, m, 'Y', y1, z1, x1, z2, x2);
				}
				else if (dz < eps) {
					rectangle_data(&nline, gline, mline, m, 'Z', z1, x1, y1, x2, y2);
				}
				// box
				else {
					rectangle_data(&nline, gline, mline, m, 'Z', z1, x1, y1, x2, y2);
					rectangle_data(&nline, gline, mline, m, 'Z', z2, x1, y1, x2, y2);
					line_data(&nline, gline, mline, m, x1, y1, z1, x1, y1, z2);
					line_data(&nline, gline, mline, m, x2, y1, z1, x2, y1, z2);
					line_data(&nline, gline, mline, m, x2, y2, z1, x2, y2, z2);
					line_data(&nline, gline, mline, m, x1, y2, z1, x1, y2, z2);
				}
			}
			else if (s == 2) {
				// sphere
				ellipse_data(&nline, gline, mline, m, 'X', x0, y1, z1, y2, z2, rdiv);
				ellipse_data(&nline, gline, mline, m, 'Y', y0, z1, x1, z2, x2, rdiv);
				ellipse_data(&nline, gline, mline, m, 'Z', z0, x1, y1, x2, y2, rdiv);
			}
			else if (s == 11) {
				// X-cylinder
				ellipse_data(&nline, gline, mline, m, 'X', x1, y1, z1, y2, z2, rdiv);
				ellipse_data(&nline, gline, mline, m, 'X', x2, y1, z1, y2, z2, rdiv);
				line_data(&nline, gline, mline, m, x1, y1, z0, x2, y1, z0);
				line_data(&nline, gline, mline, m, x1, y2, z0, x2, y2, z0);
				line_data(&nline, gline, mline, m, x1, y0, z1, x2, y0, z1);
				line_data(&nline, gline, mline, m, x1, y0, z2, x2, y0, z2);
			}
			else if (s == 12) {
				// Y-cylinder
				ellipse_data(&nline, gline, mline, m, 'Y', y1, z1, x1, z2, x2, rdiv);
				ellipse_data(&nline, gline, mline, m, 'Y', y2, z1, x1, z2, x2, rdiv);
				line_data(&nline, gline, mline, m, x0, y1, z1, x0, y2, z1);
				line_data(&nline, gline, mline, m, x0, y1, z2, x0, y2, z2);
				line_data(&nline, gline, mline, m, x1, y1, z0, x1, y2, z0);
				line_data(&nline, gline, mline, m, x2, y1, z0, x2, y2, z0);
			}
			else if (s == 13) {
				// Z-cylinder
				ellipse_data(&nline, gline, mline, m, 'Z', z1, x1, y1, x2, y2, rdiv);
				ellipse_data(&nline, gline, mline, m, 'Z', z2, x1, y1, x2, y2, rdiv);
				line_data(&nline, gline, mline, m, x1, y0, z1, x1, y0, z2);
				line_data(&nline, gline, mline, m, x2, y0, z1, x2, y0, z2);
				line_data(&nline, gline, mline, m, x0, y1, z1, x0, y1, z2);
				line_data(&nline, gline, mline, m, x0, y2, z1, x0, y2, z2);
			}
			else if (s == 31) {
				// X-pillar
				double *px = g;
				double *py = &g[2];
				double *pz = &g[5];
				for (int i = 0; i < 2; i++) {
					line_data(&nline, gline, mline, m, px[i], py[0], pz[0], px[i], py[1], pz[1]);
					line_data(&nline, gline, mline, m, px[i], py[1], pz[1], px[i], py[2], pz[2]);
					line_data(&nline, gline, mline, m, px[i], py[2], pz[2], px[i], py[0], pz[0]);
				}
				for (int i = 0; i < 3; i++) {
					line_data(&nline, gline, mline, m, px[0], py[i], pz[i], px[1], py[i], pz[i]);
				}
			}
			else if (s == 32) {
				// Y-pillar
				double *py = g;
				double *pz = &g[2];
				double *px = &g[5];
				for (int i = 0; i < 2; i++) {
					line_data(&nline, gline, mline, m, px[0], py[i], pz[0], px[1], py[i], pz[1]);
					line_data(&nline, gline, mline, m, px[1], py[i], pz[1], px[2], py[i], pz[2]);
					line_data(&nline, gline, mline, m, px[2], py[i], pz[2], px[0], py[i], pz[0]);
				}
				for (int i = 0; i < 3; i++) {
					line_data(&nline, gline, mline, m, px[i], py[0], pz[i], px[i], py[1], pz[i]);
				}
			}
			else if (s == 33) {
				// Z-pillar
				double *pz = g;
				double *px = &g[2];
				double *py = &g[5];
				for (int i = 0; i < 2; i++) {
					line_data(&nline, gline, mline, m, px[0], py[0], pz[i], px[1], py[1], pz[i]);
					line_data(&nline, gline, mline, m, px[1], py[1], pz[i], px[2], py[2], pz[i]);
					line_data(&nline, gline, mline, m, px[2], py[2], pz[i], px[0], py[0], pz[i]);
				}
				for (int i = 0; i < 3; i++) {
					line_data(&nline, gline, mline, m, px[i], py[i], pz[0], px[i], py[i], pz[1]);
				}
			}
			else if (s == 41) {
				// X-pyramid
				x1 = g[0];
				x2 = g[1];
				y0 = g[2];
				z0 = g[3];
				const double h1y = g[4] / 2;
				const double h1z = g[5] / 2;
				const double h2y = g[6] / 2;
				const double h2z = g[7] / 2;
				rectangle_data(&nline, gline, mline, m, 'X', x1, y0 - h1y, z0 - h1z, y0 + h1y, z0 + h1z);
				rectangle_data(&nline, gline, mline, m, 'X', x2, y0 - h2y, z0 - h2z, y0 + h2y, z0 + h2z);
				line_data(&nline, gline, mline, m, x1, y0 - h1y, z0 - h1z, x2, y0 - h2y, z0 - h2z);
				line_data(&nline, gline, mline, m, x1, y0 - h1y, z0 + h1z, x2, y0 - h2y, z0 + h2z);
				line_data(&nline, gline, mline, m, x1, y0 + h1y, z0 - h1z, x2, y0 + h2y, z0 - h2z);
				line_data(&nline, gline, mline, m, x1, y0 + h1y, z0 + h1z, x2, y0 + h2y, z0 + h2z);
			}
			else if (s == 42) {
				// Y-pyramid
				y1 = g[0];
				y2 = g[1];
				z0 = g[2];
				x0 = g[3];
				const double h1z = g[4] / 2;
				const double h1x = g[5] / 2;
				const double h2z = g[6] / 2;
				const double h2x = g[7] / 2;
				rectangle_data(&nline, gline, mline, m, 'Y', y1, z0 - h1z, x0 - h1x, z0 + h1z, x0 + h1x);
				rectangle_data(&nline, gline, mline, m, 'Y', y2, z0 - h2z, x0 - h2x, z0 + h2z, x0 + h2x);
				line_data(&nline, gline, mline, m, x0 - h1x, y1, z0 - h1z, x0 - h2x, y2, z0 - h2z);
				line_data(&nline, gline, mline, m, x0 - h1x, y1, z0 + h1z, x0 - h2x, y2, z0 + h2z);
				line_data(&nline, gline, mline, m, x0 + h1x, y1, z0 - h1z, x0 + h2x, y2, z0 - h2z);
				line_data(&nline, gline, mline, m, x0 + h1x, y1, z0 + h1z, x0 + h2x, y2, z0 + h2z);
			}
			else if (s == 43) {
				// Z-pyramid
				z1 = g[0];
				z2 = g[1];
				x0 = g[2];
				y0 = g[3];
				const double h1x = g[4] / 2;
				const double h1y = g[5] / 2;
				const double h2x = g[6] / 2;
				const double h2y = g[7] / 2;
				rectangle_data(&nline, gline, mline, m, 'Z', z1, x0 - h1x, y0 - h1y, x0 + h1x, y0 + h1y);
				rectangle_data(&nline, gline, mline, m, 'Z', z2, x0 - h2x, y0 - h2y, x0 + h2x, y0 + h2y);
				line_data(&nline, gline, mline, m, x0 - h1x, y0 - h1y, z1, x0 - h2x, y0 - h2y, z2);
				line_data(&nline, gline, mline, m, x0 - h1x, y0 + h1y, z1, x0 - h2x, y0 + h2y, z2);
				line_data(&nline, gline, mline, m, x0 + h1x, y0 - h1y, z1, x0 + h2x, y0 - h2y, z2);
				line_data(&nline, gline, mline, m, x0 + h1x, y0 + h1y, z1, x0 + h2x, y0 + h2y, z2);
			}
			else if (s == 51) {
				// X-cone
				x1 = g[0];
				x2 = g[1];
				y0 = g[2];
				z0 = g[3];
				const double r1y = g[4] / 2;
				const double r1z = g[5] / 2;
				const double r2y = g[6] / 2;
				const double r2z = g[7] / 2;
				ellipse_data(&nline, gline, mline, m, 'X', x1, y0 - r1y, z0 - r1z, y0 + r1y, z0 + r1z, rdiv);
				ellipse_data(&nline, gline, mline, m, 'X', x2, y0 - r2y, z0 - r2z, y0 + r2y, z0 + r2z, rdiv);
				line_data(&nline, gline, mline, m, x1, y0,       z0 + r1z, x2, y0,       z0 + r2z);
				line_data(&nline, gline, mline, m, x1, y0,       z0 - r1z, x2, y0,       z0 - r2z);
				line_data(&nline, gline, mline, m, x1, y0 + r1y, z0,       x2, y0 + r2y, z0      );
				line_data(&nline, gline, mline, m, x1, y0 - r1y, z0,       x2, y0 - r2y, z0      );
			}
			else if (s == 52) {
				// Y-cone
				y1 = g[0];
				y2 = g[1];
				z0 = g[2];
				x0 = g[3];
				const double r1z = g[4] / 2;
				const double r1x = g[5] / 2;
				const double r2z = g[6] / 2;
				const double r2x = g[7] / 2;
				ellipse_data(&nline, gline, mline, m, 'Y', y1, z0 - r1z, x0 - r1x, z0 + r1z, x0 + r1x, rdiv);
				ellipse_data(&nline, gline, mline, m, 'Y', y2, z0 - r2z, x0 - r2x, z0 + r2z, x0 + r2x, rdiv);
				line_data(&nline, gline, mline, m, x0,       y1, z0 + r1z, x0,       y2, z0 + r2z);
				line_data(&nline, gline, mline, m, x0,       y1, z0 - r1z, x0,       y2, z0 - r2z);
				line_data(&nline, gline, mline, m, x0 + r1x, y1, z0,       x0 + r2x, y2, z0      );
				line_data(&nline, gline, mline, m, x0 - r1x, y1, z0,       x0 - r2x, y2, z0      );
			}
			else if (s == 53) {
				// Z-cone
				z1 = g[0];
				z2 = g[1];
				x0 = g[2];
				y0 = g[3];
				const double r1x = g[4] / 2;
				const double r1y = g[5] / 2;
				const double r2x = g[6] / 2;
				const double r2y = g[7] / 2;
				ellipse_data(&nline, gline, mline, m, 'Z', z1, x0 - r1x, y0 - r1y, x0 + r1x, y0 + r1y, rdiv);
				ellipse_data(&nline, gline, mline, m, 'Z', z2, x0 - r2x, y0 - r2y, x0 + r2x, y0 + r2y, rdiv);
				line_data(&nline, gline, mline, m, x0,       y0 + r1y, z1, x0,       y0 + r2y, z2);
				line_data(&nline, gline, mline, m, x0,       y0 - r1y, z1, x0,       y0 - r2y, z2);
				line_data(&nline, gline, mline, m, x0 + r1x, y0,       z1, x0 + r2x, y0,       z2);
				line_data(&nline, gline, mline, m, x0 - r1x, y0,       z1, x0 - r2x, y0,       z2);
			}
		}
	}

	return nline;
}
