/*
ingeometry.c (OpenRCWA, OpenTHFD, OpenSTF)

in geometry ?
*/

#include <math.h>

// 1 : in or border, 0 : out
static int inout3(double x, double y, double tri[3][2], double eps)
{
	int    ret = 0;

	double xmin, xmax, ymin, ymax;
	xmin = xmax = tri[0][0];
	ymin = ymax = tri[0][1];
	for (int i = 1; i < 3; i++) {
		if (tri[i][0] < xmin) xmin = tri[i][0];
		if (tri[i][0] > xmax) xmax = tri[i][0];
		if (tri[i][1] < ymin) ymin = tri[i][1];
		if (tri[i][1] > ymax) ymax = tri[i][1];
	}

	const double zero = eps * (fabs(xmax - xmin) + fabs(ymax - ymin));

	const double x1 = tri[0][0];
	const double x2 = tri[1][0];
	const double x3 = tri[2][0];
	const double y1 = tri[0][1];
	const double y2 = tri[1][1];
	const double y3 = tri[2][1];

	const double det = ((x2 - x1) * (y3 - y1)) - ((x3 - x1) * (y2 - y1));

	if (fabs(det) > (zero * zero)) {
		const double a = + ((x - x1) * (y3 - y1) - (y - y1) * (x3 - x1)) / det;
		const double b = - ((x - x1) * (y2 - y1) - (y - y1) * (x2 - x1)) / det;
		if ((a > -eps) && (b > -eps) && (a + b < 1 + eps)) {
			ret = 1;
		}
	}

	return ret;
}


int ingeometry(double x, double y, double z, int shape, double *g, double eps)
{
	const double zero = 1e-6;
	const double eps2 = eps * eps;

	// rectangle
	if      (shape == 1) {
		if (((x - g[0]) * (x - g[1]) <= eps2) &&
		    ((y - g[2]) * (y - g[3]) <= eps2) &&
		    ((z - g[4]) * (z - g[5]) <= eps2)) {
			return 1;
		}
	}
	// ellipsoid
	else if (shape == 2) {
		const double x0 = (g[0] + g[1]) / 2;
		const double y0 = (g[2] + g[3]) / 2;
		const double z0 = (g[4] + g[5]) / 2;
		const double xr = fabs(g[0] - g[1]) / 2;
		const double yr = fabs(g[2] - g[3]) / 2;
		const double zr = fabs(g[4] - g[5]) / 2;
		if ((x - x0) * (x - x0) / (xr * xr)
		  + (y - y0) * (y - y0) / (yr * yr)
		  + (z - z0) * (z - z0) / (zr * zr) < 1 + zero) {
			return 1;
		}
	}
	// cylinder
	else if (shape == 11) {
		// X cylinder
		const double y0 = (g[2] + g[3]) / 2;
		const double z0 = (g[4] + g[5]) / 2;
		const double yr = fabs(g[2] - g[3]) / 2;
		const double zr = fabs(g[4] - g[5]) / 2;
		if (((x - g[0]) * (x - g[1]) <= eps * eps) &&
		    ((y - y0) * (y - y0) / (yr * yr)
		   + (z - z0) * (z - z0) / (zr * zr) < 1 + zero)) {
			return 1;
		}
	}
	else if (shape == 12) {
		// Y clinder
		const double x0 = (g[0] + g[1]) / 2;
		const double z0 = (g[4] + g[5]) / 2;
		const double xr = fabs(g[0] - g[1]) / 2;
		const double zr = fabs(g[4] - g[5]) / 2;
		if (((y - g[2]) * (y - g[3]) <= eps * eps) &&
		    ((z - z0) * (z - z0) / (zr * zr)
		   + (x - x0) * (x - x0) / (xr * xr) < 1 + zero)) {
			return 1;
		}
	}
	else if (shape == 13) {
		// Z cylinder
		const double x0 = (g[0] + g[1]) / 2;
		const double y0 = (g[2] + g[3]) / 2;
		const double xr = fabs(g[0] - g[1]) / 2;
		const double yr = fabs(g[2] - g[3]) / 2;
		if (((z - g[4]) * (z - g[5]) <= eps * eps) &&
		    ((x - x0) * (x - x0) / (xr * xr)
		   + (y - y0) * (y - y0) / (yr * yr) < 1 + zero)) {
			return 1;
		}
	}
	// pillar
	else if (shape == 31) {
		// X-pillar
		if ((x - g[0]) * (x - g[1]) > eps2) return 0;

		double tri[3][2];
		for (int i = 0; i < 3; i++) {
			tri[i][0] = g[i + 2];
			tri[i][1] = g[i + 5];
		}
		return inout3(y, z, tri, zero);
	}
	else if (shape == 32) {
		// Y-pillar
		if ((y - g[0]) * (y - g[1]) > eps2) return 0;

		double tri[3][2];
		for (int i = 0; i < 3; i++) {
			tri[i][0] = g[i + 2];
			tri[i][1] = g[i + 5];
		}
		return inout3(z, x, tri, zero);
	}
	else if (shape == 33) {
		// Z-pillar
		if ((z - g[0]) * (z - g[1]) > eps2) return 0;

		double tri[3][2];
		for (int i = 0; i < 3; i++) {
			tri[i][0] = g[i + 2];
			tri[i][1] = g[i + 5];
		}
		return inout3(x, y, tri, zero);
	}
	// pyramid
	else if (shape == 41) {
		// X-pyramid
		const double x1  = g[0];
		const double x2  = g[1];
		const double y0  = g[2];
		const double z0  = g[3];
		const double h1y = g[4] / 2;
		const double h1z = g[5] / 2;
		const double h2y = g[6] / 2;
		const double h2z = g[7] / 2;
		const double f = (fabs(x1 - x2) > eps) ? (x - x1) / (x2 - x1) : 0.5;
		const double hy = h1y + f * (h2y - h1y);
		const double hz = h1z + f * (h2z - h1z);
		if (((x - x1) * (x - x2) < eps2) &&
		    (fabs(y - y0) < hy * (1 + zero)) &&
		    (fabs(z - z0) < hz * (1 + zero))) {
			return 1;
		}
	}
	else if (shape == 42) {
		// Y-pyramid
		const double y1  = g[0];
		const double y2  = g[1];
		const double z0  = g[2];
		const double x0  = g[3];
		const double h1z = g[4] / 2;
		const double h1x = g[5] / 2;
		const double h2z = g[6] / 2;
		const double h2x = g[7] / 2;
		const double f = (fabs(y1 - y2) > eps) ? (y - y1) / (y2 - y1) : 0;
		const double hz = h1z + f * (h2z - h1z);
		const double hx = h1x + f * (h2x - h1x);
		if (((y - y1) * (y - y2) < eps2) &&
		    (fabs(z - z0) < hz * (1 + zero)) &&
		    (fabs(x - x0) < hx * (1 + zero))) {
			return 1;
		}
	}
	else if (shape == 43) {
		// Z-pyramid
		const double z1  = g[0];
		const double z2  = g[1];
		const double x0  = g[2];
		const double y0  = g[3];
		const double h1x = g[4] / 2;
		const double h1y = g[5] / 2;
		const double h2x = g[6] / 2;
		const double h2y = g[7] / 2;
		const double f = (fabs(z1 - z2) > eps) ? (z - z1) / (z2 - z1) : 0;
		const double hx = h1x + f * (h2x - h1x);
		const double hy = h1y + f * (h2y - h1y);
		if (((z - z1) * (z - z2) < eps2) &&
		    (fabs(x - x0) < hx * (1 + zero)) &&
		    (fabs(y - y0) < hy * (1 + zero))) {
			return 1;
		}
	}

	// cone
	else if (shape == 51) {
		// X-cone
		const double x1  = g[0];
		const double x2  = g[1];
		const double y0  = g[2];
		const double z0  = g[3];
		const double r1y = g[4] / 2;
		const double r1z = g[5] / 2;
		const double r2y = g[6] / 2;
		const double r2z = g[7] / 2;
		const double f = (fabs(x1 - x2) > eps) ? (x - x1) / (x2 - x1) : 0;
		const double ry = r1y + f * (r2y - r1y);
		const double rz = r1z + f * (r2z - r1z);
		if (((x - x1) * (x - x2) < eps2) &&
		    ((y - y0) * (y - y0) / (ry * ry) + (z - z0) * (z - z0) / (rz * rz) < 1 + zero)) {
			return 1;
		}
	}
	else if (shape == 52) {
		// Y-cone
		const double y1  = g[0];
		const double y2  = g[1];
		const double z0  = g[2];
		const double x0  = g[3];
		const double r1z = g[4] / 2;
		const double r1x = g[5] / 2;
		const double r2z = g[6] / 2;
		const double r2x = g[7] / 2;
		const double f = (fabs(y1 - y2) > eps) ? (y - y1) / (y2 - y1) : 0;
		const double rz = r1z + f * (r2z - r1z);
		const double rx = r1x + f * (r2x - r1x);
		if (((y - y1) * (y - y2) < eps2) &&
		    ((z - z0) * (z - z0) / (rz * rz) + (x - x0) * (x - x0) / (rx * rx) < 1 + zero)) {
			return 1;
		}
	}
	else if (shape == 53) {
		// Z-cone
		const double z1  = g[0];
		const double z2  = g[1];
		const double x0  = g[2];
		const double y0  = g[3];
		const double r1x = g[4] / 2;
		const double r1y = g[5] / 2;
		const double r2x = g[6] / 2;
		const double r2y = g[7] / 2;
		const double f = (fabs(z1 - z2) > eps) ? (z - z1) / (z2 - z1) : 0;
		const double rx = r1x + f * (r2x - r1x);
		const double ry = r1y + f * (r2y - r1y);
		if (((z - z1) * (z - z2) < eps2) &&
		    ((x - x0) * (x - x0) / (rx * rx) + (y - y0) * (y - y0) / (ry * ry) < 1 + zero)) {
			return 1;
		}
	}
/*
	// cone
	else if (shape == 54) {
		double x1 = g[0];
		double x2 = g[1];
		double f = (fabs(x2 - x1) > eps) ? (x2 - x) / (x2 - x1) : 0;
		double y0 = (g[2] + g[3]) / 2;
		double z0 = (g[4] + g[5]) / 2;
		double y1 = y0 + (f * (g[2] - y0));
		double y2 = y0 + (f * (g[3] - y0));
		double z1 = z0 + (f * (g[4] - z0));
		double z2 = z0 + (f * (g[5] - z0));
		double fy = (y - y0) / (fabs(y2 - y1) / 2);
		double fz = (z - z0) / (fabs(z2 - z1) / 2);
		if (((fy * fy) + (fz * fz) < 1 + zero) && ((x - x1) * (x - x2) < eps * eps)) {
			return 1;
		}
	}
	else if (shape == 55) {
		double y1 = g[2];
		double y2 = g[3];
		double f = (fabs(y2 - y1) > eps) ? (y2 - y) / (y2 - y1) : 0;
		double z0 = (g[4] + g[5]) / 2;
		double x0 = (g[0] + g[1]) / 2;
		double z1 = z0 + (f * (g[4] - z0));
		double z2 = z0 + (f * (g[5] - z0));
		double x1 = x0 + (f * (g[0] - x0));
		double x2 = x0 + (f * (g[1] - x0));
		double fz = (z - z0) / (fabs(z2 - z1) / 2);
		double fx = (x - x0) / (fabs(x2 - x1) / 2);
		if (((fz * fz) + (fx * fx) < 1 + zero) && ((y - y1) * (y - y2) < eps * eps)) {
			return 1;
		}
	}
	else if (shape == 56) {
		double z1 = g[4];
		double z2 = g[5];
		double f = (fabs(z2 - z1) > eps) ? (z2 - z) / (z2 - z1) : 0;
		double x0 = (g[0] + g[1]) / 2;
		double y0 = (g[2] + g[3]) / 2;
		double x1 = x0 + (f * (g[0] - x0));
		double x2 = x0 + (f * (g[1] - x0));
		double y1 = y0 + (f * (g[2] - y0));
		double y2 = y0 + (f * (g[3] - y0));
		double fx = (x - x0) / (fabs(x2 - x1) / 2);
		double fy = (y - y0) / (fabs(y2 - y1) / 2);
		if (((fx * fx) + (fy * fy) < 1 + zero) && ((z - z1) * (z - z2) < eps * eps)) {
			return 1;
		}
	}
	// triangle pillar
	else if ((shape == 61) || (shape == 62) || (shape == 63) || (shape == 64) || (shape == 65) || (shape == 66)) {
		double x1 = g[0];
		double x2 = g[1];
		double y1 = g[2];
		double y2 = g[3];
		double z1 = g[4];
		double z2 = g[5];
		double x0 = (x1 + x2) / 2;
		double y0 = (y1 + y2) / 2;
		double z0 = (z1 + z2) / 2;
		if      (shape == 61) {
			double f = (fabs(x2 - x1) > eps) ? (x2 - x) / (x2 - x1) : 0;
			z1 = z0 + (f * (z1 - z0));
			z2 = z0 + (f * (z2 - z0));
		}
		else if (shape == 62) {
			double f = (fabs(x2 - x1) > eps) ? (x2 - x) / (x2 - x1) : 0;
			y1 = y0 + (f * (y1 - y0));
			y2 = y0 + (f * (y2 - y0));
		}
		else if (shape == 63) {
			double f = (fabs(y2 - y1) > eps) ? (y2 - y) / (y2 - y1) : 0;
			x1 = x0 + (f * (x1 - x0));
			x2 = x0 + (f * (x2 - x0));
		}
		else if (shape == 64) {
			double f = (fabs(y2 - y1) > eps) ? (y2 - y) / (y2 - y1) : 0;
			z1 = z0 + (f * (z1 - z0));
			z2 = z0 + (f * (z2 - z0));
		}
		else if (shape == 65) {
			double f = (fabs(z2 - z1) > eps) ? (z2 - z) / (z2 - z1) : 0;
			y1 = y0 + (f * (y1 - y0));
			y2 = y0 + (f * (y2 - y0));
		}
		else if (shape == 66) {
			double f = (fabs(z2 - z1) > eps) ? (z2 - z) / (z2 - z1) : 0;
			x1 = x0 + (f * (x1 - x0));
			x2 = x0 + (f * (x2 - x0));
		}
		if (((x - x1) * (x - x2) < eps * eps) &&
		    ((y - y1) * (y - y2) < eps * eps) &&
		    ((z - z1) * (z - z2) < eps * eps)) {
			return 1;
		}
	}
*/
	return 0;
}
