/*
fieldnode.cu

calculate field and DFT at a node
*/

__host__ __device__
static void fieldnode(
	int nx, int ny, int nz, int imin, int imax,
	int i, int j, int k,
	real_t *ex, real_t *ey, real_t *ez,
	real_t *hx, real_t *hy, real_t *hz,
	d_complex_t *cex, d_complex_t *cey, d_complex_t *cez,
	d_complex_t *chx, d_complex_t *chy, d_complex_t *chz,
	d_complex_t fe, d_complex_t fh, param_t *p)
{
	real_t e_x, e_y, e_z;
	real_t h_x, h_y, h_z;

	nx = nx;
	ny = ny;
	nz = nz;
	if      (i <= imin) {
		i = imin;
		e_x = (ex[LA(p, i + 0, j, k)] * 3 - ex[LA(p, i + 1, j, k)] * 1) / 2;
	}
	else if (i >= imax) {
		i = imax;
		e_x = (ex[LA(p, i - 1, j, k)] * 3 - ex[LA(p, i - 2, j, k)] * 1) / 2;
	}
	else {
		e_x = (ex[LA(p, i, j, k)] + ex[LA(p, i - 1, j, k)]) / 2;
	}

	if      (j <=    0) {
		j =    0;
		e_y = (ey[LA(p, i, j + 0, k)] * 3 - ey[LA(p, i, j + 1, k)] * 1) / 2;
	}
	else if (j >=   ny) {
		j =   ny;
		e_y = (ey[LA(p, i, j - 1, k)] * 3 - ey[LA(p, i, j - 2, k)] * 1) / 2;
	}
	else {
		e_y = (ey[LA(p, i, j, k)] + ey[LA(p, i, j - 1, k)]) / 2;
	}

	if      (k <=    0) {
		k =    0;
		e_z = (ez[LA(p, i, j, k + 0)] * 3 - ez[LA(p, i, j, k + 1)] * 1) / 2;
	}
	else if (k >=   nz) {
		k =   nz;
		e_z = (ez[LA(p, i, j, k - 1)] * 3 - ez[LA(p, i, j, k - 2)] * 1) / 2;
	}
	else {
		e_z = (ez[LA(p, i, j, k)] + ez[LA(p, i, j, k - 1)]) / 2;
	}

	h_x = (hx[LA(p, i, j, k)] + hx[LA(p, i, j - 1, k    )] + hx[LA(p, i,     j,     k - 1)] + hx[LA(p, i,     j - 1, k - 1)]) / 4;
	h_y = (hy[LA(p, i, j, k)] + hy[LA(p, i, j,     k - 1)] + hy[LA(p, i - 1, j,     k    )] + hy[LA(p, i - 1, j,     k - 1)]) / 4;
	h_z = (hz[LA(p, i, j, k)] + hz[LA(p, i - 1, j, k    )] + hz[LA(p, i,     j - 1, k    )] + hz[LA(p, i - 1, j - 1, k    )]) / 4;

	// E
	cex->r += e_x * fe.r;
	cex->i += e_x * fe.i;
	cey->r += e_y * fe.r;
	cey->i += e_y * fe.i;
	cez->r += e_z * fe.r;
	cez->i += e_z * fe.i;

	// H
	chx->r += h_x * fh.r;
	chx->i += h_x * fh.i;
	chy->r += h_y * fh.r;
	chy->i += h_y * fh.i;
	chz->r += h_z * fh.r;
	chz->i += h_z * fh.i;
}
