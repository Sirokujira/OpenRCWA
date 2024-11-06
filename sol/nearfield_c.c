/*
nearfield_c.c

near field at node (i, j, k) (complex)
*/

#include "orcwa.h"
#include "complex.h"

// E at node
void NodeE_c(int ifreq, int i, int j, int k, d_complex_t *cex, d_complex_t *cey, d_complex_t *cez)
{
	d_complex_t c1, c2;

	if      (i <= 0) {
		i = 0;
		c1 = d_complex(cEx_r[(ifreq * NN) + NA(i + 0, j, k)],
		               cEx_i[(ifreq * NN) + NA(i + 0, j, k)]);
		c2 = d_complex(cEx_r[(ifreq * NN) + NA(i + 1, j, k)],
		               cEx_i[(ifreq * NN) + NA(i + 1, j, k)]);
		*cex = d_sub(d_rmul(1.5, c1), d_rmul(0.5, c2));
	}
	else if (i >= Nx) {
		i = Nx;
		c1 = d_complex(cEx_r[(ifreq * NN) + NA(i - 1, j, k)],
		               cEx_i[(ifreq * NN) + NA(i - 1, j, k)]);
		c2 = d_complex(cEx_r[(ifreq * NN) + NA(i - 2, j, k)],
		               cEx_i[(ifreq * NN) + NA(i - 2, j, k)]);
		*cex = d_sub(d_rmul(1.5, c1), d_rmul(0.5, c2));
	}
	else {
		c1 = d_complex(cEx_r[(ifreq * NN) + NA(i + 0, j, k)],
		               cEx_i[(ifreq * NN) + NA(i + 0, j, k)]);
		c2 = d_complex(cEx_r[(ifreq * NN) + NA(i - 1, j, k)],
		               cEx_i[(ifreq * NN) + NA(i - 1, j, k)]);
		*cex = d_rmul(0.5, d_add(c1, c2));
	}

	if      (j <= 0) {
		j = 0;
		c1 = d_complex(cEy_r[(ifreq * NN) + NA(i, j + 0, k)],
		               cEy_i[(ifreq * NN) + NA(i, j + 0, k)]);
		c2 = d_complex(cEy_r[(ifreq * NN) + NA(i, j + 1, k)],
		               cEy_i[(ifreq * NN) + NA(i, j + 1, k)]);
		*cey = d_sub(d_rmul(1.5, c1), d_rmul(0.5, c2));
	}
	else if (j >= Ny) {
		j = Ny;
		c1 = d_complex(cEy_r[(ifreq * NN) + NA(i, j - 1, k)],
		               cEy_i[(ifreq * NN) + NA(i, j - 1, k)]);
		c2 = d_complex(cEy_r[(ifreq * NN) + NA(i, j - 2, k)],
		               cEy_i[(ifreq * NN) + NA(i, j - 2, k)]);
		*cey = d_sub(d_rmul(1.5, c1), d_rmul(0.5, c2));
	}
	else {
		c1 = d_complex(cEy_r[(ifreq * NN) + NA(i, j + 0, k)],
		               cEy_i[(ifreq * NN) + NA(i, j + 0, k)]);
		c2 = d_complex(cEy_r[(ifreq * NN) + NA(i, j - 1, k)],
		               cEy_i[(ifreq * NN) + NA(i, j - 1, k)]);
		*cey = d_rmul(0.5, d_add(c1, c2));
	}

	if      (k <= 0) {
		k = 0;
		c1 = d_complex(cEz_r[(ifreq * NN) + NA(i, j, k + 0)],
		               cEz_i[(ifreq * NN) + NA(i, j, k + 0)]);
		c2 = d_complex(cEz_r[(ifreq * NN) + NA(i, j, k + 1)],
		               cEz_i[(ifreq * NN) + NA(i, j, k + 1)]);
		*cez = d_sub(d_rmul(1.5, c1), d_rmul(0.5, c2));
	}
	else if (k >= Nz) {
		k = Nz;
		c1 = d_complex(cEz_r[(ifreq * NN) + NA(i, j, k - 1)],
		               cEz_i[(ifreq * NN) + NA(i, j, k - 1)]);
		c2 = d_complex(cEz_r[(ifreq * NN) + NA(i, j, k - 2)],
		               cEz_i[(ifreq * NN) + NA(i, j, k - 2)]);
		*cez = d_sub(d_rmul(1.5, c1), d_rmul(0.5, c2));
	}
	else {
		c1 = d_complex(cEz_r[(ifreq * NN) + NA(i, j, k + 0)],
		               cEz_i[(ifreq * NN) + NA(i, j, k + 0)]);
		c2 = d_complex(cEz_r[(ifreq * NN) + NA(i, j, k - 1)],
		               cEz_i[(ifreq * NN) + NA(i, j, k - 1)]);
		*cez = d_rmul(0.5, d_add(c1, c2));
	}
}


// H at node
void NodeH_c(int ifreq, int i, int j, int k, d_complex_t *chx, d_complex_t *chy, d_complex_t *chz)
{
	d_complex_t c1, c2, c3, c4;

	c1 = d_complex(cHx_r[(ifreq * NN) + NA(i,     j,     k    )],
	               cHx_i[(ifreq * NN) + NA(i,     j,     k    )]);
	c2 = d_complex(cHx_r[(ifreq * NN) + NA(i,     j - 1, k    )],
	               cHx_i[(ifreq * NN) + NA(i,     j - 1, k    )]);
	c3 = d_complex(cHx_r[(ifreq * NN) + NA(i,     j,     k - 1)],
	               cHx_i[(ifreq * NN) + NA(i,     j,     k - 1)]);
	c4 = d_complex(cHx_r[(ifreq * NN) + NA(i,     j - 1, k - 1)],
	               cHx_i[(ifreq * NN) + NA(i,     j - 1, k - 1)]);
	*chx = d_rmul(0.25, d_add4(c1, c2, c3, c4));

	c1 = d_complex(cHy_r[(ifreq * NN) + NA(i,     j,     k    )],
	               cHy_i[(ifreq * NN) + NA(i,     j,     k    )]);
	c2 = d_complex(cHy_r[(ifreq * NN) + NA(i,     j,     k - 1)],
	               cHy_i[(ifreq * NN) + NA(i,     j,     k - 1)]);
	c3 = d_complex(cHy_r[(ifreq * NN) + NA(i - 1, j,     k    )],
	               cHy_i[(ifreq * NN) + NA(i - 1, j,     k    )]);
	c4 = d_complex(cHy_r[(ifreq * NN) + NA(i - 1, j,     k - 1)],
	               cHy_i[(ifreq * NN) + NA(i - 1, j,     k - 1)]);
	*chy = d_rmul(0.25, d_add4(c1, c2, c3, c4));

	c1 = d_complex(cHz_r[(ifreq * NN) + NA(i,     j,     k    )],
	               cHz_i[(ifreq * NN) + NA(i,     j,     k    )]);
	c2 = d_complex(cHz_r[(ifreq * NN) + NA(i - 1, j,     k    )],
	               cHz_i[(ifreq * NN) + NA(i - 1, j,     k    )]);
	c3 = d_complex(cHz_r[(ifreq * NN) + NA(i,     j - 1, k    )],
	               cHz_i[(ifreq * NN) + NA(i,     j - 1, k    )]);
	c4 = d_complex(cHz_r[(ifreq * NN) + NA(i - 1, j - 1, k    )],
	               cHz_i[(ifreq * NN) + NA(i - 1, j - 1, k    )]);
	*chz = d_rmul(0.25, d_add4(c1, c2, c3, c4));
}
