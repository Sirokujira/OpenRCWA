/*
setup_mpi.c (MPI)
*/

#include "orcwa.h"

void setup_mpi(void)
{
	const int la = (iABC == 0) ? 1 : cPML.l;

	// X boundary
	if (Npx > 1) {
		Bid.jhy_x[0] = jMin + ((Ipy == 0      ) ? (- la + 1) : 0);
		Bid.jhy_x[1] = jMax + ((Ipy == Npy - 1) ? (+ la - 1) : 0);
		Bid.khy_x[0] = kMin + ((Ipz == 0      ) ? (- la    ) : 0);
		Bid.khy_x[1] = kMax + ((Ipz == Npz - 1) ? (+ la - 1) : -1);
		Bid.jhz_x[0] = jMin + ((Ipy == 0      ) ? (- la    ) : 0);
		Bid.jhz_x[1] = jMax + ((Ipy == Npy - 1) ? (+ la - 1) : -1);
		Bid.khz_x[0] = kMin + ((Ipz == 0      ) ? (- la + 1) : 0);
		Bid.khz_x[1] = kMax + ((Ipz == Npz - 1) ? (+ la - 1) : 0);
		Bid.numhy_x = (Bid.jhy_x[1] - Bid.jhy_x[0] + 1)
		            * (Bid.khy_x[1] - Bid.khy_x[0] + 1);
		Bid.numhz_x = (Bid.jhz_x[1] - Bid.jhz_x[0] + 1)
		            * (Bid.khz_x[1] - Bid.khz_x[0] + 1);
		size_t size_hy_x = Bid.numhy_x * sizeof(real_t);
		size_t size_hz_x = Bid.numhz_x * sizeof(real_t);
		Sendbuf_hy_x = (real_t *)malloc(size_hy_x);
		Recvbuf_hy_x = (real_t *)malloc(size_hy_x);
		Sendbuf_hz_x = (real_t *)malloc(size_hz_x);
		Recvbuf_hz_x = (real_t *)malloc(size_hz_x);
		Bid.ip[0] = iMin - 1;  // - boundary - 1 (recv)
		Bid.ip[1] = iMin;      // - boundary     (send)
		Bid.ip[2] = iMax - 1;  // + boundary - 1 (send)
		Bid.ip[3] = iMax;      // + boundary     (recv)
	}

	// Y boundary
	if (Npy > 1) {
		Bid.khz_y[0] = kMin + ((Ipz == 0      ) ? (- la + 1) : 0);
		Bid.khz_y[1] = kMax + ((Ipz == Npz - 1) ? (+ la - 1) : 0);
		Bid.ihz_y[0] = iMin + ((Ipx == 0      ) ? (- la    ) : 0);
		Bid.ihz_y[1] = iMax + ((Ipx == Npx - 1) ? (+ la - 1) : -1);
		Bid.khx_y[0] = kMin + ((Ipz == 0      ) ? (- la    ) : 0);
		Bid.khx_y[1] = kMax + ((Ipz == Npz - 1) ? (+ la - 1) : -1);
		Bid.ihx_y[0] = iMin + ((Ipx == 0      ) ? (- la + 1) : 0);
		Bid.ihx_y[1] = iMax + ((Ipx == Npx - 1) ? (+ la - 1) : 0);
		Bid.numhz_y = (Bid.khz_y[1] - Bid.khz_y[0] + 1)
		            * (Bid.ihz_y[1] - Bid.ihz_y[0] + 1);
		Bid.numhx_y = (Bid.khx_y[1] - Bid.khx_y[0] + 1)
		            * (Bid.ihx_y[1] - Bid.ihx_y[0] + 1);
		size_t size_hz_y = Bid.numhz_y * sizeof(real_t);
		size_t size_hx_y = Bid.numhx_y * sizeof(real_t);
		Sendbuf_hz_y = (real_t *)malloc(size_hz_y);
		Recvbuf_hz_y = (real_t *)malloc(size_hz_y);
		Sendbuf_hx_y = (real_t *)malloc(size_hx_y);
		Recvbuf_hx_y = (real_t *)malloc(size_hx_y);
		Bid.jp[0] = jMin - 1;  // - boundary - 1 (recv)
		Bid.jp[1] = jMin;      // - boundary     (send)
		Bid.jp[2] = jMax - 1;  // + boundary - 1 (send)
		Bid.jp[3] = jMax;      // + boundary     (recv)
	}

	// Z boundary
	if (Npz > 1) {
		Bid.ihx_z[0] = iMin + ((Ipx == 0      ) ? (- la + 1) : 0);
		Bid.ihx_z[1] = iMax + ((Ipx == Npx - 1) ? (+ la - 1) : 0);
		Bid.jhx_z[0] = jMin + ((Ipy == 0      ) ? (- la    ) : 0);
		Bid.jhx_z[1] = jMax + ((Ipy == Npy - 1) ? (+ la - 1) : -1);
		Bid.ihy_z[0] = iMin + ((Ipx == 0      ) ? (- la    ) : 0);
		Bid.ihy_z[1] = iMax + ((Ipx == Npx - 1) ? (+ la - 1) : -1);
		Bid.jhy_z[0] = jMin + ((Ipy == 0      ) ? (- la + 1) : 0);
		Bid.jhy_z[1] = jMax + ((Ipy == Npy - 1) ? (+ la - 1) : 0);
		Bid.numhx_z = (Bid.ihx_z[1] - Bid.ihx_z[0] + 1)
		            * (Bid.jhx_z[1] - Bid.jhx_z[0] + 1);
		Bid.numhy_z = (Bid.ihy_z[1] - Bid.ihy_z[0] + 1)
		            * (Bid.jhy_z[1] - Bid.jhy_z[0] + 1);
		size_t size_hx_z = Bid.numhx_z * sizeof(real_t);
		size_t size_hy_z = Bid.numhy_z * sizeof(real_t);
		Sendbuf_hx_z = (real_t *)malloc(size_hx_z);
		Recvbuf_hx_z = (real_t *)malloc(size_hx_z);
		Sendbuf_hy_z = (real_t *)malloc(size_hy_z);
		Recvbuf_hy_z = (real_t *)malloc(size_hy_z);
		Bid.kp[0] = kMin - 1;  // - boundary - 1 (recv)
		Bid.kp[1] = kMin;      // - boundary     (send)
		Bid.kp[2] = kMax - 1;  // + boundary - 1 (send)
		Bid.kp[3] = kMax;      // + boundary     (recv)
	}
}
