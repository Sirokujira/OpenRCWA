/*
comm_Y.c (MPI)
share Y plane Hz and Hx
mode = 0/1 : boundary/PBC
*/

#ifdef _MPI
#include <mpi.h>
#endif

#include "orcwa.h"

// H to buffer
static void h_to_buffer(int j, const real_t h[], real_t buf[], const int krange[], const int irange[])
{
	for (int k = krange[0]; k <= krange[1]; k++) {
	for (int i = irange[0]; i <= irange[1]; i++) {
		const int64_t m = (k - krange[0]) * (irange[1] - irange[0] + 1) + (i - irange[0]);
		buf[m] = h[NA(i, j, k)];
	}
	}
}


// buffer to H
static void buffer_to_h(int j, const real_t buf[], real_t h[], const int krange[], const int irange[])
{
	for (int k = krange[0]; k <= krange[1]; k++) {
	for (int i = irange[0]; i <= irange[1]; i++) {
		const int64_t m = (k - krange[0]) * (irange[1] - irange[0] + 1) + (i - irange[0]);
		h[NA(i, j, k)] = buf[m];
	}
	}
}


void comm_Y(int mode)
{
#ifdef _MPI
	MPI_Status status;
	const int tag = 0;
	int by[][2] = {{Ipy > 0, Ipy < Npy - 1}, {Ipy == 0, Ipy == Npy - 1}};
	int py[][2] = {{Ipy - 1, Ipy + 1}, {Npy - 1, 0}};
	int jrecv[] = {jMin - 1, jMax + 0};
	int jsend[] = {jMin + 0, jMax - 1};
	int j;

	// -Y/+Y boundary
	for (int side = 0; side < 2; side++) {
		if (by[mode][side]) {
			// H to buffer
			j = jsend[side];
			h_to_buffer(j, Hz, Sendbuf_hz_y, Bid.khz_y, Bid.ihz_y);
			h_to_buffer(j, Hx, Sendbuf_hx_y, Bid.khx_y, Bid.ihx_y);

			// MPI
			const int ipy = py[mode][side];
			const int dst = (Ipx * Npy * Npz) + (ipy * Npz) + Ipz;
			MPI_Sendrecv(Sendbuf_hz_y, Bid.numhz_y, MPI_REAL_T, dst, tag,
			             Recvbuf_hz_y, Bid.numhz_y, MPI_REAL_T, dst, tag, MPI_COMM_WORLD, &status);
			MPI_Sendrecv(Sendbuf_hx_y, Bid.numhx_y, MPI_REAL_T, dst, tag,
			             Recvbuf_hx_y, Bid.numhx_y, MPI_REAL_T, dst, tag, MPI_COMM_WORLD, &status);

			// buffer to H
			j = jrecv[side];
			buffer_to_h(j, Recvbuf_hz_y, Hz, Bid.khz_y, Bid.ihz_y);
			buffer_to_h(j, Recvbuf_hx_y, Hx, Bid.khx_y, Bid.ihx_y);
		}
	}
#else
	mode = 0;  // dummy
#endif
}
