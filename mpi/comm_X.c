/*
comm_X.c (MPI)
share X plane Hy and Hz
mode = 0/1 : boundary/PBC
*/

#ifdef _MPI
#include <mpi.h>
#endif

#include "orcwa.h"

// H to buffer
static void h_to_buffer(int i, const real_t h[], real_t buf[], const int jrange[], const int krange[])
{
	for (int j = jrange[0]; j <= jrange[1]; j++) {
	for (int k = krange[0]; k <= krange[1]; k++) {
		const int64_t m = (j - jrange[0]) * (krange[1] - krange[0] + 1) + (k - krange[0]);
		buf[m] = h[NA(i, j, k)];
	}
	}
}


// buffer to H
static void buffer_to_h(int i, const real_t buf[], real_t h[], const int jrange[], const int krange[])
{
	for (int j = jrange[0]; j <= jrange[1]; j++) {
	for (int k = krange[0]; k <= krange[1]; k++) {
		const int64_t m = (j - jrange[0]) * (krange[1] - krange[0] + 1) + (k - krange[0]);
		h[NA(i, j, k)] = buf[m];
	}
	}
}


void comm_X(int mode)
{
#ifdef _MPI
	MPI_Status status;
	const int tag = 0;
	int bx[][2] = {{Ipx > 0, Ipx < Npx - 1}, {Ipx == 0, Ipx == Npx - 1}};
	int px[][2] = {{Ipx - 1, Ipx + 1}, {Npx - 1, 0}};
	int irecv[] = {iMin - 1, iMax + 0};
	int isend[] = {iMin + 0, iMax - 1};
	int i;

	// -X/+X boundary
	for (int side = 0; side < 2; side++) {
		if (bx[mode][side]) {
			// H to buffer
			i = isend[side];
			h_to_buffer(i, Hy, Sendbuf_hy_x, Bid.jhy_x, Bid.khy_x);
			h_to_buffer(i, Hz, Sendbuf_hz_x, Bid.jhz_x, Bid.khz_x);

			// MPI
			const int ipx = px[mode][side];
			const int dst = (ipx * Npy * Npz) + (Ipy * Npz) + Ipz;
			MPI_Sendrecv(Sendbuf_hy_x, Bid.numhy_x, MPI_REAL_T, dst, tag,
			             Recvbuf_hy_x, Bid.numhy_x, MPI_REAL_T, dst, tag, MPI_COMM_WORLD, &status);
			MPI_Sendrecv(Sendbuf_hz_x, Bid.numhz_x, MPI_REAL_T, dst, tag,
			             Recvbuf_hz_x, Bid.numhz_x, MPI_REAL_T, dst, tag, MPI_COMM_WORLD, &status);

			// buffer to H
			i = irecv[side];
			buffer_to_h(i, Recvbuf_hy_x, Hy, Bid.jhy_x, Bid.khy_x);
			buffer_to_h(i, Recvbuf_hz_x, Hz, Bid.jhz_x, Bid.khz_x);
		}
	}
#else
	mode = 0;  // dummy
#endif
}
