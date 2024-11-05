/*
comm_Z.c (MPI)
share Z plane Hx and Hy
mode = 0/1 : boundary/PBC
*/

#ifdef _MPI
#include <mpi.h>
#endif

#include "orcwa.h"

// H to buffer
static void h_to_buffer(int k, const real_t h[], real_t buf[], const int irange[], const int jrange[])
{
	for (int i = irange[0]; i <= irange[1]; i++) {
	for (int j = jrange[0]; j <= jrange[1]; j++) {
		const int64_t m = (i - irange[0]) * (jrange[1] - jrange[0] + 1) + (j - jrange[0]);
		buf[m] = h[NA(i, j, k)];
	}
	}
}


// buffer to H
static void buffer_to_h(int k, const real_t buf[], real_t h[], const int irange[], const int jrange[])
{
	for (int i = irange[0]; i <= irange[1]; i++) {
	for (int j = jrange[0]; j <= jrange[1]; j++) {
		const int64_t m = (i - irange[0]) * (jrange[1] - jrange[0] + 1) + (j - jrange[0]);
		h[NA(i, j, k)] = buf[m];
	}
	}
}


void comm_Z(int mode)
{
#ifdef _MPI
	MPI_Status status;
	const int tag = 0;
	int bz[][2] = {{Ipz > 0, Ipz < Npz - 1}, {Ipz == 0, Ipz == Npz - 1}};
	int pz[][2] = {{Ipz - 1, Ipz + 1}, {Npz - 1, 0}};
	int krecv[] = {kMin - 1, kMax + 0};
	int ksend[] = {kMin + 0, kMax - 1};
	int k;

	// -Z/+Z boundary
	for (int side = 0; side < 2; side++) {
		if (bz[mode][side]) {
			// H to buffer
			k = ksend[side];
			h_to_buffer(k, Hx, Sendbuf_hx_z, Bid.ihx_z, Bid.jhx_z);
			h_to_buffer(k, Hy, Sendbuf_hy_z, Bid.ihy_z, Bid.jhy_z);

			// MPI
			const int ipz = pz[mode][side];
			const int dst = (Ipx * Npy * Npz) + (Ipy * Npz) + ipz;
			MPI_Sendrecv(Sendbuf_hx_z, Bid.numhx_z, MPI_REAL_T, dst, tag,
			             Recvbuf_hx_z, Bid.numhx_z, MPI_REAL_T, dst, tag, MPI_COMM_WORLD, &status);
			MPI_Sendrecv(Sendbuf_hy_z, Bid.numhy_z, MPI_REAL_T, dst, tag,
			             Recvbuf_hy_z, Bid.numhy_z, MPI_REAL_T, dst, tag, MPI_COMM_WORLD, &status);

			// buffer to H
			k = krecv[side];
			buffer_to_h(k, Recvbuf_hx_z, Hx, Bid.ihx_z, Bid.jhx_z);
			buffer_to_h(k, Recvbuf_hy_z, Hy, Bid.ihy_z, Bid.jhy_z);
		}
	}
#else
	mode = 0;  // dummy
#endif
}
