# -*- coding: utf-8 -*-
"""
comm_Z.py (MPI)
Z境界のHx/Hyを共有する
"""

from mpi4py import MPI
from numba import jit, prange

def share(mode, Hx, Hy,
    SendBuf_Bz_hx, SendBuf_Bz_hy, RecvBuf_Bz_hx, RecvBuf_Bz_hy,
    Bz_ihx, Bz_ihy, Bz_jhx, Bz_jhy,
    Npx, Npy, Npz, Ipx, Ipy, Ipz,
    kMin, kMax, Ni, Nj, Nk, N0):

    assert((mode == 0) or (mode == 1))

    bz = [[Ipz > 0, Ipz < Npz - 1], [Ipz == 0, Ipz == Npz - 1]]
    pz = [[Ipz - 1, Ipz + 1], [Npz - 1, 0]]
    ksend = [kMin + 0, kMax - 1]
    krecv = [kMin - 1, kMax + 0]

    for side in range(2):
        if bz[mode][side]:
            # from H to buffer
            k = ksend[side]
            h_to_buffer(k, Hx, SendBuf_Bz_hx, Bz_ihx, Bz_jhx, Ni, Nj, Nk, N0)
            h_to_buffer(k, Hy, SendBuf_Bz_hy, Bz_ihy, Bz_jhy, Ni, Nj, Nk, N0)

            # MPI
            ipz = pz[mode][side]
            dst = (Ipx * Npy * Npz) + (Ipy * Npz) + ipz
            MPI.COMM_WORLD.Sendrecv(SendBuf_Bz_hx, dst, recvbuf=RecvBuf_Bz_hx, source=dst)
            MPI.COMM_WORLD.Sendrecv(SendBuf_Bz_hy, dst, recvbuf=RecvBuf_Bz_hy, source=dst)

            # from buffer to H
            k = krecv[side]
            buffer_to_h(k, Hx, RecvBuf_Bz_hx, Bz_ihx, Bz_jhx, Ni, Nj, Nk, N0)
            buffer_to_h(k, Hy, RecvBuf_Bz_hy, Bz_ihy, Bz_jhy, Ni, Nj, Nk, N0)

# H to buffer
@jit(cache=True, nopython=True, nogil=True, parallel=True)
def h_to_buffer(k, h, buf, ih, jh, Ni, Nj, Nk, N0):

    for i in prange(ih[0], ih[1] + 1):
        for j in range(jh[0], jh[1] + 1):
            m = (i - ih[0]) * (jh[1] - jh[0] + 1) + (j - jh[0])
            n = (Ni * i) + (Nj * j) + (Nk * k) + N0
            buf[m] = h[n]

# H to buffer
@jit(cache=True, nopython=True, nogil=True, parallel=True)
def buffer_to_h(k, h, buf, ih, jh, Ni, Nj, Nk, N0):

    for i in prange(ih[0], ih[1] + 1):
        for j in range(jh[0], jh[1] + 1):
            m = (i - ih[0]) * (jh[1] - jh[0] + 1) + (j - jh[0])
            n = (Ni * i) + (Nj * j) + (Nk * k) + N0
            h[n] = buf[m]
