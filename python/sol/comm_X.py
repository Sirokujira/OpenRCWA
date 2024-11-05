# -*- coding: utf-8 -*-
"""
comm_X.py (MPI)
X境界のHy/Hzを共有する
"""

from mpi4py import MPI
from numba import jit, prange

def share(mode, Hy, Hz,
    SendBuf_Bx_hy, SendBuf_Bx_hz, RecvBuf_Bx_hy, RecvBuf_Bx_hz,
    Bx_jhy, Bx_jhz, Bx_khy, Bx_khz,
    Npx, Npy, Npz, Ipx, Ipy, Ipz,
    iMin, iMax, Ni, Nj, Nk, N0):

    assert((mode == 0) or (mode == 1))

    jhy0, jhy1 = Bx_jhy
    khy0, khy1 = Bx_khy
    jhz0, jhz1 = Bx_jhz
    khz0, khz1 = Bx_khz

    bx = [[Ipx > 0, Ipx < Npx - 1], [Ipx == 0, Ipx == Npx - 1]]
    px = [[Ipx - 1, Ipx + 1], [Npx - 1, 0]]
    isend = [iMin + 0, iMax - 1]
    irecv = [iMin - 1, iMax + 0]

    for side in range(2):
        if bx[mode][side]:
            # from H to buffer
            i = isend[side]
            h_to_buffer(i, Hy, SendBuf_Bx_hy, jhy0, jhy1, khy0, khy1, Ni, Nj, Nk, N0)
            h_to_buffer(i, Hz, SendBuf_Bx_hz, jhz0, jhz1, khz0, khz1, Ni, Nj, Nk, N0)

            # MPI
            ipx = px[mode][side]
            dst = (ipx * Npy * Npz) + (Ipy * Npz) + Ipz
            MPI.COMM_WORLD.Sendrecv(SendBuf_Bx_hy, dst, recvbuf=RecvBuf_Bx_hy, source=dst)
            MPI.COMM_WORLD.Sendrecv(SendBuf_Bx_hz, dst, recvbuf=RecvBuf_Bx_hz, source=dst)

            # from buffer to H
            i = irecv[side]
            buffer_to_h(i, Hy, RecvBuf_Bx_hy, jhy0, jhy1, khy0, khy1, Ni, Nj, Nk, N0)
            buffer_to_h(i, Hz, RecvBuf_Bx_hz, jhz0, jhz1, khz0, khz1, Ni, Nj, Nk, N0)

# H to buffer
@jit(cache=True, nopython=True, nogil=True, parallel=True)
def h_to_buffer(i, h, buf, j0, j1, k0, k1, Ni, Nj, Nk, N0):

    for j in prange(j0, j1 + 1):
        for k in range(k0, k1 + 1):
            m = (j - j0) * (k1 - k0 + 1) + (k - k0)
            n = (Ni * i) + (Nj * j) + (Nk * k) + N0
            buf[m] = h[n]

# H to buffer
@jit(cache=True, nopython=True, nogil=True, parallel=True)
def buffer_to_h(i, h, buf, j0, j1, k0, k1, Ni, Nj, Nk, N0):

    for j in prange(j0, j1 + 1):
        for k in range(k0, k1 + 1):
            m = (j - j0) * (k1 - k0 + 1) + (k - k0)
            n = (Ni * i) + (Nj * j) + (Nk * k) + N0
            h[n] = buf[m]
