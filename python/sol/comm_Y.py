# -*- coding: utf-8 -*-
"""
comm_Y.py (MPI)
Y境界のHz/Hxを共有する
"""

from mpi4py import MPI
from numba import jit, prange

def share(mode, Hz, Hx,
    SendBuf_By_hz, SendBuf_By_hx, RecvBuf_By_hz, RecvBuf_By_hx,
    By_khz, By_khx, By_ihz, By_ihx,
    Npx, Npy, Npz, Ipx, Ipy, Ipz,
    jMin, jMax, Ni, Nj, Nk, N0):

    assert((mode == 0) or (mode == 1))

    by = [[Ipy > 0, Ipy < Npy - 1], [Ipy == 0, Ipy == Npy - 1]]
    py = [[Ipy - 1, Ipy + 1], [Npy - 1, 0]]
    jsend = [jMin + 0, jMax - 1]
    jrecv = [jMin - 1, jMax + 0]

    for side in range(2):
        if by[mode][side]:
            # from H to buffer
            j = jsend[side]
            h_to_buffer(j, Hz, SendBuf_By_hz, By_khz, By_ihz, Ni, Nj, Nk, N0)
            h_to_buffer(j, Hx, SendBuf_By_hx, By_khx, By_ihx, Ni, Nj, Nk, N0)

            # MPI
            ipy = py[mode][side]
            dst = (Ipx * Npy * Npz) + (ipy * Npz) + Ipz
            MPI.COMM_WORLD.Sendrecv(SendBuf_By_hz, dst, recvbuf=RecvBuf_By_hz, source=dst)
            MPI.COMM_WORLD.Sendrecv(SendBuf_By_hx, dst, recvbuf=RecvBuf_By_hx, source=dst)

            # from buffer to H
            j = jrecv[side]
            buffer_to_h(j, Hz, RecvBuf_By_hz, By_khz, By_ihz, Ni, Nj, Nk, N0)
            buffer_to_h(j, Hx, RecvBuf_By_hx, By_khx, By_ihx, Ni, Nj, Nk, N0)

# H to buffer
@jit(cache=True, nopython=True, nogil=True, parallel=True)
def h_to_buffer(j, h, buf, kh, ih, Ni, Nj, Nk, N0):

    for k in prange(kh[0], kh[1] + 1):
        for i in range(ih[0], ih[1] + 1):
            m = (k - kh[0]) * (ih[1] - ih[0] + 1) + (i - ih[0])
            n = (Ni * i) + (Nj * j) + (Nk * k) + N0
            buf[m] = h[n]

# H to buffer
@jit(cache=True, nopython=True, nogil=True, parallel=True)
def buffer_to_h(j, h, buf, kh, ih, Ni, Nj, Nk, N0):

    for k in prange(kh[0], kh[1] + 1):
        for i in range(ih[0], ih[1] + 1):
            m = (k - kh[0]) * (ih[1] - ih[0] + 1) + (i - ih[0])
            n = (Ni * i) + (Nj * j) + (Nk * k) + N0
            h[n] = buf[m]
