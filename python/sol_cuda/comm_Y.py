# -*- coding: utf-8 -*-
"""
comm_Y.py (CUDA + MPI)
Y境界のHz/Hxを共有する
"""

import math
from numba import cuda
from mpi4py import MPI

def share(mode, d_Hz, d_Hx,
    SendBuf_By_hz, SendBuf_By_hx, RecvBuf_By_hz, RecvBuf_By_hx,
    d_SendBuf_By_hz, d_SendBuf_By_hx, d_RecvBuf_By_hz, d_RecvBuf_By_hx,
    By_khz, By_khx, By_ihz, By_ihx,
    Npx, Npy, Npz, Ipx, Ipy, Ipz,
    jMin, jMax, Ni, Nj, Nk, N0):

    assert((mode == 0) or (mode == 1))

    khz0, khz1 = By_khz
    ihz0, ihz1 = By_ihz
    khx0, khx1 = By_khx
    ihx0, ihx1 = By_ihx

    # grid, block
    block = (16, 16)
    grid_hz = (math.ceil((khz1 - khz0 + 1) / block[0]),
               math.ceil((ihz1 - ihz0 + 1) / block[1]))
    grid_hx = (math.ceil((khx1 - khx0 + 1) / block[0]),
               math.ceil((ihx1 - ihx0 + 1) / block[1]))

    by = [[Ipy > 0, Ipy < Npy - 1], [Ipy == 0, Ipy == Npy - 1]]
    py = [[Ipy - 1, Ipy + 1], [Npy - 1, 0]]
    jsend = [jMin + 0, jMax - 1]
    jrecv = [jMin - 1, jMax + 0]

    for side in range(2):
        if by[mode][side]:
            # from device memory to host buffer
            j = jsend[side]
            _d2h_gpu[grid_hz, block](j, d_Hz, d_SendBuf_By_hz, khz0, khz1, ihz0, ihz1, Ni, Nj, Nk, N0)
            _d2h_gpu[grid_hx, block](j, d_Hx, d_SendBuf_By_hx, khx0, khx1, ihx0, ihx1, Ni, Nj, Nk, N0)
            d_SendBuf_By_hz.copy_to_host(SendBuf_By_hz)
            d_SendBuf_By_hx.copy_to_host(SendBuf_By_hx)

            # MPI
            ipy = py[mode][side]
            dst = (Ipx * Npy * Npz) + (ipy * Npz) + Ipz
            MPI.COMM_WORLD.Sendrecv(SendBuf_By_hz, dst, recvbuf=RecvBuf_By_hz, source=dst)
            MPI.COMM_WORLD.Sendrecv(SendBuf_By_hx, dst, recvbuf=RecvBuf_By_hx, source=dst)

            # from host buffer to device memory
            j = jrecv[side]
            d_RecvBuf_By_hz = cuda.to_device(RecvBuf_By_hz)
            d_RecvBuf_By_hx = cuda.to_device(RecvBuf_By_hx)
            _h2d_gpu[grid_hz, block](j, d_Hz, d_RecvBuf_By_hz, khz0, khz1, ihz0, ihz1, Ni, Nj, Nk, N0)
            _h2d_gpu[grid_hx, block](j, d_Hx, d_RecvBuf_By_hx, khx0, khx1, ihx0, ihx1, Ni, Nj, Nk, N0)

# (private) (kernel関数)
@cuda.jit(cache=True)
def _d2h_gpu(j, h, buf, k0, k1, i0, i1, Ni, Nj, Nk, N0):

    k, i = cuda.grid(2)
    k += k0
    i += i0
    if (k < k1 + 1) and \
       (i < i1 + 1):
        n = (Ni * i) + (Nj * j) + (Nk * k) + N0
        m = (k - k0) * (i1 - i0 + 1) + (i - i0)
        buf[m] = h[n]

# (private) (kernel関数)
@cuda.jit(cache=True)
def _h2d_gpu(j, h, buf, k0, k1, i0, i1, Ni, Nj, Nk, N0):

    k, i = cuda.grid(2)
    k += k0
    i += i0
    if (k < k1 + 1) and \
       (i < i1 + 1):
        n = (Ni * i) + (Nj * j) + (Nk * k) + N0
        m = (k - k0) * (i1 - i0 + 1) + (i - i0)
        h[n] = buf[m]
