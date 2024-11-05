# -*- coding: utf-8 -*-
"""
comm_X.py (CUDA + MPI)
X境界のHy/Hzを共有する
"""

import math
from numba import cuda
from mpi4py import MPI

def share(mode, d_Hy, d_Hz,
    SendBuf_Bx_hy, SendBuf_Bx_hz, RecvBuf_Bx_hy, RecvBuf_Bx_hz,
    d_SendBuf_Bx_hy, d_SendBuf_Bx_hz, d_RecvBuf_Bx_hy, d_RecvBuf_Bx_hz,
    Bx_jhy, Bx_jhz, Bx_khy, Bx_khz,
    Npx, Npy, Npz, Ipx, Ipy, Ipz,
    iMin, iMax, Ni, Nj, Nk, N0):

    assert((mode == 0) or (mode == 1))

    jhy0, jhy1 = Bx_jhy
    khy0, khy1 = Bx_khy
    jhz0, jhz1 = Bx_jhz
    khz0, khz1 = Bx_khz

    # grid, block
    block = (16, 16)
    grid_hy = (math.ceil((jhy1 - jhy0 + 1) / block[0]),
               math.ceil((khy1 - khy0 + 1) / block[1]))
    grid_hz = (math.ceil((jhz1 - jhz0 + 1) / block[0]),
               math.ceil((khz1 - khz0 + 1) / block[1]))

    bx = [[Ipx > 0, Ipx < Npx - 1], [Ipx == 0, Ipx == Npx - 1]]
    px = [[Ipx - 1, Ipx + 1], [Npx - 1, 0]]
    isend = [iMin + 0, iMax - 1]
    irecv = [iMin - 1, iMax + 0]

    for side in range(2):
        if bx[mode][side]:
            # from device memory to host buffer
            i = isend[side]
            _d2h_gpu[grid_hy, block](i, d_Hy, d_SendBuf_Bx_hy, jhy0, jhy1, khy0, khy1, Ni, Nj, Nk, N0)
            _d2h_gpu[grid_hz, block](i, d_Hz, d_SendBuf_Bx_hz, jhz0, jhz1, khz0, khz1, Ni, Nj, Nk, N0)
            d_SendBuf_Bx_hy.copy_to_host(SendBuf_Bx_hy)
            d_SendBuf_Bx_hz.copy_to_host(SendBuf_Bx_hz)

            # MPI
            ipx = px[mode][side]
            dst = (ipx * Npy * Npz) + (Ipy * Npz) + Ipz
            MPI.COMM_WORLD.Sendrecv(SendBuf_Bx_hy, dst, recvbuf=RecvBuf_Bx_hy, source=dst)
            MPI.COMM_WORLD.Sendrecv(SendBuf_Bx_hz, dst, recvbuf=RecvBuf_Bx_hz, source=dst)

            # from host buffer to device memory
            i = irecv[side]
            d_RecvBuf_Bx_hy = cuda.to_device(RecvBuf_Bx_hy)
            d_RecvBuf_Bx_hz = cuda.to_device(RecvBuf_Bx_hz)
            _h2d_gpu[grid_hy, block](i, d_Hy, d_RecvBuf_Bx_hy, jhy0, jhy1, khy0, khy1, Ni, Nj, Nk, N0)
            _h2d_gpu[grid_hz, block](i, d_Hz, d_RecvBuf_Bx_hz, jhz0, jhz1, khz0, khz1, Ni, Nj, Nk, N0)

# (private) (kernel関数)
@cuda.jit(cache=True)
def _d2h_gpu(i, h, buf, j0, j1, k0, k1, Ni, Nj, Nk, N0):

    j, k = cuda.grid(2)
    j += j0
    k += k0
    if (j < j1 + 1) and \
       (k < k1 + 1):
        n = (Ni * i) + (Nj * j) + (Nk * k) + N0
        m = (j - j0) * (k1 - k0 + 1) + (k - k0)
        buf[m] = h[n]

# (private) (kernel関数)
@cuda.jit(cache=True)
def _h2d_gpu(i, h, buf, j0, j1, k0, k1, Ni, Nj, Nk, N0):

    j, k = cuda.grid(2)
    j += j0
    k += k0
    if (j < j1 + 1) and \
       (k < k1 + 1):
        n = (Ni * i) + (Nj * j) + (Nk * k) + N0
        m = (j - j0) * (k1 - k0 + 1) + (k - k0)
        h[n] = buf[m]
