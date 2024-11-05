# -*- coding: utf-8 -*-
"""
comm_Z.py (CUDA + MPI)
Z境界のHx/Hyを共有する
"""

import math
from numba import cuda
from mpi4py import MPI

def share(mode, d_Hx, d_Hy,
    SendBuf_Bz_hx, SendBuf_Bz_hy, RecvBuf_Bz_hx, RecvBuf_Bz_hy,
    d_SendBuf_Bz_hx, d_SendBuf_Bz_hy, d_RecvBuf_Bz_hx, d_RecvBuf_Bz_hy,
    Bz_ihx, Bz_ihy, Bz_jhx, Bz_jhy,
    Npx, Npy, Npz, Ipx, Ipy, Ipz,
    kMin, kMax, Ni, Nj, Nk, N0):

    assert((mode == 0) or (mode == 1))

    ihx0, ihx1 = Bz_ihx
    jhx0, jhx1 = Bz_jhx
    ihy0, ihy1 = Bz_ihy
    jhy0, jhy1 = Bz_jhy

    # grid, block
    block = (16, 16)
    grid_hx = (math.ceil((ihx1 - ihx0 + 1) / block[0]),
               math.ceil((jhx1 - jhx0 + 1) / block[1]))
    grid_hy = (math.ceil((ihy1 - ihy0 + 1) / block[0]),
               math.ceil((jhy1 - jhy0 + 1) / block[1]))

    bz = [[Ipz > 0, Ipz < Npz - 1], [Ipz == 0, Ipz == Npz - 1]]
    pz = [[Ipz - 1, Ipz + 1], [Npz - 1, 0]]
    ksend = [kMin + 0, kMax - 1]
    krecv = [kMin - 1, kMax + 0]

    for side in range(2):
        if bz[mode][side]:
            # from device memory to host buffer
            k = ksend[side]
            _d2h_gpu[grid_hx, block](k, d_Hx, d_SendBuf_Bz_hx, ihx0, ihx1, jhx0, jhx1, Ni, Nj, Nk, N0)
            _d2h_gpu[grid_hy, block](k, d_Hy, d_SendBuf_Bz_hy, ihy0, ihy1, jhy0, jhy1, Ni, Nj, Nk, N0)
            d_SendBuf_Bz_hx.copy_to_host(SendBuf_Bz_hx)
            d_SendBuf_Bz_hy.copy_to_host(SendBuf_Bz_hy)

            # MPI
            ipz = pz[mode][side]
            dst = (Ipx * Npy * Npz) + (Ipy * Npz) + ipz
            MPI.COMM_WORLD.Sendrecv(SendBuf_Bz_hx, dst, recvbuf=RecvBuf_Bz_hx, source=dst)
            MPI.COMM_WORLD.Sendrecv(SendBuf_Bz_hy, dst, recvbuf=RecvBuf_Bz_hy, source=dst)

            # from host buffer to device memory
            k = krecv[side]
            d_RecvBuf_Bz_hx = cuda.to_device(RecvBuf_Bz_hx)
            d_RecvBuf_Bz_hy = cuda.to_device(RecvBuf_Bz_hy)
            _h2d_gpu[grid_hx, block](k, d_Hx, d_RecvBuf_Bz_hx, ihx0, ihx1, jhx0, jhx1, Ni, Nj, Nk, N0)
            _h2d_gpu[grid_hy, block](k, d_Hy, d_RecvBuf_Bz_hy, ihy0, ihy1, jhy0, jhy1, Ni, Nj, Nk, N0)

# (private) (kernel関数)
@cuda.jit(cache=True)
def _d2h_gpu(k, h, buf, i0, i1, j0, j1, Ni, Nj, Nk, N0):

    i, j = cuda.grid(2)
    i += i0
    j += j0
    if (i < i1 + 1) and \
       (j < j1 + 1):
        n = (Ni * i) + (Nj * j) + (Nk * k) + N0
        m = (i - i0) * (j1 - j0 + 1) + (j - j0)
        buf[m] = h[n]

# (private) (kernel関数)
@cuda.jit(cache=True)
def _h2d_gpu(k, h, buf, i0, i1, j0, j1, Ni, Nj, Nk, N0):

    i, j = cuda.grid(2)
    i += i0
    j += j0
    if (i < i1 + 1) and \
       (j < j1 + 1):
        n = (Ni * i) + (Nj * j) + (Nk * k) + N0
        m = (i - i0) * (j1 - j0 + 1) + (j - j0)
        h[n] = buf[m]
