# -*- coding: utf-8 -*-
"""
PbcX.py (CUDA)
X方向周期境界条件
"""

import math
from numba import cuda

def x(block2d, Hy, Hz, Nx, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0):

    gridHy = (math.ceil((jMax - jMin + 1) / block2d[0]),
              math.ceil((kMax - kMin + 2) / block2d[1]))
    _Hy_gpu[gridHy, block2d](Hy, Nx, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)

    gridHz = (math.ceil((jMax - jMin + 2) / block2d[0]),
              math.ceil((kMax - kMin + 1) / block2d[1]))
    _Hz_gpu[gridHz, block2d](Hz, Nx, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)

@cuda.jit(cache=True)
def _Hy_gpu(Hy, Nx, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0):

    j, k = cuda.grid(2)
    j += jMin - 0
    k += kMin - 1

    if (j <= jMax) and \
       (k <= kMax):
        id1 = -1
        id2 = 0
        id3 = Nx - 1
        id4 = Nx
        n1 = (Ni * id1) + (Nj * j) + (Nk * k) + N0
        n2 = (Ni * id2) + (Nj * j) + (Nk * k) + N0
        n3 = (Ni * id3) + (Nj * j) + (Nk * k) + N0
        n4 = (Ni * id4) + (Nj * j) + (Nk * k) + N0
        Hy[n1] = Hy[n3]
        Hy[n4] = Hy[n2]

@cuda.jit(cache=True)
def _Hz_gpu(Hz, Nx, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0):

    j, k = cuda.grid(2)
    j += jMin - 1
    k += kMin - 0

    if (j <= jMax) and \
       (k <= kMax):
        id1 = -1
        id2 = 0
        id3 = Nx - 1
        id4 = Nx
        n1 = (Ni * id1) + (Nj * j) + (Nk * k) + N0
        n2 = (Ni * id2) + (Nj * j) + (Nk * k) + N0
        n3 = (Ni * id3) + (Nj * j) + (Nk * k) + N0
        n4 = (Ni * id4) + (Nj * j) + (Nk * k) + N0
        Hz[n1] = Hz[n3]
        Hz[n4] = Hz[n2]
