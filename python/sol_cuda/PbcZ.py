# -*- coding: utf-8 -*-
"""
PbcZ.py (CUDA)
Z方向周期境界条件
"""

import math
from numba import cuda

def z(block2d, Hx, Hy, Nz, iMin, iMax, jMin, jMax, Ni, Nj, Nk, N0):

    gridHx = (math.ceil((iMax - iMin + 1) / block2d[0]),
              math.ceil((jMax - jMin + 2) / block2d[1]))
    _Hx_gpu[gridHx, block2d](Hx, Nz, iMin, iMax, jMin, jMax, Ni, Nj, Nk, N0)

    gridHy = (math.ceil((iMax - iMin + 2) / block2d[0]),
              math.ceil((jMax - jMin + 1) / block2d[1]))
    _Hy_gpu[gridHy, block2d](Hy, Nz, iMin, iMax, jMin, jMax, Ni, Nj, Nk, N0)

@cuda.jit(cache=True)
def _Hx_gpu(Hx, Nz, iMin, iMax, jMin, jMax, Ni, Nj, Nk, N0):

    i, j = cuda.grid(2)
    i += iMin - 0
    j += jMin - 1

    if (i <= iMax) and \
       (j <= jMax):
        id1 = -1
        id2 = 0
        id3 = Nz - 1
        id4 = Nz
        n1 = (Ni * i) + (Nj * j) + (Nk * id1) + N0
        n2 = (Ni * i) + (Nj * j) + (Nk * id2) + N0
        n3 = (Ni * i) + (Nj * j) + (Nk * id3) + N0
        n4 = (Ni * i) + (Nj * j) + (Nk * id4) + N0
        Hx[n1] = Hx[n3]
        Hx[n4] = Hx[n2]

@cuda.jit(cache=True)
def _Hy_gpu(Hy, Nz, iMin, iMax, jMin, jMax, Ni, Nj, Nk, N0):

    i, j = cuda.grid(2)
    i += iMin - 1
    j += jMin - 0

    if (i <= iMax) and \
       (j <= jMax):
        id1 = -1
        id2 = 0
        id3 = Nz - 1
        id4 = Nz
        n1 = (Ni * i) + (Nj * j) + (Nk * id1) + N0
        n2 = (Ni * i) + (Nj * j) + (Nk * id2) + N0
        n3 = (Ni * i) + (Nj * j) + (Nk * id3) + N0
        n4 = (Ni * i) + (Nj * j) + (Nk * id4) + N0
        Hy[n1] = Hy[n3]
        Hy[n4] = Hy[n2]
