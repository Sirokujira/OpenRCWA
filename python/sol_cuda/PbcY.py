# -*- coding: utf-8 -*-
"""
PbcY.py (CUDA)
Y方向周期境界条件
"""

import math
from numba import cuda

def y(block2d, Hz, Hx, Ny, kMin, kMax, iMin, iMax, Ni, Nj, Nk, N0):

    gridHz = (math.ceil((kMax - kMin + 1) / block2d[0]),
              math.ceil((iMax - iMin + 2) / block2d[1]))
    _Hz_gpu[gridHz, block2d](Hz, Ny, kMin, kMax, iMin, iMax, Ni, Nj, Nk, N0)

    gridHx = (math.ceil((kMax - kMin + 2) / block2d[0]),
              math.ceil((iMax - iMin + 1) / block2d[1]))
    _Hx_gpu[gridHx, block2d](Hx, Ny, kMin, kMax, iMin, iMax, Ni, Nj, Nk, N0)

@cuda.jit(cache=True)
def _Hz_gpu(Hz, Ny, kMin, kMax, iMin, iMax, Ni, Nj, Nk, N0):

    k, i = cuda.grid(2)
    k += kMin - 0
    i += iMin - 1

    if (k <= kMax) and \
       (i <= iMax):
        id1 = -1
        id2 = 0
        id3 = Ny - 1
        id4 = Ny
        n1 = (Ni * i) + (Nj * id1) + (Nk * k) + N0
        n2 = (Ni * i) + (Nj * id2) + (Nk * k) + N0
        n3 = (Ni * i) + (Nj * id3) + (Nk * k) + N0
        n4 = (Ni * i) + (Nj * id4) + (Nk * k) + N0
        Hz[n1] = Hz[n3]
        Hz[n4] = Hz[n2]

@cuda.jit(cache=True)
def _Hx_gpu(Hx, Ny, kMin, kMax, iMin, iMax, Ni, Nj, Nk, N0):

    k, i = cuda.grid(2)
    k += kMin - 1
    i += iMin - 0

    if (k <= kMax) and \
       (i <= iMax):
        id1 = -1
        id2 = 0
        id3 = Ny - 1
        id4 = Ny
        n1 = (Ni * i) + (Nj * id1) + (Nk * k) + N0
        n2 = (Ni * i) + (Nj * id2) + (Nk * k) + N0
        n3 = (Ni * i) + (Nj * id3) + (Nk * k) + N0
        n4 = (Ni * i) + (Nj * id4) + (Nk * k) + N0
        Hx[n1] = Hx[n3]
        Hx[n4] = Hx[n2]
