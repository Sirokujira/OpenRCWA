# -*- coding: utf-8 -*-
"""
PmlHz.py (CUDA)
"""

import math
from numba import cuda

def calHz(block1d,
    Nx, Ny, Hz, Ex, Ey, Hzx, Hzy, RXc, RYc,
    iPmlHz, gPmlXc, gPmlYc, rPmlH, rPml, pml_l,
    Ni, Nj, Nk, N0):

    numHz = iPmlHz.shape[0]
    grid = math.ceil(numHz / block1d)

    _calcHz_gpu[grid, block1d](
        Nx, Ny, Hz, Ex, Ey, Hzx, Hzy, RXc, RYc,
        numHz, iPmlHz, gPmlXc, gPmlYc, rPmlH, rPml, pml_l,
        Ni, Nj, Nk, N0)

# (private)
@cuda.jit(cache=True)
def _calcHz_gpu(
    Nx, Ny, Hz, Ex, Ey, Hzx, Hzy, RXc, RYc,
    numHz, iPmlHz, gPmlXc, gPmlYc, rPmlH, rPml, pml_l,
    Ni, Nj, Nk, N0):

    num = cuda.grid(1)

    if num < numHz:
        i, j, k, m = iPmlHz[num]

        n   = Ni * (i    ) + Nj * (j    ) + Nk * (k    ) + N0
        ni1 = Ni * (i + 1) + Nj * (j    ) + Nk * (k    ) + N0
        nj1 = Ni * (i    ) + Nj * (j + 1) + Nk * (k    ) + N0

        ic = min(max(i, 0), Nx - 1)
        dey = Ey[ni1] - Ey[n]
        rx = RXc[ic] * rPmlH[m]
        Hzx[num] = (Hzx[num] - (rx * dey)) / (1 + (gPmlXc[i + pml_l] * rPml[m]))

        jc = min(max(j, 0), Ny - 1)
        dex = Ex[nj1] - Ex[n]
        ry = RYc[jc] * rPmlH[m]
        Hzy[num] = (Hzy[num] + (ry * dex)) / (1 + (gPmlYc[j + pml_l] * rPml[m]))

        Hz[n] = Hzx[num] + Hzy[num]
