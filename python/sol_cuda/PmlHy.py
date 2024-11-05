# -*- coding: utf-8 -*-
"""
PmlHy.py (CUDA)
"""

import math
from numba import cuda

def calHy(block1d,
    Nz, Nx, Hy, Ez, Ex, Hyz, Hyx, RZc, RXc,
    iPmlHy, gPmlZc, gPmlXc, rPmlH, rPml, pml_l,
    Ni, Nj, Nk, N0):

    numHy = iPmlHy.shape[0]
    grid = math.ceil(numHy / block1d)

    _calcHy_gpu[grid, block1d](
        Nz, Nx, Hy, Ez, Ex, Hyz, Hyx, RZc, RXc,
        numHy, iPmlHy, gPmlZc, gPmlXc, rPmlH, rPml, pml_l,
        Ni, Nj, Nk, N0)

# (private)
@cuda.jit(cache=True)
def _calcHy_gpu(
    Nz, Nx, Hy, Ez, Ex, Hyz, Hyx, RZc, RXc,
    numHy, iPmlHy, gPmlZc, gPmlXc, rPmlH, rPml, pml_l,
    Ni, Nj, Nk, N0):

    num = cuda.grid(1)

    if num < numHy:
        i, j, k, m = iPmlHy[num]

        n   = Ni * (i    ) + Nj * (j    ) + Nk * (k    ) + N0
        nk1 = Ni * (i    ) + Nj * (j    ) + Nk * (k + 1) + N0
        ni1 = Ni * (i + 1) + Nj * (j    ) + Nk * (k    ) + N0

        kc = min(max(k, 0), Nz - 1)
        dex = Ex[nk1] - Ex[n]
        rz = RZc[kc] * rPmlH[m]
        Hyz[num] = (Hyz[num] - (rz * dex)) / (1 + (gPmlZc[k + pml_l] * rPml[m]))

        ic = min(max(i, 0), Nx - 1)
        dez = Ez[ni1] - Ez[n]
        rx = RXc[ic] * rPmlH[m]
        Hyx[num] = (Hyx[num] + (rx * dez)) / (1 + (gPmlXc[i + pml_l] * rPml[m]))

        Hy[n] = Hyz[num] + Hyx[num]
