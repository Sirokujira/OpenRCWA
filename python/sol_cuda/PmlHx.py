# -*- coding: utf-8 -*-
"""
PmlHx.py (CUDA)
"""

import math
from numba import cuda

def calHx(block1d,
    Ny, Nz, Hx, Ey, Ez, Hxy, Hxz, RYc, RZc,
    iPmlHx, gPmlYc, gPmlZc, rPmlH, rPml, pml_l,
    Ni, Nj, Nk, N0):

    numHx = iPmlHx.shape[0]
    grid = math.ceil(numHx / block1d)

    _calcHx_gpu[grid, block1d](
        Ny, Nz, Hx, Ey, Ez, Hxy, Hxz, RYc, RZc,
        numHx, iPmlHx, gPmlYc, gPmlZc, rPmlH, rPml, pml_l,
        Ni, Nj, Nk, N0)

# (private)
@cuda.jit(cache=True)
def _calcHx_gpu(
    Ny, Nz, Hx, Ey, Ez, Hxy, Hxz, RYc, RZc,
    numHx, iPmlHx, gPmlYc, gPmlZc, rPmlH, rPml, pml_l,
    Ni, Nj, Nk, N0):

    num = cuda.grid(1)

    if num < numHx:
        i, j, k, m = iPmlHx[num]

        n   = Ni * (i    ) + Nj * (j    ) + Nk * (k    ) + N0
        nj1 = Ni * (i    ) + Nj * (j + 1) + Nk * (k    ) + N0
        nk1 = Ni * (i    ) + Nj * (j    ) + Nk * (k + 1) + N0

        jc = min(max(j, 0), Ny - 1)
        dez = Ez[nj1] - Ez[n]
        ry = RYc[jc] * rPmlH[m]
        Hxy[num] = (Hxy[num] - (ry * dez)) / (1 + (gPmlYc[j + pml_l] * rPml[m]))

        kc = min(max(k, 0), Nz - 1)
        dey = Ey[nk1] - Ey[n]
        rz = RZc[kc] * rPmlH[m]
        Hxz[num] = (Hxz[num] + (rz * dey)) / (1 + (gPmlZc[k + pml_l] * rPml[m]))

        Hx[n] = Hxy[num] + Hxz[num]
