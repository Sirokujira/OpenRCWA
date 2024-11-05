# -*- coding: utf-8 -*-
"""
PmlEz.py (CUDA)
"""

import math
from numba import cuda

def calEz(block1d,
    Nx, Ny, Ez, Hx, Hy, Ezx, Ezy, RXn, RYn,
    iPmlEz, gPmlXn, gPmlYn, rPmlE, rPml, pml_l,
    Ni, Nj, Nk, N0):

    numEz = iPmlEz.shape[0]
    grid = math.ceil(numEz / block1d)

    _calcEz_gpu[grid, block1d](
        Nx, Ny, Ez, Hx, Hy, Ezx, Ezy, RXn, RYn,
        numEz, iPmlEz, gPmlXn, gPmlYn, rPmlE, rPml, pml_l,
        Ni, Nj, Nk, N0)

# (private)
@cuda.jit(cache=True)
def _calcEz_gpu(
    Nx, Ny, Ez, Hx, Hy, Ezx, Ezy, RXn, RYn,
    numEz, iPmlEz, gPmlXn, gPmlYn, rPmlE, rPml, pml_l,
    Ni, Nj, Nk, N0):

    num = cuda.grid(1)

    if num < numEz:
        i, j, k, m = iPmlEz[num]

        n   = Ni * (i    ) + Nj * (j    ) + Nk * (k    ) + N0
        ni1 = Ni * (i - 1) + Nj * (j    ) + Nk * (k    ) + N0
        nj1 = Ni * (i    ) + Nj * (j - 1) + Nk * (k    ) + N0

        ic = min(max(i, 0), Nx    )
        dhy = Hy[n] - Hy[ni1]
        rx = RXn[ic] * rPmlE[m]
        Ezx[num] = (Ezx[num] + (rx * dhy)) / (1 + (gPmlXn[i + pml_l] * rPml[m]))

        jc = min(max(j, 0), Ny    )
        dhx = Hx[n] - Hx[nj1]
        ry = RYn[jc] * rPmlE[m]
        Ezy[num] = (Ezy[num] - (ry * dhx)) / (1 + (gPmlYn[j + pml_l] * rPml[m]))

        Ez[n] = Ezx[num] + Ezy[num]
