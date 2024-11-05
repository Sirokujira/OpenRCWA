# -*- coding: utf-8 -*-
"""
PmlEy.py (CUDA)
"""

import math
from numba import cuda

def calEy(block1d,
    Nz, Nx, Ey, Hz, Hx, Eyz, Eyx, RZn, RXn,
    iPmlEy, gPmlZn, gPmlXn, rPmlE, rPml, pml_l,
    Ni, Nj, Nk, N0):

    numEy = iPmlEy.shape[0]
    grid = math.ceil(numEy / block1d)

    _calcEy_gpu[grid, block1d](
        Nz, Nx, Ey, Hz, Hx, Eyz, Eyx, RZn, RXn,
        numEy, iPmlEy, gPmlZn, gPmlXn, rPmlE, rPml, pml_l,
        Ni, Nj, Nk, N0)

# (private)
@cuda.jit(cache=True)
def _calcEy_gpu(
    Nz, Nx, Ey, Hz, Hx, Eyz, Eyx, RZn, RXn,
    numEy, iPmlEy, gPmlZn, gPmlXn, rPmlE, rPml, pml_l,
    Ni, Nj, Nk, N0):

    num = cuda.grid(1)

    if num < numEy:
        i, j, k, m = iPmlEy[num]

        n   = Ni * (i    ) + Nj * (j    ) + Nk * (k    ) + N0
        nk1 = Ni * (i    ) + Nj * (j    ) + Nk * (k - 1) + N0
        ni1 = Ni * (i - 1) + Nj * (j    ) + Nk * (k    ) + N0

        kc = min(max(k, 0), Nz    )
        dhx = Hx[n] - Hx[nk1]
        rz = RZn[kc] * rPmlE[m]
        Eyz[num] = (Eyz[num] + (rz * dhx)) / (1 + (gPmlZn[k + pml_l] * rPml[m]))

        ic = min(max(i, 0), Nx    )
        dhz = Hz[n] - Hz[ni1]
        rx = RXn[ic] * rPmlE[m]
        Eyx[num] = (Eyx[num] - (rx * dhz)) / (1 + (gPmlXn[i + pml_l] * rPml[m]))

        Ey[n] = Eyz[num] + Eyx[num]
