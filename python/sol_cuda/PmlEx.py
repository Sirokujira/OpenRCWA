# -*- coding: utf-8 -*-
"""
PmlEx.py (CUDA)
"""

import math
from numba import cuda

def calEx(block1d,
    Ny, Nz, Ex, Hy, Hz, Exy, Exz, RYn, RZn,
    iPmlEx, gPmlYn, gPmlZn, rPmlE, rPml, pml_l,
    Ni, Nj, Nk, N0):

    numEx = iPmlEx.shape[0]
    grid = math.ceil(numEx / block1d)

    _calcEx_gpu[grid, block1d](
        Ny, Nz, Ex, Hy, Hz, Exy, Exz, RYn, RZn,
        numEx, iPmlEx, gPmlYn, gPmlZn, rPmlE, rPml, pml_l,
        Ni, Nj, Nk, N0)

# (private)
@cuda.jit(cache=True)
def _calcEx_gpu(
    Ny, Nz, Ex, Hy, Hz, Exy, Exz, RYn, RZn,
    numEx, iPmlEx, gPmlYn, gPmlZn, rPmlE, rPml, pml_l,
    Ni, Nj, Nk, N0):

    num = cuda.grid(1)

    if num < numEx:
        i, j, k, m = iPmlEx[num]

        n   = Ni * (i    ) + Nj * (j    ) + Nk * (k    ) + N0
        nj1 = Ni * (i    ) + Nj * (j - 1) + Nk * (k    ) + N0
        nk1 = Ni * (i    ) + Nj * (j    ) + Nk * (k - 1) + N0

        jc = min(max(j, 0), Ny    )
        dhz = Hz[n] - Hz[nj1]
        ry = RYn[jc] * rPmlE[m]
        Exy[num] = (Exy[num] + (ry * dhz)) / (1 + (gPmlYn[j + pml_l] * rPml[m]))

        kc = min(max(k, 0), Nz    )
        dhy = Hy[n] - Hy[nk1]
        rz = RZn[kc] * rPmlE[m]
        Exz[num] = (Exz[num] - (rz * dhy)) / (1 + (gPmlZn[k + pml_l] * rPml[m]))

        Ex[n] = Exy[num] + Exz[num]
