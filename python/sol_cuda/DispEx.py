# -*- coding: utf-8 -*-
"""
DispEx.py (CUDA)
分散性媒質Ex
"""

import math
from numba import cuda
import sol_cuda.planewave

def calEx(block1d,
    source, t, fPlanewave, Xc, Yn, Zn, Ex, iDispEx, fDispEx,
    Ni, Nj, Nk, N0):

    ndisp = iDispEx.shape[0]

    if ndisp == 0:
        return

    grid = math.ceil(ndisp / block1d)

    _calEx_gpu[grid, block1d](ndisp,
        source, t, fPlanewave, Xc, Yn, Zn, Ex, iDispEx, fDispEx,
        Ni, Nj, Nk, N0)

@cuda.jit(cache=True)
def _calEx_gpu(ndisp,
    source, t, fPlanewave, Xc, Yn, Zn, Ex, iDispEx, fDispEx,
    Ni, Nj, Nk, N0):

    idisp = cuda.grid(1)

    if idisp < ndisp:
        i, j, k = iDispEx[idisp, 0:3]
        f1, f2, f3 = fDispEx[idisp, 1:4]

        n = (Ni * i) + (Nj * j) + (Nk * k) + N0

        fi = 0
        if source == 1:
            fi, _ = sol_cuda.planewave.f(Xc[i], Yn[j], Zn[k], t, fPlanewave[0], fPlanewave)

        Ex[n] += f1 * fDispEx[idisp, 0]

        fDispEx[idisp, 0] = f2 * (Ex[n] + fi) + f3 * fDispEx[idisp, 0]
