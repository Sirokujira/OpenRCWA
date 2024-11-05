# -*- coding: utf-8 -*-
"""
DispEy.py (CUDA)
分散性媒質Ey
"""

import math
from numba import cuda
import sol_cuda.planewave

def calEy(block1d,
    source, t, fPlanewave, Yc, Zn, Xn, Ey, iDispEy, fDispEy,
    Ni, Nj, Nk, N0):

    ndisp = iDispEy.shape[0]

    if ndisp == 0:
        return

    grid = math.ceil(ndisp / block1d)

    _calEy_gpu[grid, block1d](ndisp,
        source, t, fPlanewave, Yc, Zn, Xn, Ey, iDispEy, fDispEy,
        Ni, Nj, Nk, N0)

@cuda.jit(cache=True)
def _calEy_gpu(ndisp,
    source, t, fPlanewave, Yc, Zn, Xn, Ey, iDispEy, fDispEy,
    Ni, Nj, Nk, N0):

    idisp = cuda.grid(1)

    if idisp < ndisp:
        i, j, k = iDispEy[idisp, 0:3]
        f1, f2, f3 = fDispEy[idisp, 1:4]

        n = (Ni * i) + (Nj * j) + (Nk * k) + N0

        fi = 0
        if source == 1:
            fi, _ = sol_cuda.planewave.f(Xn[i], Yc[j], Zn[k], t, fPlanewave[1], fPlanewave)

        Ey[n] += f1 * fDispEy[idisp, 0]

        fDispEy[idisp, 0] = f2 * (Ey[n] + fi) + f3 * fDispEy[idisp, 0]
