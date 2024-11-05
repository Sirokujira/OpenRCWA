# -*- coding: utf-8 -*-
"""
DispEz.py (CUDA)
分散性媒質Ez
"""

import math
from numba import cuda
import sol_cuda.planewave

def calEz(block1d,
    source, t, fPlanewave, Zc, Xn, Yn, Ez, iDispEz, fDispEz,
    Ni, Nj, Nk, N0):

    ndisp = iDispEz.shape[0]

    if ndisp == 0:
        return

    grid = math.ceil(ndisp / block1d)

    _calEz_gpu[grid, block1d](ndisp,
        source, t, fPlanewave, Zc, Xn, Yn, Ez, iDispEz, fDispEz,
        Ni, Nj, Nk, N0)

@cuda.jit(cache=True)
def _calEz_gpu(ndisp,
    source, t, fPlanewave, Zc, Xn, Yn, Ez, iDispEz, fDispEz,
    Ni, Nj, Nk, N0):

    idisp = cuda.grid(1)

    if idisp < ndisp:
        i, j, k = iDispEz[idisp, 0:3]
        f1, f2, f3 = fDispEz[idisp, 1:4]

        n = (Ni * i) + (Nj * j) + (Nk * k) + N0

        fi = 0
        if source == 1:
            fi, _ = sol_cuda.planewave.f(Xn[i], Yn[j], Zc[k], t, fPlanewave[2], fPlanewave)

        Ez[n] += f1 * fDispEz[idisp, 0]

        fDispEz[idisp, 0] = f2 * (Ez[n] + fi) + f3 * fDispEz[idisp, 0]
