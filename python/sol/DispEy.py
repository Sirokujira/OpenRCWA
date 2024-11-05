# -*- coding: utf-8 -*-
"""
DispEy.py
分散性媒質Ey
"""

from numba import jit, prange
import sol.planewave

@jit(cache=True, nopython=True, nogil=True, parallel=True)
def calEy(
    source, t, fPlanewave, Yc, Zn, Xn, Ey, iDispEy, fDispEy,
    Ni, Nj, Nk, N0):

    for n in prange(iDispEy.shape[0]):
        i  = iDispEy[n, 0]
        j  = iDispEy[n, 1]
        k  = iDispEy[n, 2]
        f1 = fDispEy[n, 1]
        f2 = fDispEy[n, 2]
        f3 = fDispEy[n, 3]
        n0 = (Ni * i) + (Nj * j) + (Nk * k) + N0

        fi = 0
        if source == 1:
            fi, _ = sol.planewave.f(Xn[i], Yc[j], Zn[k], t, fPlanewave[1], fPlanewave)

        Ey[n0] += f1 * fDispEy[n, 0]

        fDispEy[n, 0] = f2 * (Ey[n0] + fi) + f3 * fDispEy[n, 0]
