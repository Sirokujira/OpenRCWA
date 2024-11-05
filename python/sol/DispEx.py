# -*- coding: utf-8 -*-
"""
DispEx.py
分散性媒質Ex
"""

from numba import jit, prange
import sol.planewave

@jit(cache=True, nopython=True, nogil=True, parallel=True)
def calEx(
    source, t, fPlanewave, Xc, Yn, Zn, Ex, iDispEx, fDispEx,
    Ni, Nj, Nk, N0):

    for n in prange(iDispEx.shape[0]):
        i  = iDispEx[n, 0]
        j  = iDispEx[n, 1]
        k  = iDispEx[n, 2]
        f1 = fDispEx[n, 1]
        f2 = fDispEx[n, 2]
        f3 = fDispEx[n, 3]
        n0 = (Ni * i) + (Nj * j) + (Nk * k) + N0

        fi = 0
        if source == 1:
            fi, _ = sol.planewave.f(Xc[i], Yn[j], Zn[k], t, fPlanewave[0], fPlanewave)

        Ex[n0] += f1 * fDispEx[n, 0]

        fDispEx[n, 0] = f2 * (Ex[n0] + fi) + f3 * fDispEx[n, 0]
