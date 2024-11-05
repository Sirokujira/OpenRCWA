# -*- coding: utf-8 -*-
"""
DispEz.py
分散性媒質Ez
"""

from numba import jit, prange
import sol.planewave

@jit(cache=True, nopython=True, nogil=True, parallel=True)
def calEz(
    source, t, fPlanewave, Zc, Xn, Yn, Ez, iDispEz, fDispEz,
    Ni, Nj, Nk, N0):

    for n in prange(iDispEz.shape[0]):
        i  = iDispEz[n, 0]
        j  = iDispEz[n, 1]
        k  = iDispEz[n, 2]
        f1 = fDispEz[n, 1]
        f2 = fDispEz[n, 2]
        f3 = fDispEz[n, 3]
        n0 = (Ni * i) + (Nj * j) + (Nk * k) + N0

        fi = 0
        if source == 1:
            fi, _ = sol.planewave.f(Xn[i], Yn[j], Zc[k], t, fPlanewave[2], fPlanewave)

        Ez[n0] += f1 * fDispEz[n, 0]

        fDispEz[n, 0] = f2 * (Ez[n0] + fi) + f3 * fDispEz[n, 0]
