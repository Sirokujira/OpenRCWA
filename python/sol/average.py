# -*- coding: utf-8 -*-
"""
average.py
平均電磁界
"""

from numba import jit, prange

@jit(cache=True, nogil=True, parallel=True, nopython=True)
def calcA(Ex, Ey, Ez, Hx, Hy, Hz,
    iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0):

    sume = 0
    sumh = 0
    for i in prange(iMin, iMax):
        for j in range(jMin, jMax):
            for k in range(kMin, kMax):
                n   = (Ni * (i    )) + (Nj * (j    )) + (Nk * (k    )) + N0
                ni  = (Ni * (i + 1)) + (Nj * (j    )) + (Nk * (k    )) + N0
                nj  = (Ni * (i    )) + (Nj * (j + 1)) + (Nk * (k    )) + N0
                nk  = (Ni * (i    )) + (Nj * (j    )) + (Nk * (k + 1)) + N0
                njk = (Ni * (i    )) + (Nj * (j + 1)) + (Nk * (k + 1)) + N0
                nki = (Ni * (i + 1)) + (Nj * (j    )) + (Nk * (k + 1)) + N0
                nij = (Ni * (i + 1)) + (Nj * (j + 1)) + (Nk * (k    )) + N0
                sume += abs(Ex[n] + Ex[nj] + Ex[nk] + Ex[njk]) \
                      + abs(Ey[n] + Ey[nk] + Ey[ni] + Ey[nki]) \
                      + abs(Ez[n] + Ez[ni] + Ez[nj] + Ez[nij])
                sumh += abs(Hx[n] + Hx[ni]) \
                      + abs(Hy[n] + Hy[nj]) \
                      + abs(Hz[n] + Hz[nk])

    #return sume / (4.0 * Nx * Ny * Nz), sumh / (2.0 * Nx * Ny * Nz)
    return sume, sumh
