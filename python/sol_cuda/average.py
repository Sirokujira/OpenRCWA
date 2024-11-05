# -*- coding: utf-8 -*-
"""
average.py (CUDA)
平均電磁界
"""

import math
#import numpy as np
from numba import cuda

# reduction sum
@cuda.reduce
def sum_reduce(a, b):
    return a + b

# E/H average
def calcA(ublock, average_array,
    Ex, Ey, Ez, Hx, Hy, Hz,
    iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0):

    area = (kMax - kMin, jMax - jMin, iMax - iMin)

    ugrid = (math.ceil(area[0] / ublock[0]),
             math.ceil(area[1] / ublock[1]),
             math.ceil(area[2] / ublock[2]))

    _calcE_gpu[ugrid, ublock](average_array,
        Ex, Ey, Ez,
        iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)
    esum = sum_reduce(average_array)

    _calcH_gpu[ugrid, ublock](average_array,
        Hx, Hy, Hz,
        iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)
    hsum = sum_reduce(average_array)

    return esum, hsum

@cuda.jit(cache=True)
def _calcE_gpu(average_array,
    Ex, Ey, Ez,
    iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0):

    k, j, i = cuda.grid(3)

    idx = i * (jMax - jMin) * (kMax - kMin) + j * (kMax - kMin) + k

    i += iMin
    j += jMin
    k += kMin
    if (i < iMax) and (j < jMax) and (k < kMax):
        n   = (Ni * (i    )) + (Nj * (j    )) + (Nk * (k    )) + N0
        ni  = (Ni * (i + 1)) + (Nj * (j    )) + (Nk * (k    )) + N0
        nj  = (Ni * (i    )) + (Nj * (j + 1)) + (Nk * (k    )) + N0
        nk  = (Ni * (i    )) + (Nj * (j    )) + (Nk * (k + 1)) + N0
        njk = (Ni * (i    )) + (Nj * (j + 1)) + (Nk * (k + 1)) + N0
        nki = (Ni * (i + 1)) + (Nj * (j    )) + (Nk * (k + 1)) + N0
        nij = (Ni * (i + 1)) + (Nj * (j + 1)) + (Nk * (k    )) + N0

        average_array[idx] = \
            abs(Ex[n] + Ex[nj] + Ex[nk] + Ex[njk]) + \
            abs(Ey[n] + Ey[nk] + Ey[ni] + Ey[nki]) + \
            abs(Ez[n] + Ez[ni] + Ez[nj] + Ez[nij])

@cuda.jit(cache=True)
def _calcH_gpu(average_array,
    Hx, Hy, Hz,
    iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0):

    k, j, i = cuda.grid(3)

    idx = i * (jMax - jMin)* (kMax - kMin) + j * (kMax - kMin) + k

    i += iMin
    j += jMin
    k += kMin
    if (i < iMax) and (j < jMax) and (k < kMax):
        n   = (Ni * (i    )) + (Nj * (j    )) + (Nk * (k    )) + N0
        ni  = (Ni * (i + 1)) + (Nj * (j    )) + (Nk * (k    )) + N0
        nj  = (Ni * (i    )) + (Nj * (j + 1)) + (Nk * (k    )) + N0
        nk  = (Ni * (i    )) + (Nj * (j    )) + (Nk * (k + 1)) + N0

        average_array[idx] = \
            abs(Hx[n] + Hx[ni]) + \
            abs(Hy[n] + Hy[nj]) + \
            abs(Hz[n] + Hz[nk])

"""
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

    return [sume / (4.0 * Nx * Ny * Nz), sumh / (2.0 * Nx * Ny * Nz)]
"""
"""
def calc(Nx, Ny, Nz, Ex, Ey, Ez, Hx, Hy, Hz):
    #return \
    #[(np.sum(np.abs(Ex)) + np.sum(np.abs(Ey)) + np.sum(np.abs(Ez))) / (Nx * Ny * Nz), \
    # (np.sum(np.abs(Hx)) + np.sum(np.abs(Hy)) + np.sum(np.abs(Hz))) / (Nx * Ny * Nz)]

    esum = np.sum(np.abs(Ex)) + np.sum(np.abs(Ey)) + np.sum(np.abs(Ez))
    hsum = np.sum(np.abs(Hx)) + np.sum(np.abs(Hy)) + np.sum(np.abs(Hz))
    return np.array([esum, hsum]) / (Nx * Ny *Nz)                                                        
"""