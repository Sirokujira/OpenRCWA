# -*- coding: utf-8 -*-
"""
updateEy.py (CUDA)
"""

import math
from numba import cuda
import sol_cuda.planewave

# Ey更新
def calEy(block3d,
    Parm, t, fPlanewave, Yc, Zn, Xn, Ey, Hz, Hx, iEy,
    C1E, C2E, K1Ey, K2Ey, RZn, RXn,
    iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0):

    grid = (math.ceil((kMax - kMin + 1) / block3d[0]),
            math.ceil((jMax - jMin + 0) / block3d[1]),
            math.ceil((iMax - iMin + 1) / block3d[2]))

    if   Parm['source'] == 0:
        # 給電点
        if Parm['vector']:
            # (1)
            _Ey_f_vector[grid, block3d](
                Ey, Hz, Hx, K1Ey, K2Ey, RZn, RXn,
                iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)
        else:
            # (2)
            _Ey_f_no_vector[grid, block3d](
                Ey, Hz, Hx, iEy, C1E, C2E, RZn, RXn,
                iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)
    elif Parm['source'] == 1:
        # 平面波入射
        if Parm['vector']:
            # (3)
            _Ey_p_vector[grid, block3d](
                t, fPlanewave, Yc, Zn, Xn,
                Ey, Hz, Hx, iEy, K1Ey, K2Ey, RZn, RXn,
                iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)
        else:
            # (4)
            _Ey_p_no_vector[grid, block3d](
                t, fPlanewave, Yc, Zn, Xn,
                Ey, Hz, Hx, iEy, C1E, C2E, RZn, RXn,
                iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)

# 給電点
# (1)
@cuda.jit(cache=True)
def _Ey_f_vector(
    Ey, Hz, Hx, K1Ey, K2Ey, RZn, RXn,
    iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0):

    k, j, i = cuda.grid(3)
    i += iMin
    j += jMin
    k += kMin

    if (i < iMax + 1) and \
       (j < jMax + 0) and \
       (k < kMax + 1):
        n = (Ni * i) + (Nj * j) + (Nk * k) + N0
        Ey[n] = K1Ey[n] * Ey[n] \
              + K2Ey[n] * (RZn[k] * (Hx[n] - Hx[n - Nk])
                         - RXn[i] * (Hz[n] - Hz[n - Ni]))

# (2)
@cuda.jit(cache=True)
def _Ey_f_no_vector(
    Ey, Hz, Hx, iEy, C1E, C2E, RZn, RXn,
    iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0):

    k, j, i = cuda.grid(3)
    i += iMin
    j += jMin
    k += kMin

    if (i < iMax + 1) and \
       (j < jMax + 0) and \
       (k < kMax + 1):
        n = (Ni * i) + (Nj * j) + (Nk * k) + N0
        m = iEy[n]
        Ey[n] = C1E[m] * Ey[n] \
              + C2E[m] * (RZn[k] * (Hx[n] - Hx[n - Nk])
                        - RXn[i] * (Hz[n] - Hz[n - Ni]))

# 平面波入射
# (3)
@cuda.jit(cache=True)
def _Ey_p_vector(
    t, fPlanewave, Yc, Zn, Xn,
    Ey, Hz, Hx, iEy, K1Ey, K2Ey, RZn, RXn,
    iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0):

    k, j, i = cuda.grid(3)
    i += iMin
    j += jMin
    k += kMin

    if (i < iMax + 1) and \
       (j < jMax + 0) and \
       (k < kMax + 1):
        n = (Ni * i) + (Nj * j) + (Nk * k) + N0
        m = iEy[n]
        if m == 0:
            Ey[n] += \
                  + RZn[k] * (Hx[n] - Hx[n - Nk]) \
                  - RXn[i] * (Hz[n] - Hz[n - Ni])
        else:
            fi, dfi = sol_cuda.planewave.f(Xn[i], Yc[j], Zn[k], t, fPlanewave[1], fPlanewave)
            Ey[n] = K1Ey[n] * Ey[n] \
                  + K2Ey[n] * (RZn[k] * (Hx[n] - Hx[n - Nk])
                             - RXn[i] * (Hz[n] - Hz[n - Ni])) \
                  - (K1Ey[n] - K2Ey[n]) * dfi \
                  - (1 - K1Ey[n]) * fi

# (4)
@cuda.jit(cache=True)
def _Ey_p_no_vector(
    t, fPlanewave, Yc, Zn, Xn,
    Ey, Hz, Hx, iEy, C1E, C2E, RZn, RXn,
    iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0):

    k, j, i = cuda.grid(3)
    i += iMin
    j += jMin
    k += kMin

    if (i < iMax + 1) and \
       (j < jMax + 0) and \
       (k < kMax + 1):
        n = (Ni * i) + (Nj * j) + (Nk * k) + N0
        m = iEy[n]
        if m == 0:
            Ey[n] += \
                  + RZn[k] * (Hx[n] - Hx[n - Nk]) \
                  - RXn[i] * (Hz[n] - Hz[n - Ni])
        else:
            fi, dfi = sol_cuda.planewave.f(Xn[i], Yc[j], Zn[k], t, fPlanewave[1], fPlanewave)
            Ey[n] = C1E[m] * Ey[n] \
                  + C2E[m] * (RZn[k] * (Hx[n] - Hx[n - Nk])
                            - RXn[i] * (Hz[n] - Hz[n - Ni])) \
                  - (C1E[m] - C2E[m]) * dfi \
                  - (1 - C1E[m]) * fi
