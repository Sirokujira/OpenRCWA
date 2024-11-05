# -*- coding: utf-8 -*-
"""
updateHy.py (CUDA)
"""

import math
from numba import cuda
import sol_cuda.planewave

# Hy更新
def calHy(block3d,
    Parm, t, fPlanewave, Yn, Zc, Xc, Hy, Ez, Ex, iHy,
    C1H, C2H, K1Hy, K2Hy, RZc, RXc,
    iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0):

    grid = (math.ceil((kMax - kMin + 0) / block3d[0]),
            math.ceil((jMax - jMin + 1) / block3d[1]),
            math.ceil((iMax - iMin + 0) / block3d[2]))

    if   Parm['source'] == 0:
        # 給電点
        if Parm['vector']:
            # (1)
            _Hy_f_vector[grid, block3d](
                Hy, Ez, Ex, K1Hy, K2Hy, RZc, RXc,
                iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)
        else:
            # (2)
            _Hy_f_no_vector[grid, block3d](
                Hy, Ez, Ex, iHy, C1H, C2H, RZc, RXc,
                iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)
    elif Parm['source'] == 1:
        # 平面波入射
        if Parm['vector']:
            # (3)
            _Hy_p_vector[grid, block3d](
                t, fPlanewave, Yn, Zc, Xc,
                Hy, Ez, Ex, iHy, K1Hy, K2Hy, RZc, RXc,
                iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)
        else:
            # (4)
            _Hy_p_no_vector[grid, block3d](
                t, fPlanewave, Yn, Zc, Xc,
                Hy, Ez, Ex, iHy, C1H, C2H, RZc, RXc,
                iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)

# 給電点
# (1)
@cuda.jit(cache=True)
def _Hy_f_vector(
    Hy, Ez, Ex, K1Hy, K2Hy, RZc, RXc,
    iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0):

    k, j, i = cuda.grid(3)
    i += iMin
    j += jMin
    k += kMin

    if (i < iMax + 0) and \
       (j < jMax + 1) and \
       (k < kMax + 0):
        n = (Ni * i) + (Nj * j) + (Nk * k) + N0
        Hy[n] = K1Hy[n] * Hy[n] \
              - K2Hy[n] * (RZc[k] * (Ex[n + Nk] - Ex[n])
                         - RXc[i] * (Ez[n + Ni] - Ez[n]))

# (2)
@cuda.jit(cache=True)
def _Hy_f_no_vector(
    Hy, Ez, Ex, iHy, C1H, C2H, RZc, RXc,
    iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0):

    k, j, i = cuda.grid(3)
    i += iMin
    j += jMin
    k += kMin

    if (i < iMax + 0) and \
       (j < jMax + 1) and \
       (k < kMax + 0):
        n = (Ni * i) + (Nj * j) + (Nk * k) + N0
        m = iHy[n]
        Hy[n] = C1H[m] * Hy[n] \
              - C2H[m] * (RZc[k] * (Ex[n + Nk] - Ex[n])
                        - RXc[i] * (Ez[n + Ni] - Ez[n]))

# 平面波入射
# (3)
@cuda.jit(cache=True)
def _Hy_p_vector(
    t, fPlanewave, Yn, Zc, Xc,
    Hy, Ez, Ex, iHy, K1Hy, K2Hy, RZc, RXc,
    iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0):

    k, j, i = cuda.grid(3)
    i += iMin
    j += jMin
    k += kMin

    if (i < iMax + 0) and \
       (j < jMax + 1) and \
       (k < kMax + 0):
        n = (Ni * i) + (Nj * j) + (Nk * k) + N0
        m = iHy[n]
        if m == 0:
            Hy[n] += \
                  - RZc[k] * (Ex[n + Nk] - Ex[n]) \
                  + RXc[i] * (Ez[n + Ni] - Ez[n])
        else:
            fi, dfi = sol_cuda.planewave.f(Xc[i], Yn[j], Zc[k], t, fPlanewave[4], fPlanewave)
            Hy[n] = K1Hy[n] * Hy[n] \
                  - K2Hy[n] * (RZc[k] * (Ex[n + Nk] - Ex[n])
                             - RXc[i] * (Ez[n + Ni] - Ez[n])) \
                  - (K1Hy[n] - K2Hy[n]) * dfi \
                  - (1 - K1Hy[n]) * fi

# (4)
@cuda.jit(cache=True)
def _Hy_p_no_vector(
    t, fPlanewave, Yn, Zc, Xc,
    Hy, Ez, Ex, iHy, C1H, C2H, RZc, RXc,
    iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0):

    k, j, i = cuda.grid(3)
    i += iMin
    j += jMin
    k += kMin

    if (i < iMax + 0) and \
       (j < jMax + 1) and \
       (k < kMax + 0):
        n = (Ni * i) + (Nj * j) + (Nk * k) + N0
        m = iHy[n]
        if m == 0:
            Hy[n] += \
                  - RZc[k] * (Ex[n + Nk] - Ex[n]) \
                  + RXc[i] * (Ez[n + Ni] - Ez[n])
        else:
            fi, dfi = sol_cuda.planewave.f(Xc[i], Yn[j], Zc[k], t, fPlanewave[4], fPlanewave)
            Hy[n] = C1H[m] * Hy[n] \
                  - C2H[m] * (RZc[k] * (Ex[n + Nk] - Ex[n])
                            - RXc[i] * (Ez[n + Ni] - Ez[n])) \
                  - (C1H[m] - C2H[m]) * dfi \
                  - (1 - C1H[m]) * fi
