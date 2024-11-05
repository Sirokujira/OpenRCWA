# -*- coding: utf-8 -*-
"""
updateHx.py (CUDA)
"""

import math
from numba import cuda
import sol_cuda.planewave

# Hx更新
def calHx(block3d,
    Parm, t, fPlanewave, Xn, Yc, Zc, Hx, Ey, Ez, iHx,
    C1H, C2H, K1Hx, K2Hx, RYc, RZc,
    iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0):

    grid = (math.ceil((kMax - kMin + 0) / block3d[0]),
            math.ceil((jMax - jMin + 0) / block3d[1]),
            math.ceil((iMax - iMin + 1) / block3d[2]))

    if   Parm['source'] == 0:
        # 給電点
        if Parm['vector']:
            # (1)
            _Hx_f_vector[grid, block3d](
                Hx, Ey, Ez, K1Hx, K2Hx, RYc, RZc,
                iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)
        else:
            # (2)
            _Hx_f_no_vector[grid, block3d](
                Hx, Ey, Ez, iHx, C1H, C2H, RYc, RZc,
                iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)
    elif Parm['source'] == 1:
        # 平面波入射
        if Parm['vector']:
            # (3)
            _Hx_p_vector[grid, block3d](
                t, fPlanewave, Xn, Yc, Zc,
                Hx, Ey, Ez, iHx, K1Hx, K2Hx, RYc, RZc,
                iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)
        else:
            # (4)
            _Hx_p_no_vector[grid, block3d](
                t, fPlanewave, Xn, Yc, Zc,
                Hx, Ey, Ez, iHx, C1H, C2H, RYc, RZc,
                iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)

# 給電点
# (1)
@cuda.jit(cache=True)
def _Hx_f_vector(
    Hx, Ey, Ez, K1Hx, K2Hx, RYc, RZc,
    iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0):

    k, j, i = cuda.grid(3)
    i += iMin
    j += jMin
    k += kMin

    if (i < iMax + 1) and \
       (j < jMax + 0) and \
       (k < kMax + 0):
        n = (Ni * i) + (Nj * j) + (Nk * k) + N0
        Hx[n] = K1Hx[n] * Hx[n] \
              - K2Hx[n] * (RYc[j] * (Ez[n + Nj] - Ez[n])
                         - RZc[k] * (Ey[n + Nk] - Ey[n]))

# (2)
@cuda.jit(cache=True)
def _Hx_f_no_vector(
    Hx, Ey, Ez, iHx, C1H, C2H, RYc, RZc,
    iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0):

    k, j, i = cuda.grid(3)
    i += iMin
    j += jMin
    k += kMin

    if (i < iMax + 1) and \
       (j < jMax + 0) and \
       (k < kMax + 0):
        n = (Ni * i) + (Nj * j) + (Nk * k) + N0
        m = iHx[n]
        Hx[n] = C1H[m] * Hx[n] \
              - C2H[m] * (RYc[j] * (Ez[n + Nj] - Ez[n])
                        - RZc[k] * (Ey[n + Nk] - Ey[n]))

# 平面波入射
# (3)
@cuda.jit(cache=True)
def _Hx_p_vector(
    t, fPlanewave, Xn, Yc, Zc,
    Hx, Ey, Ez, iHx, K1Hx, K2Hx, RYc, RZc,
    iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0):

    k, j, i = cuda.grid(3)
    i += iMin
    j += jMin
    k += kMin

    if (i < iMax + 1) and \
       (j < jMax + 0) and \
       (k < kMax + 0):
        n = (Ni * i) + (Nj * j) + (Nk * k) + N0
        m = iHx[n]
        if m == 0:
            Hx[n] += \
                  - RYc[j] * (Ez[n + Nj] - Ez[n]) \
                  + RZc[k] * (Ey[n + Nk] - Ey[n])
        else:
            fi, dfi = sol_cuda.planewave.f(Xn[i], Yc[j], Zc[k], t, fPlanewave[3], fPlanewave)
            Hx[n] = K1Hx[n] * Hx[n] \
                  - K2Hx[n] * (RYc[j] * (Ez[n + Nj] - Ez[n])
                             - RZc[k] * (Ey[n + Nk] - Ey[n])) \
                  - (K1Hx[n] - K2Hx[n]) * dfi \
                  - (1 - K1Hx[n]) * fi

# (4)
@cuda.jit(cache=True)
def _Hx_p_no_vector(
    t, fPlanewave, Xn, Yc, Zc,
    Hx, Ey, Ez, iHx, C1H, C2H, RYc, RZc,
    iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0):

    k, j, i = cuda.grid(3)
    i += iMin
    j += jMin
    k += kMin

    if (i < iMax + 1) and \
       (j < jMax + 0) and \
       (k < kMax + 0):
        n = (Ni * i) + (Nj * j) + (Nk * k) + N0
        m = iHx[n]
        if m == 0:
            Hx[n] += \
                  - RYc[j] * (Ez[n + Nj] - Ez[n]) \
                  + RZc[k] * (Ey[n + Nk] - Ey[n])
        else:
            fi, dfi = sol_cuda.planewave.f(Xn[i], Yc[j], Zc[k], t, fPlanewave[3], fPlanewave)
            Hx[n] = C1H[m] * Hx[n] \
                  - C2H[m] * (RYc[j] * (Ez[n + Nj] - Ez[n])
                            - RZc[k] * (Ey[n + Nk] - Ey[n])) \
                  - (C1H[m] - C2H[m]) * dfi \
                  - (1 - C1H[m]) * fi
