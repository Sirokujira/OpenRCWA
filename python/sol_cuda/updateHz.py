# -*- coding: utf-8 -*-
"""
updateHz.py (CUDA)
"""

import math
from numba import cuda
import sol_cuda.planewave

# Hz更新
def calHz(block3d,
    Parm, t, fPlanewave, Zn, Xc, Yc, Hz, Ex, Ey, iHz,
    C1H, C2H, K1Hz, K2Hz, RXc, RYc,
    iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0):

    grid = (math.ceil((kMax - kMin + 1) / block3d[0]),
            math.ceil((jMax - jMin + 0) / block3d[1]),
            math.ceil((iMax - iMin + 0) / block3d[2]))

    if   Parm['source'] == 0:
        # 給電点
        if Parm['vector']:
            # (1)
            _Hz_f_vector[grid, block3d](
                Hz, Ex, Ey, K1Hz, K2Hz, RXc, RYc,
                iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)
        else:
            # (2)
            _Hz_f_no_vector[grid, block3d](
                Hz, Ex, Ey, iHz, C1H, C2H, RXc, RYc,
                iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)
    elif Parm['source'] == 1:
        # 平面波入射
        if Parm['vector']:
            # (3)
            _Hz_p_vector[grid, block3d](
                t, fPlanewave, Zn, Xc, Yc,
                Hz, Ex, Ey, iHz, K1Hz, K2Hz, RXc, RYc,
                iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)
        else:
            # (4)
            _Hz_p_no_vector[grid, block3d](
                t, fPlanewave, Zn, Xc, Yc,
                Hz, Ex, Ey, iHz, C1H, C2H, RXc, RYc,
                iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)

# 給電点
# (1)
@cuda.jit(cache=True)
def _Hz_f_vector(
    Hz, Ex, Ey, K1Hz, K2Hz, RXc, RYc,
    iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0):

    k, j, i = cuda.grid(3)
    i += iMin
    j += jMin
    k += kMin

    if (i < iMax + 0) and \
       (j < jMax + 0) and \
       (k < kMax + 1):
        n = (Ni * i) + (Nj * j) + (Nk * k) + N0
        Hz[n] = K1Hz[n] * Hz[n] \
              - K2Hz[n] * (RXc[i] * (Ey[n + Ni] - Ey[n])
                         - RYc[j] * (Ex[n + Nj] - Ex[n]))

# (2)
@cuda.jit(cache=True)
def _Hz_f_no_vector(
    Hz, Ex, Ey, iHz, C1H, C2H, RXc, RYc,
    iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0):

    k, j, i = cuda.grid(3)
    i += iMin
    j += jMin
    k += kMin

    if (i < iMax + 0) and \
       (j < jMax + 0) and \
       (k < kMax + 1):
        n = (Ni * i) + (Nj * j) + (Nk * k) + N0
        m = iHz[n]
        Hz[n] = C1H[m] * Hz[n] \
              - C2H[m] * (RXc[i] * (Ey[n + Ni] - Ey[n])
                        - RYc[j] * (Ex[n + Nj] - Ex[n]))

# 平面波入射
# (3)
@cuda.jit(cache=True)
def _Hz_p_vector(
    t, fPlanewave, Zn, Xc, Yc,
    Hz, Ex, Ey, iHz, K1Hz, K2Hz, RXc, RYc,
    iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0):

    k, j, i = cuda.grid(3)
    i += iMin
    j += jMin
    k += kMin

    if (i < iMax + 0) and \
       (j < jMax + 0) and \
       (k < kMax + 1):
        n = (Ni * i) + (Nj * j) + (Nk * k) + N0
        m = iHz[n]
        if m == 0:
            Hz[n] += \
                  - RXc[i] * (Ey[n + Ni] - Ey[n]) \
                  + RYc[j] * (Ex[n + Nj] - Ex[n])
        else:
            fi, dfi = sol_cuda.planewave.f(Xc[i], Yc[j], Zn[k], t, fPlanewave[5], fPlanewave)
            Hz[n] = K1Hz[n] * Hz[n] \
                  - K2Hz[n] * (RXc[i] * (Ey[n + Ni] - Ey[n])
                             - RYc[j] * (Ex[n + Nj] - Ex[n])) \
                  - (K1Hz[n] - K2Hz[n]) * dfi \
                  - (1 - K1Hz[n]) * fi

# (4)
@cuda.jit(cache=True)
def _Hz_p_no_vector(
    t, fPlanewave, Zn, Xc, Yc,
    Hz, Ex, Ey, iHz, C1H, C2H, RXc, RYc,
    iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0):

    k, j, i = cuda.grid(3)
    i += iMin
    j += jMin
    k += kMin

    if (i < iMax + 0) and \
       (j < jMax + 0) and \
       (k < kMax + 1):
        n = (Ni * i) + (Nj * j) + (Nk * k) + N0
        m = iHz[n]
        if m == 0:
            Hz[n] += \
                  - RXc[i] * (Ey[n + Ni] - Ey[n]) \
                  + RYc[j] * (Ex[n + Nj] - Ex[n])
        else:
            fi, dfi = sol_cuda.planewave.f(Xc[i], Yc[j], Zn[k], t, fPlanewave[5], fPlanewave)
            Hz[n] = C1H[m] * Hz[n] \
                  - C2H[m] * (RXc[i] * (Ey[n + Ni] - Ey[n])
                            - RYc[j] * (Ex[n + Nj] - Ex[n])) \
                  - (C1H[m] - C2H[m]) * dfi \
                  - (1 - C1H[m]) * fi
