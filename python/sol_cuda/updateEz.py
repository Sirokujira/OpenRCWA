# -*- coding: utf-8 -*-
"""
updateEz.py (CUDA)
"""

import math
from numba import cuda
import sol_cuda.planewave

# Ez更新
def calEz(block3d,
    Parm, t, fPlanewave, Zc, Xn, Yn, Ez, Hx, Hy, iEz,
    C1E, C2E, K1Ez, K2Ez, RXn, RYn,
    iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0):

    grid = (math.ceil((kMax - kMin + 0) / block3d[0]),
            math.ceil((jMax - jMin + 1) / block3d[1]),
            math.ceil((iMax - iMin + 1) / block3d[2]))

    if   Parm['source'] == 0:
        # 給電点
        if Parm['vector']:
            # (1)
            _Ez_f_vector[grid, block3d](
                Ez, Hx, Hy, K1Ez, K2Ez, RXn, RYn,
                iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)
        else:
            # (2)
            _Ez_f_no_vector[grid, block3d](
                Ez, Hx, Hy, iEz, C1E, C2E, RXn, RYn,
                iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)
    elif Parm['source'] == 1:
        # 平面波入射
        if Parm['vector']:
            # (3)
            _Ez_p_vector[grid, block3d](
                t, fPlanewave, Zc, Xn, Yn,
                Ez, Hx, Hy, iEz, K1Ez, K2Ez, RXn, RYn,
                iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)
        else:
            # (4)
            _Ez_p_no_vector[grid, block3d](
                t, fPlanewave, Zc, Xn, Yn,
                Ez, Hx, Hy, iEz, C1E, C2E, RXn, RYn,
                iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)

# 給電点
# (1)
@cuda.jit(cache=True)
def _Ez_f_vector(
    Ez, Hx, Hy, K1Ez, K2Ez, RXn, RYn,
    iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0):

    k, j, i = cuda.grid(3)
    i += iMin
    j += jMin
    k += kMin

    if (i < iMax + 1) and \
       (j < jMax + 1) and \
       (k < kMax + 0):
        n = (Ni * i) + (Nj * j) + (Nk * k) + N0
        Ez[n] = K1Ez[n] * Ez[n] \
              + K2Ez[n] * (RXn[i] * (Hy[n] - Hy[n - Ni])
                         - RYn[j] * (Hx[n] - Hx[n - Nj]))

# (2)
@cuda.jit(cache=True)
def _Ez_f_no_vector(
    Ez, Hx, Hy, iEz, C1E, C2E, RXn, RYn,
    iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0):

    k, j, i = cuda.grid(3)
    i += iMin
    j += jMin
    k += kMin

    if (i < iMax + 1) and \
       (j < jMax + 1) and \
       (k < kMax + 0):
        n = (Ni * i) + (Nj * j) + (Nk * k) + N0
        m = iEz[n]
        Ez[n] = C1E[m] * Ez[n] \
              + C2E[m] * (RXn[i] * (Hy[n] - Hy[n - Ni])
                        - RYn[j] * (Hx[n] - Hx[n - Nj]))

# 平面波入射
# (3)
@cuda.jit(cache=True)
def _Ez_p_vector(
    t, fPlanewave, Zc, Xn, Yn,
    Ez, Hx, Hy, iEz, K1Ez, K2Ez, RXn, RYn,
    iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0):

    k, j, i = cuda.grid(3)
    i += iMin
    j += jMin
    k += kMin

    if (i < iMax + 1) and \
       (j < jMax + 1) and \
       (k < kMax + 0):
        n = (Ni * i) + (Nj * j) + (Nk * k) + N0
        m = iEz[n]
        if m == 0:
            Ez[n] += \
                  + RXn[i] * (Hy[n] - Hy[n - Ni]) \
                  - RYn[j] * (Hx[n] - Hx[n - Nj])
        else:
            fi, dfi = sol_cuda.planewave.f(Xn[i], Yn[j], Zc[k], t, fPlanewave[2], fPlanewave)
            Ez[n] = K1Ez[n] * Ez[n] \
                  + K2Ez[n] * (RXn[i] * (Hy[n] - Hy[n - Ni])
                             - RYn[j] * (Hx[n] - Hx[n - Nj])) \
                  - (K1Ez[n] - K2Ez[n]) * dfi \
                  - (1 - K1Ez[n]) * fi

# (4)
@cuda.jit(cache=True)
def _Ez_p_no_vector(
    t, fPlanewave, Zc, Xn, Yn,
    Ez, Hx, Hy, iEz, C1E, C2E, RXn, RYn,
    iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0):

    k, j, i = cuda.grid(3)
    i += iMin
    j += jMin
    k += kMin

    if (i < iMax + 1) and \
       (j < jMax + 1) and \
       (k < kMax + 0):
        n = (Ni * i) + (Nj * j) + (Nk * k) + N0
        m = iEz[n]
        if m == 0:
            Ez[n] += \
                  + RXn[i] * (Hy[n] - Hy[n - Ni]) \
                  - RYn[j] * (Hx[n] - Hx[n - Nj])
        else:
            fi, dfi = sol_cuda.planewave.f(Xn[i], Yn[j], Zc[k], t, fPlanewave[2], fPlanewave)
            Ez[n] = C1E[m] * Ez[n] \
                  + C2E[m] * (RXn[i] * (Hy[n] - Hy[n - Ni])
                            - RYn[j] * (Hx[n] - Hx[n - Nj])) \
                  - (C1E[m] - C2E[m]) * dfi \
                  - (1 - C1E[m]) * fi
