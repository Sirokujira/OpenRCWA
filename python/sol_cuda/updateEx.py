# -*- coding: utf-8 -*-
"""
updateEx.py (CUDA)
"""

import math
from numba import cuda
import sol_cuda.planewave

# Ex更新
def calEx(block3d,
    Parm, t, fPlanewave, Xc, Yn, Zn, Ex, Hy, Hz, iEx,
    C1E, C2E, K1Ex, K2Ex, RYn, RZn,
    iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0):

    grid = (math.ceil((kMax - kMin + 1) / block3d[0]),
            math.ceil((jMax - jMin + 1) / block3d[1]),
            math.ceil((iMax - iMin + 0) / block3d[2]))

    if   Parm['source'] == 0:
        # 給電点
        if Parm['vector']:
            # (1)
            _Ex_f_vector[grid, block3d](
                Ex, Hy, Hz, K1Ex, K2Ex, RYn, RZn,
                iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)
        else:
            # (2)
            _Ex_f_no_vector[grid, block3d](
                Ex, Hy, Hz, iEx, C1E, C2E, RYn, RZn,
                iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)
    elif Parm['source'] == 1:
        # 平面波入射
        if Parm['vector']:
            # (3)
            _Ex_p_vector[grid, block3d](
                t, fPlanewave, Xc, Yn, Zn,
                Ex, Hy, Hz, iEx, K1Ex, K2Ex, RYn, RZn,
                iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)
        else:
            # (4)
            _Ex_p_no_vector[grid, block3d](
                t, fPlanewave, Xc, Yn, Zn,
                Ex, Hy, Hz, iEx, C1E, C2E, RYn, RZn,
                iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)

# 給電点
# (1)
@cuda.jit(cache=True)
def _Ex_f_vector(
    Ex, Hy, Hz, K1Ex, K2Ex, RYn, RZn,
    iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0):

    k, j, i = cuda.grid(3)
    i += iMin
    j += jMin
    k += kMin

    if (i < iMax + 0) and \
       (j < jMax + 1) and \
       (k < kMax + 1):
        n = (Ni * i) + (Nj * j) + (Nk * k) + N0
        Ex[n] = K1Ex[n] * Ex[n] \
              + K2Ex[n] * (RYn[j] * (Hz[n] - Hz[n - Nj])
                         - RZn[k] * (Hy[n] - Hy[n - Nk]))

# (2)
@cuda.jit(cache=True)
def _Ex_f_no_vector(
    Ex, Hy, Hz, iEx, C1E, C2E, RYn, RZn,
    iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0):

    k, j, i = cuda.grid(3)
    i += iMin
    j += jMin
    k += kMin

    if (i < iMax + 0) and \
       (j < jMax + 1) and \
       (k < kMax + 1):
        n = (Ni * i) + (Nj * j) + (Nk * k) + N0
        m = iEx[n]
        Ex[n] = C1E[m] * Ex[n] \
              + C2E[m] * (RYn[j] * (Hz[n] - Hz[n - Nj])
                        - RZn[k] * (Hy[n] - Hy[n - Nk]))

# 平面波入射
# (3)
@cuda.jit(cache=True)
def _Ex_p_vector(
    t, fPlanewave, Xc, Yn, Zn,
    Ex, Hy, Hz, iEx, K1Ex, K2Ex, RYn, RZn,
    iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0):

    k, j, i = cuda.grid(3)
    i += iMin
    j += jMin
    k += kMin

    if (i < iMax + 0) and \
       (j < jMax + 1) and \
       (k < kMax + 1):
        n = (Ni * i) + (Nj * j) + (Nk * k) + N0
        m = iEx[n]
        if m == 0:
            Ex[n] += \
                  + RYn[j] * (Hz[n] - Hz[n - Nj]) \
                  - RZn[k] * (Hy[n] - Hy[n - Nk])
        else:
            fi, dfi = sol_cuda.planewave.f(Xc[i], Yn[j], Zn[k], t, fPlanewave[0], fPlanewave)
            Ex[n] = K1Ex[n] * Ex[n] \
                  + K2Ex[n] * (RYn[j] * (Hz[n] - Hz[n - Nj])
                             - RZn[k] * (Hy[n] - Hy[n - Nk])) \
                  - (K1Ex[n] - K2Ex[n]) * dfi \
                  - (1 - K1Ex[n]) * fi

# (4)
@cuda.jit(cache=True)
def _Ex_p_no_vector(
    t, fPlanewave, Xc, Yn, Zn,
    Ex, Hy, Hz, iEx, C1E, C2E, RYn, RZn,
    iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0):

    k, j, i = cuda.grid(3)
    i += iMin
    j += jMin
    k += kMin

    if (i < iMax + 0) and \
       (j < jMax + 1) and \
       (k < kMax + 1):
        n = (Ni * i) + (Nj * j) + (Nk * k) + N0
        m = iEx[n]
        if m == 0:
            Ex[n] += \
                  + RYn[j] * (Hz[n] - Hz[n - Nj]) \
                  - RZn[k] * (Hy[n] - Hy[n - Nk])
        else:
            fi, dfi = sol_cuda.planewave.f(Xc[i], Yn[j], Zn[k], t, fPlanewave[0], fPlanewave)
            Ex[n] = C1E[m] * Ex[n] \
                  + C2E[m] * (RYn[j] * (Hz[n] - Hz[n - Nj])
                            - RZn[k] * (Hy[n] - Hy[n - Nk])) \
                  - (C1E[m] - C2E[m]) * dfi \
                  - (1 - C1E[m]) * fi
