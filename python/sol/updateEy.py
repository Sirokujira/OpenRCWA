# -*- coding: utf-8 -*-
"""
updateEy.py
"""

from numba import jit, prange
import sol.planewave

# Ey更新
def calEy(
    Parm, t, fPlanewave, Yc, Zn, Xn, Ey, Hz, Hx, iEy,
    C1E, C2E, K1Ey, K2Ey, RZn, RXn,
    iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0):

    if   Parm['source'] == 0:
        # 給電点
        if Parm['vector']:
            _Ey_f_vector(
                Ey, Hz, Hx, K1Ey, K2Ey, RZn, RXn,
                iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)
        else:
            _Ey_f_no_vector(
                Ey, Hz, Hx, iEy, C1E, C2E, RZn, RXn,
                iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)
    elif Parm['source'] == 1:
        # 平面波入射
        if Parm['vector']:
            _Ey_p_vector(
                t, fPlanewave[1], fPlanewave, Yc, Zn, Xn,
                Ey, Hz, Hx, iEy, K1Ey, K2Ey, RZn, RXn,
                iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)
        else:
            _Ey_p_no_vector(
                t, fPlanewave[1], fPlanewave, Yc, Zn, Xn,
                Ey, Hz, Hx, iEy, C1E, C2E, RZn, RXn,
                iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)

# 給電点
@jit(cache=True, nogil=True, parallel=True, nopython=True)
def _Ey_f_vector(
    Ey, Hz, Hx, K1Ey, K2Ey, RZn, RXn,
    iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0):

    for i in prange(iMin, iMax + 1):
        for j in range(jMin, jMax + 0):
            for k in range(kMin, kMax + 1):
                n = (Ni * i) + (Nj * j) + (Nk * k) + N0
                Ey[n] = K1Ey[n] * Ey[n] \
                      + K2Ey[n] * (RZn[k] * (Hx[n] - Hx[n - Nk])
                                 - RXn[i] * (Hz[n] - Hz[n - Ni]))

@jit(cache=True, nogil=True, parallel=True, nopython=True)
def _Ey_f_no_vector(
    Ey, Hz, Hx, iEy, C1E, C2E, RZn, RXn,
    iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0):

    for i in prange(iMin, iMax + 1):
        for j in range(jMin, jMax + 0):
            for k in range(kMin, kMax + 1):
                n = (Ni * i) + (Nj * j) + (Nk * k) + N0
                m = iEy[n]
                Ey[n] = C1E[m] * Ey[n] \
                      + C2E[m] * (RZn[k] * (Hx[n] - Hx[n - Nk])
                                - RXn[i] * (Hz[n] - Hz[n - Ni]))

# 平面波入射
@jit(cache=True, nogil=True, parallel=True, nopython=True)
def _Ey_p_vector(
    t, ey, fPlanewave, Yc, Zn, Xn,
    Ey, Hz, Hx, iEy, K1Ey, K2Ey, RZn, RXn,
    iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0):

    for i in prange(iMin, iMax + 1):
        for j in range(jMin, jMax + 0):
            for k in range(kMin, kMax + 1):
                n = (Ni * i) + (Nj * j) + (Nk * k) + N0
                m = iEy[n]
                if m == 0:
                    Ey[n] += \
                          + RZn[k] * (Hx[n] - Hx[n - Nk]) \
                          - RXn[i] * (Hz[n] - Hz[n - Ni])
                else:
                    fi, dfi = sol.planewave.f(Xn[i], Yc[j], Zn[k], t, ey, fPlanewave)
                    Ey[n] = K1Ey[n] * Ey[n] \
                          + K2Ey[n] * (RZn[k] * (Hx[n] - Hx[n - Nk])
                                     - RXn[i] * (Hz[n] - Hz[n - Ni])) \
                          - (K1Ey[n] - K2Ey[n]) * dfi \
                          - (1 - K1Ey[n]) * fi

@jit(cache=True, nogil=True, parallel=True, nopython=True)
def _Ey_p_no_vector(
    t, ey, fPlanewave, Yc, Zn, Xn,
    Ey, Hz, Hx, iEy, C1E, C2E, RZn, RXn,
    iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0):

    for i in prange(iMin, iMax + 1):
        for j in range(jMin, jMax + 0):
            for k in range(kMin, kMax + 1):
                n = (Ni * i) + (Nj * j) + (Nk * k) + N0
                m = iEy[n]
                if m == 0:
                    Ey[n] += \
                          + RZn[k] * (Hx[n] - Hx[n - Nk]) \
                          - RXn[i] * (Hz[n] - Hz[n - Ni])
                else:
                    fi, dfi = sol.planewave.f(Xn[i], Yc[j], Zn[k], t, ey, fPlanewave)
                    Ey[n] = C1E[m] * Ey[n] \
                          + C2E[m] * (RZn[k] * (Hx[n] - Hx[n - Nk])
                                    - RXn[i] * (Hz[n] - Hz[n - Ni])) \
                          - (C1E[m] - C2E[m]) * dfi \
                          - (1 - C1E[m]) * fi
