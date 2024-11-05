# -*- coding: utf-8 -*-
"""
updateHy.py
"""

from numba import jit, prange
import sol.planewave

# Hy更新
def calHy(
    Parm, t, fPlanewave, Yn, Zc, Xc, Hy, Ez, Ex, iHy,
    C1H, C2H, K1Hy, K2Hy, RZc, RXc,
    iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0):

    if   Parm['source'] == 0:
        # 給電点
        if Parm['vector']:
            _Hy_f_vector(
                Hy, Ez, Ex, K1Hy, K2Hy, RZc, RXc,
                iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)
        else:
            _Hy_f_no_vector(
                Hy, Ez, Ex, iHy, C1H, C2H, RZc, RXc,
                iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)
    elif Parm['source'] == 1:
        # 平面波入射
        if Parm['vector']:
            _Hy_p_vector(
                t, fPlanewave[4], fPlanewave, Yn, Zc, Xc,
                Hy, Ez, Ex, iHy, K1Hy, K2Hy, RZc, RXc,
                iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)
        else:
            _Hy_p_no_vector(
                t, fPlanewave[4], fPlanewave, Yn, Zc, Xc,
                Hy, Ez, Ex, iHy, C1H, C2H, RZc, RXc,
                iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)

# 給電点
@jit(cache=True, nogil=True, parallel=True, nopython=True)
def _Hy_f_vector(
    Hy, Ez, Ex, K1Hy, K2Hy, RZc, RXc,
    iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0):

    for i in prange(iMin, iMax + 0):
        for j in range(jMin, jMax + 1):
            for k in range(kMin, kMax + 0):
                n = (Ni * i) + (Nj * j) + (Nk * k) + N0
                Hy[n] = K1Hy[n] * Hy[n] \
                      - K2Hy[n] * (RZc[k] * (Ex[n + Nk] - Ex[n])
                                 - RXc[i] * (Ez[n + Ni] - Ez[n]))

@jit(cache=True, nogil=True, parallel=True, nopython=True)
def _Hy_f_no_vector(
    Hy, Ez, Ex, iHy, C1H, C2H, RZc, RXc,
    iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0):

    for i in prange(iMin, iMax + 0):
        for j in range(jMin, jMax + 1):
            for k in range(kMin, kMax + 0):
                n = (Ni * i) + (Nj * j) + (Nk * k) + N0
                m = iHy[n]
                Hy[n] = C1H[m] * Hy[n] \
                      - C2H[m] * (RZc[k] * (Ex[n + Nk] - Ex[n])
                                - RXc[i] * (Ez[n + Ni] - Ez[n]))

# 平面波入射
@jit(cache=True, nogil=True, parallel=True, nopython=True)
def _Hy_p_vector(
    t, hy, fPlanewave, Yn, Zc, Xc,
    Hy, Ez, Ex, iHy, K1Hy, K2Hy, RZc, RXc,
    iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0):

    for i in prange(iMin, iMax + 0):
        for j in range(jMin, jMax + 1):
            for k in range(kMin, kMax + 0):
                n = (Ni * i) + (Nj * j) + (Nk * k) + N0
                m = iHy[n]
                if m == 0:
                    Hy[n] += \
                          - RZc[k] * (Ex[n + Nk] - Ex[n]) \
                          + RXc[i] * (Ez[n + Ni] - Ez[n])
                else:
                    fi, dfi = sol.planewave.f(Xc[i], Yn[j], Zc[k], t, hy, fPlanewave)
                    Hy[n] = K1Hy[n] * Hy[n] \
                          - K2Hy[n] * (RZc[k] * (Ex[n + Nk] - Ex[n])
                                     - RXc[i] * (Ez[n + Ni] - Ez[n])) \
                          - (K1Hy[n] - K2Hy[n]) * dfi \
                          - (1 - K1Hy[n]) * fi

@jit(cache=True, nogil=True, parallel=True, nopython=True)
def _Hy_p_no_vector(
    t, hy, fPlanewave, Yn, Zc, Xc,
    Hy, Ez, Ex, iHy, C1H, C2H, RZc, RXc,
    iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0):

    for i in prange(iMin, iMax + 0):
        for j in range(jMin, jMax + 1):
            for k in range(kMin, kMax + 0):
                n = (Ni * i) + (Nj * j) + (Nk * k) + N0
                m = iHy[n]
                if m == 0:
                    Hy[n] += \
                          - RZc[k] * (Ex[n + Nk] - Ex[n]) \
                          + RXc[i] * (Ez[n + Ni] - Ez[n])
                else:
                    fi, dfi = sol.planewave.f(Xc[i], Yn[j], Zc[k], t, hy, fPlanewave)
                    Hy[n] = C1H[m] * Hy[n] \
                          - C2H[m] * (RZc[k] * (Ex[n + Nk] - Ex[n])
                                    - RXc[i] * (Ez[n + Ni] - Ez[n])) \
                          - (C1H[m] - C2H[m]) * dfi \
                          - (1 - C1H[m]) * fi
