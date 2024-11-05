# -*- coding: utf-8 -*-
"""
updateHx.py
"""

from numba import jit, prange
import sol.planewave

# Hx更新
def calHx(
    Parm, t, fPlanewave, Xn, Yc, Zc, Hx, Ey, Ez, iHx,
    C1H, C2H, K1Hx, K2Hx, RYc, RZc,
    iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0):

    if   Parm['source'] == 0:
        # 給電点
        if Parm['vector']:
            _Hx_f_vector(
                Hx, Ey, Ez, K1Hx, K2Hx, RYc, RZc,
                iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)
        else:
            _Hx_f_no_vector(
                Hx, Ey, Ez, iHx, C1H, C2H, RYc, RZc,
                iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)
    elif Parm['source'] == 1:
        # 平面波入射
        if Parm['vector']:
            _Hx_p_vector(
                t, fPlanewave[3], fPlanewave, Xn, Yc, Zc,
                Hx, Ey, Ez, iHx, K1Hx, K2Hx, RYc, RZc,
                iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)
        else:
            _Hx_p_no_vector(
                t, fPlanewave[3], fPlanewave, Xn, Yc, Zc,
                Hx, Ey, Ez, iHx, C1H, C2H, RYc, RZc,
                iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)

# 給電点
@jit(cache=True, nogil=True, parallel=True, nopython=True)
def _Hx_f_vector(
    Hx, Ey, Ez, K1Hx, K2Hx, RYc, RZc,
    iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0):

    for i in prange(iMin, iMax + 1):
        for j in range(jMin, jMax + 0):
            for k in range(kMin, kMax + 0):
                n = (Ni * i) + (Nj * j) + (Nk * k) + N0
                Hx[n] = K1Hx[n] * Hx[n] \
                      - K1Hx[n] * (RYc[j] * (Ez[n + Nj] - Ez[n])
                                 - RZc[k] * (Ey[n + Nk] - Ey[n]))

@jit(cache=True, nogil=True, parallel=True, nopython=True)
def _Hx_f_no_vector(
    Hx, Ey, Ez, iHx, C1H, C2H, RYc, RZc,
    iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0):

    for i in prange(iMin, iMax + 1):
        for j in range(jMin, jMax + 0):
            for k in range(kMin, kMax + 0):
                n = (Ni * i) + (Nj * j) + (Nk * k) + N0
                m = iHx[n]
                Hx[n] = C1H[m] * Hx[n] \
                      - C2H[m] * (RYc[j] * (Ez[n + Nj] - Ez[n])
                                - RZc[k] * (Ey[n + Nk] - Ey[n]))

# 平面波入射
@jit(cache=True, nogil=True, parallel=True, nopython=True)
def _Hx_p_vector(
    t, hx, fPlanewave, Xn, Yc, Zc,
    Hx, Ey, Ez, iHx, K1Hx, K2Hx, RYc, RZc,
    iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0):

    for i in prange(iMin, iMax + 1):
        for j in range(jMin, jMax + 0):
            for k in range(kMin, kMax + 0):
                n = (Ni * i) + (Nj * j) + (Nk * k) + N0
                m = iHx[n]
                if m == 0:
                    Hx[n] += \
                          - RYc[j] * (Ez[n + Nj] - Ez[n]) \
                          + RZc[k] * (Ey[n + Nk] - Ey[n])
                else:
                    fi, dfi = sol.planewave.f(Xn[i], Yc[j], Zc[k], t, hx, fPlanewave)
                    Hx[n] = K1Hx[n] * Hx[n] \
                          - K2Hx[n] * (RYc[j] * (Ez[n + Nj] - Ez[n])
                                     - RZc[k] * (Ey[n + Nk] - Ey[n])) \
                          - (K1Hx[n] - K2Hx[n]) * dfi \
                          - (1 - K1Hx[n]) * fi

@jit(cache=True, nogil=True, parallel=True, nopython=True)
def _Hx_p_no_vector(
    t, hx, fPlanewave, Xn, Yc, Zc,
    Hx, Ey, Ez, iHx, C1H, C2H, RYc, RZc,
    iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0):

    for i in prange(iMin, iMax + 1):
        for j in range(jMin, jMax + 0):
            for k in range(kMin, kMax + 0):
                n = (Ni * i) + (Nj * j) + (Nk * k) + N0
                m = iHx[n]
                if m == 0:
                    Hx[n] += \
                          - RYc[j] * (Ez[n + Nj] - Ez[n]) \
                          + RZc[k] * (Ey[n + Nk] - Ey[n])
                else:
                    fi, dfi = sol.planewave.f(Xn[i], Yc[j], Zc[k], t, hx, fPlanewave)
                    Hx[n] = C1H[m] * Hx[n] \
                          - C2H[m] * (RYc[j] * (Ez[n + Nj] - Ez[n])
                                    - RZc[k] * (Ey[n + Nk] - Ey[n])) \
                          - (C1H[m] - C2H[m]) * dfi \
                          - (1 - C1H[m]) * fi
