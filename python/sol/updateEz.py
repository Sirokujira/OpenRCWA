# -*- coding: utf-8 -*-
"""
updateEz.py
"""

from numba import jit, prange
import sol.planewave

# Ez更新
def calEz(
    Parm, t, fPlanewave, Zc, Xn, Yn, Ez, Hx, Hy, iEz,
    C1E, C2E, K1Ez, K2Ez, RXn, RYn,
    iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0):

    if   Parm['source'] == 0:
        # 給電点
        if Parm['vector']:
            _Ez_f_vector(
                Ez, Hx, Hy, K1Ez, K2Ez, RXn, RYn,
                iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)
        else:
            _Ez_f_no_vector(
                Ez, Hx, Hy, iEz, C1E, C2E, RXn, RYn,
                iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)
    elif Parm['source'] == 1:
        # 平面波入射
        if Parm['vector']:
            _Ez_p_vector(
                t, fPlanewave[2], fPlanewave, Zc, Xn, Yn,
                Ez, Hx, Hy, iEz, K1Ez, K2Ez, RXn, RYn,
                iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)
        else:
            _Ez_p_no_vector(
                t, fPlanewave[2], fPlanewave, Zc, Xn, Yn,
                Ez, Hx, Hy, iEz, C1E, C2E, RXn, RYn,
                iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)

# 給電点
@jit(cache=True, nogil=True, parallel=True, nopython=True)
def _Ez_f_vector(
    Ez, Hx, Hy, K1Ez, K2Ez, RXn, RYn,
    iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0):

    for i in prange(iMin, iMax + 1):
        for j in range(jMin, jMax + 1):
            for k in range(kMin, kMax + 0):
                n = (Ni * i) + (Nj * j) + (Nk * k) + N0
                Ez[n] = K1Ez[n] * Ez[n] \
                      + K2Ez[n] * (RXn[i] * (Hy[n] - Hy[n - Ni])
                                 - RYn[j] * (Hx[n] - Hx[n - Nj]))

@jit(cache=True, nogil=True, parallel=True, nopython=True)
def _Ez_f_no_vector(
    Ez, Hx, Hy, iEz, C1E, C2E, RXn, RYn,
    iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0):

    for i in prange(iMin, iMax + 1):
        for j in range(jMin, jMax + 1):
            for k in range(kMin, kMax + 0):
                n = (Ni * i) + (Nj * j) + (Nk * k) + N0
                m = iEz[n]
                Ez[n] = C1E[m] * Ez[n] \
                      + C2E[m] * (RXn[i] * (Hy[n] - Hy[n - Ni])
                                - RYn[j] * (Hx[n] - Hx[n - Nj]))

# 平面波入射
@jit(cache=True, nogil=True, parallel=True, nopython=True)
def _Ez_p_vector(
    t, ez, fPlanewave, Zc, Xn, Yn,
    Ez, Hx, Hy, iEz, K1Ez, K2Ez, RXn, RYn,
    iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0):

    for i in prange(iMin, iMax + 1):
        for j in range(jMin, jMax + 1):
            for k in range(kMin, kMax + 0):
                n = (Ni * i) + (Nj * j) + (Nk * k) + N0
                m = iEz[n]
                if m == 0:
                    Ez[n] += \
                          + RXn[i] * (Hy[n] - Hy[n - Ni]) \
                          - RYn[j] * (Hx[n] - Hx[n - Nj])
                else:
                    fi, dfi = sol.planewave.f(Xn[i], Yn[j], Zc[k], t, ez, fPlanewave)
                    Ez[n] = K1Ez[n] * Ez[n] \
                          + K2Ez[n] * (RXn[i] * (Hy[n] - Hy[n - Ni])
                                     - RYn[j] * (Hx[n] - Hx[n - Nj])) \
                          - (K1Ez[n] - K2Ez[n]) * dfi \
                          - (1 - K1Ez[n]) * fi

@jit(cache=True, nogil=True, parallel=True, nopython=True)
def _Ez_p_no_vector(
    t, ez, fPlanewave, Zc, Xn, Yn,
    Ez, Hx, Hy, iEz, C1E, C2E, RXn, RYn,
    iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0):

    for i in prange(iMin, iMax + 1):
        for j in range(jMin, jMax + 1):
            for k in range(kMin, kMax + 0):
                n = (Ni * i) + (Nj * j) + (Nk * k) + N0
                m = iEz[n]
                if m == 0:
                    Ez[n] += \
                          + RXn[i] * (Hy[n] - Hy[n - Ni]) \
                          - RYn[j] * (Hx[n] - Hx[n - Nj])
                else:
                    fi, dfi = sol.planewave.f(Xn[i], Yn[j], Zc[k], t, ez, fPlanewave)
                    Ez[n] = C1E[m] * Ez[n] \
                          + C2E[m] * (RXn[i] * (Hy[n] - Hy[n - Ni])
                                    - RYn[j] * (Hx[n] - Hx[n - Nj])) \
                          - (C1E[m] - C2E[m]) * dfi \
                          - (1 - C1E[m]) * fi
