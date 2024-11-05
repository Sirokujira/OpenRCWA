# -*- coding: utf-8 -*-
"""
updateHz.py
"""

from numba import jit, prange
import sol.planewave

# Hz更新
def calHz(
    Parm, t, fPlanewave, Zn, Xc, Yc, Hz, Ex, Ey, iHz,
    C1H, C2H, K1Hz, K2Hz, RXc, RYc,
    iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0):

    if   Parm['source'] == 0:
        # 給電点
        if Parm['vector']:
            _Hz_f_vector(
                Hz, Ex, Ey, K1Hz, K2Hz, RXc, RYc,
                iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)
        else:
            _Hz_f_no_vector(
                Hz, Ex, Ey, iHz, C1H, C2H, RXc, RYc,
                iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)
    elif Parm['source'] == 1:
        # 平面波入射
        if Parm['vector']:
            _Hz_p_vector(
                t, fPlanewave[5], fPlanewave, Zn, Xc, Yc,
                Hz, Ex, Ey, iHz, K1Hz, K2Hz, RXc, RYc,
                iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)
        else:
            _Hz_p_no_vector(
                t, fPlanewave[5], fPlanewave, Zn, Xc, Yc,
                Hz, Ex, Ey, iHz, C1H, C2H, RXc, RYc,
                iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)

# 給電点
@jit(cache=True, nogil=True, parallel=True, nopython=True)
def _Hz_f_vector(
    Hz, Ex, Ey, K1Hz, K2Hz, RXc, RYc,
    iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0):

    for i in prange(iMin, iMax + 0):
        for j in range(jMin, jMax + 0):
            for k in range(kMin, kMax + 1):
                n = (Ni * i) + (Nj * j) + (Nk * k) + N0
                Hz[n] = K1Hz[n] * Hz[n] \
                      - K2Hz[n] * (RXc[i] * (Ey[n + Ni] - Ey[n])
                                 - RYc[j] * (Ex[n + Nj] - Ex[n]))

@jit(cache=True, nogil=True, parallel=True, nopython=True)
def _Hz_f_no_vector(
    Hz, Ex, Ey, iHz, C1H, C2H, RXc, RYc,
    iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0):

    for i in prange(iMin, iMax + 0):
        for j in range(jMin, jMax + 0):
            for k in range(kMin, kMax + 1):
                n = (Ni * i) + (Nj * j) + (Nk * k) + N0
                m = iHz[n]
                Hz[n] = C1H[m] * Hz[n] \
                      - C2H[m] * (RXc[i] * (Ey[n + Ni] - Ey[n])
                                - RYc[j] * (Ex[n + Nj] - Ex[n]))

# 平面波入射
@jit(cache=True, nogil=True, parallel=True, nopython=True)
def _Hz_p_vector(
    t, hz, fPlanewave, Zn, Xc, Yc,
    Hz, Ex, Ey, iHz, K1Hz, K2Hz, RXc, RYc,
    iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0):

    for i in prange(iMin, iMax + 0):
        for j in range(jMin, jMax + 0):
            for k in range(kMin, kMax + 1):
                n = (Ni * i) + (Nj * j) + (Nk * k) + N0
                m = iHz[n]
                if m == 0:
                    Hz[n] += \
                          - RXc[i] * (Ey[n + Ni] - Ey[n]) \
                          + RYc[j] * (Ex[n + Nj] - Ex[n])
                else:
                    fi, dfi = sol.planewave.f(Xc[i], Yc[j], Zn[k], t, hz, fPlanewave)
                    Hz[n] = K1Hz[n] * Hz[n] \
                          - K2Hz[n] * (RXc[i] * (Ey[n + Ni] - Ey[n])
                                     - RYc[j] * (Ex[n + Nj] - Ex[n])) \
                          - (K1Hz[n] - K2Hz[n]) * dfi \
                          - (1 - K1Hz[n]) * fi

@jit(cache=True, nogil=True, parallel=True, nopython=True)
def _Hz_p_no_vector(
    t, hz, fPlanewave, Zn, Xc, Yc,
    Hz, Ex, Ey, iHz, C1H, C2H, RXc, RYc,
    iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0):

    for i in prange(iMin, iMax + 0):
        for j in range(jMin, jMax + 0):
            for k in range(kMin, kMax + 1):
                n = (Ni * i) + (Nj * j) + (Nk * k) + N0
                m = iHz[n]
                if m == 0:
                    Hz[n] += \
                          - RXc[i] * (Ey[n + Ni] - Ey[n]) \
                          + RYc[j] * (Ex[n + Nj] - Ex[n])
                else:
                    fi, dfi = sol.planewave.f(Xc[i], Yc[j], Zn[k], t, hz, fPlanewave)
                    Hz[n] = C1H[m] * Hz[n] \
                          - C2H[m] * (RXc[i] * (Ey[n + Ni] - Ey[n])
                                    - RYc[j] * (Ex[n + Nj] - Ex[n])) \
                          - (C1H[m] - C2H[m]) * dfi \
                          - (1 - C1H[m]) * fi
