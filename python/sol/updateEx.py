# -*- coding: utf-8 -*-
"""
updateEx.py
"""

from numba import jit, prange
import sol.planewave

# Ex更新
def calEx(
    Parm, t, fPlanewave, Xc, Yn, Zn, Ex, Hy, Hz, iEx,
    C1E, C2E, K1Ex, K2Ex, RYn, RZn,
    iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0):

    if   Parm['source'] == 0:
        # 給電点
        if Parm['vector']:
            _Ex_f_vector(
                Ex, Hy, Hz, K1Ex, K2Ex, RYn, RZn,
                iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)
        else:
            _Ex_f_no_vector(
                Ex, Hy, Hz, iEx, C1E, C2E, RYn, RZn,
                iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)
    elif Parm['source'] == 1:
        # 平面波入射
        if Parm['vector']:
            _Ex_p_vector(
                t, fPlanewave[0], fPlanewave, Xc, Yn, Zn,
                Ex, Hy, Hz, iEx, K1Ex, K2Ex, RYn, RZn,
                iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)
        else:
            _Ex_p_no_vector(
                t, fPlanewave[0], fPlanewave, Xc, Yn, Zn,
                Ex, Hy, Hz, iEx, C1E, C2E, RYn, RZn,
                iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)

# 給電点
@jit(cache=True, nogil=True, parallel=True, nopython=True)
def _Ex_f_vector(
    Ex, Hy, Hz, K1Ex, K2Ex, RYn, RZn,
    iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0):

    for i in prange(iMin, iMax + 0):
        for j in range(jMin, jMax + 1):
            for k in range(kMin, kMax + 1):
                n = (Ni * i) + (Nj * j) + (Nk * k) + N0
                Ex[n] = K1Ex[n] * Ex[n] \
                      + K2Ex[n] * (RYn[j] * (Hz[n] - Hz[n - Nj])
                                 - RZn[k] * (Hy[n] - Hy[n - Nk]))

@jit(cache=True, nogil=True, parallel=True, nopython=True)
def _Ex_f_no_vector(
    Ex, Hy, Hz, iEx, C1E, C2E, RYn, RZn,
    iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0):

    for i in prange(iMin, iMax + 0):
        for j in range(jMin, jMax + 1):
            for k in range(kMin, kMax + 1):
                n = (Ni * i) + (Nj * j) + (Nk * k) + N0
                m = iEx[n]
                Ex[n] = C1E[m] * Ex[n] \
                      + C2E[m] * (RYn[j] * (Hz[n] - Hz[n - Nj])
                                - RZn[k] * (Hy[n] - Hy[n - Nk]))

# 平面波入射
@jit(cache=True, nogil=True, parallel=True, nopython=True)
def _Ex_p_vector(
    t, ex, fPlanewave, Xc, Yn, Zn,
    Ex, Hy, Hz, iEx, K1Ex, K2Ex, RYn, RZn,
    iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0):

    for i in prange(iMin, iMax + 0):
        for j in range(jMin, jMax + 1):
            for k in range(kMin, kMax + 1):
                n = (Ni * i) + (Nj * j) + (Nk * k) + N0
                m = iEx[n]
                if m == 0:
                    Ex[n] += \
                          + RYn[j] * (Hz[n] - Hz[n - Nj]) \
                          - RZn[k] * (Hy[n] - Hy[n - Nk])
                else:
                    fi, dfi = sol.planewave.f(Xc[i], Yn[j], Zn[k], t, ex, fPlanewave)
                    Ex[n] = K1Ex[n] * Ex[n] \
                          + K2Ex[n] * (RYn[j] * (Hz[n] - Hz[n - Nj])
                                     - RZn[k] * (Hy[n] - Hy[n - Nk])) \
                          - (K1Ex[n] - K2Ex[n]) * dfi \
                          - (1 - K1Ex[n]) * fi

@jit(cache=True, nogil=True, parallel=True, nopython=True)
def _Ex_p_no_vector(
    t, ex, fPlanewave, Xc, Yn, Zn,
    Ex, Hy, Hz, iEx, C1E, C2E, RYn, RZn,
    iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0):

    for i in prange(iMin, iMax + 0):
        for j in range(jMin, jMax + 1):
            for k in range(kMin, kMax + 1):
                n = (Ni * i) + (Nj * j) + (Nk * k) + N0
                m = iEx[n]
                if m == 0:
                    Ex[n] += \
                          + RYn[j] * (Hz[n] - Hz[n - Nj]) \
                          - RZn[k] * (Hy[n] - Hy[n - Nk])
                else:
                    fi, dfi = sol.planewave.f(Xc[i], Yn[j], Zn[k], t, ex, fPlanewave)
                    Ex[n] = C1E[m] * Ex[n] \
                          + C2E[m] * (RYn[j] * (Hz[n] - Hz[n - Nj])
                                    - RZn[k] * (Hy[n] - Hy[n - Nk])) \
                          - (C1E[m] - C2E[m]) * dfi \
                          - (1 - C1E[m]) * fi
