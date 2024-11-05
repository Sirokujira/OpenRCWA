# -*- coding: utf-8 -*-
"""
Near3d.py
"""

#import numpy as np
from numba import jit, prange

# 全領域の電磁界のDFT(毎時刻ごと、1周波数で計算時間の4割を占める)
@jit(cache=True, nogil=True, parallel=True, nopython=True)
def dft(
    itime, Freq2, cEdft, cHdft,
    Ex, Ey, Ez, Hx, Hy, Hz,
    cEx, cEy, cEz, cHx, cHy, cHz,
    iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0, NN):

    for ifreq in range(len(Freq2)):

        n0 = ifreq * NN
        idx = (itime * len(Freq2)) + ifreq
        ef = cEdft[idx]
        hf = cHdft[idx]

        # Ex
        for i in prange(iMin - 0, iMax + 0):
            for j in range(jMin - 0, jMax + 1):
                n = (Ni * i) + (Nj * j) + (Nk * (kMin - 0)) + N0
                for _ in range(kMin - 0, kMax + 1):
                    #n = (Ni * i) + (Nj * j) + (Nk * k) + N0
                    cEx[n0 + n] += Ex[n] * ef
                    n += 1

        # Ey
        for i in prange(iMin - 0, iMax + 1):
            for j in range(jMin - 0, jMax + 0):
                n = (Ni * i) + (Nj * j) + (Nk * (kMin - 0)) + N0
                for _ in range(kMin - 0, kMax + 1):
                    #n = (Ni * i) + (Nj * j) + (Nk * k) + N0
                    cEy[n0 + n] += Ey[n] * ef
                    n += 1

        # Ez
        for i in prange(iMin - 0, iMax + 1):
            for j in range(jMin - 0, jMax + 1):
                n = (Ni * i) + (Nj * j) + (Nk * (kMin - 0)) + N0
                for _ in range(kMin - 0, kMax + 0):
                    #n = (Ni * i) + (Nj * j) + (Nk * k) + N0
                    cEz[n0 + n] += Ez[n] * ef
                    n += 1

        # Hx
        for i in prange(iMin - 1, iMax + 2):
            for j in range(jMin - 1, jMax + 1):
                n = (Ni * i) + (Nj * j) + (Nk * (kMin - 1)) + N0
                for _ in range(kMin - 1, kMax + 1):
                    #n = (Ni * i) + (Nj * j) + (Nk * k) + N0
                    cHx[n0 + n] += Hx[n] * hf
                    n += 1

        # Hy
        for i in prange(iMin - 1, iMax + 1):
            for j in range(jMin - 1, jMax + 2):
                n = (Ni * i) + (Nj * j) + (Nk * (kMin - 1)) + N0
                for _ in range(kMin - 1, kMax + 1):
                    #n = (Ni * i) + (Nj * j) + (Nk * k) + N0
                    cHy[n0 + n] += Hy[n] * hf
                    n += 1

        # Hz
        for i in prange(iMin - 1, iMax + 1):
            for j in range(jMin - 1, jMax + 1):
                n = (Ni * i) + (Nj * j) + (Nk * (kMin - 1)) + N0
                for _ in range(kMin - 1, kMax + 2):
                    #n = (Ni * i) + (Nj * j) + (Nk * k) + N0
                    cHz[n0 + n] += Hz[n] * hf
                    n += 1
