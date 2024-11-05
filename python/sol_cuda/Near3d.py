# -*- coding: utf-8 -*-
"""
Near3d.py (CUDA)
"""

import math
from numba import cuda

# 全領域の電磁界のDFT(毎時刻ごと、1周波数で計算時間の4割を占める)
def dft(block3d,
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
        gridEx = (math.ceil((kMax - kMin + 1) / block3d[0]),
                  math.ceil((jMax - jMin + 1) / block3d[1]),
                  math.ceil((iMax - iMin + 0) / block3d[2]))
        _dftEx_gpu[gridEx, block3d](
            Ex, cEx, ef, n0,
            iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)

        # Ey
        gridEy = (math.ceil((kMax - kMin + 1) / block3d[0]),
                  math.ceil((jMax - jMin + 0) / block3d[1]),
                  math.ceil((iMax - iMin + 1) / block3d[2]))
        _dftEy_gpu[gridEy, block3d](
            Ey, cEy, ef, n0,
            iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)

        # Ez
        gridEz = (math.ceil((kMax - kMin + 0) / block3d[0]),
                  math.ceil((jMax - jMin + 1) / block3d[1]),
                  math.ceil((iMax - iMin + 1) / block3d[2]))
        _dftEz_gpu[gridEz, block3d](
            Ez, cEz, ef, n0,
            iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)

        # Hx
        gridHx = (math.ceil((kMax - kMin + 2) / block3d[0]),
                  math.ceil((jMax - jMin + 2) / block3d[1]),
                  math.ceil((iMax - iMin + 1) / block3d[2]))
        _dftHx_gpu[gridHx, block3d](
            Hx, cHx, hf, n0,
            iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)

        # Hy
        gridHy = (math.ceil((kMax - kMin + 2) / block3d[0]),
                  math.ceil((jMax - jMin + 1) / block3d[1]),
                  math.ceil((iMax - iMin + 2) / block3d[2]))
        _dftHy_gpu[gridHy, block3d](
            Hy, cHy, hf, n0,
            iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)

        # Hz
        gridHz = (math.ceil((kMax - kMin + 1) / block3d[0]),
                  math.ceil((jMax - jMin + 2) / block3d[1]),
                  math.ceil((iMax - iMin + 2) / block3d[2]))
        _dftHz_gpu[gridHz, block3d](
            Hz, cHz, hf, n0,
            iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)

# Ex
@cuda.jit(cache=True)
def _dftEx_gpu(
    Ex, cEx, ef, n0,
    iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0):

    k, j, i = cuda.grid(3)
    i += iMin
    j += jMin
    k += kMin
    if (i <  iMax) and \
       (j <= jMax) and \
       (k <= kMax):
        n = (Ni * i) + (Nj * j) + (Nk * k) + N0
        cEx[n0 + n] += Ex[n] * ef
        #cEx[n0 + n] += Ex[n].real * ef + 1j * Ex[n].imag * ef

# Ey
@cuda.jit(cache=True)
def _dftEy_gpu(
    Ey, cEy, ef, n0,
    iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0):

    k, j, i = cuda.grid(3)
    i += iMin
    j += jMin
    k += kMin
    if (i <= iMax) and \
       (j <  jMax) and \
       (k <= kMax):
        n = (Ni * i) + (Nj * j) + (Nk * k) + N0
        cEy[n0 + n] += Ey[n] * ef
        #cEy[n0 + n] += Ey[n].real * ef + 1j * Ey[n].imag * ef

# Ez
@cuda.jit(cache=True)
def _dftEz_gpu(
    Ez, cEz, ef, n0,
    iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0):

    k, j, i = cuda.grid(3)
    i += iMin
    j += jMin
    k += kMin
    if (i <= iMax) and \
       (j <= jMax) and \
       (k <  kMax):
        n = (Ni * i) + (Nj * j) + (Nk * k) + N0
        cEz[n0 + n] += Ez[n] * ef
        #cEz[n0 + n] += Ez[n].real * ef + 1j * Ez[n].imag * ef

# Hx
@cuda.jit(cache=True)
def _dftHx_gpu(
    Hx, cHx, hf, n0,
    iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0):

    k, j, i = cuda.grid(3)
    i += iMin - 0
    j += jMin - 1
    k += kMin - 1
    if (i <= iMax) and \
       (j <= jMax) and \
       (k <= kMax):
        n = (Ni * i) + (Nj * j) + (Nk * k) + N0
        cHx[n0 + n] += Hx[n] * hf
        #cHx[n0 + n] += Hx[n].real * hf + 1j * Hx[n].imag * hf

# Hy
@cuda.jit(cache=True)
def _dftHy_gpu(
    Hy, cHy, hf, n0,
    iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0):

    k, j, i = cuda.grid(3)
    i += iMin - 1
    j += jMin - 0
    k += kMin - 1
    if (i <= iMax) and \
       (j <= jMax) and \
       (k <= kMax):
        n = (Ni * i) + (Nj * j) + (Nk * k) + N0
        cHy[n0 + n] += Hy[n] * hf
        #cHy[n0 + n] += Hy[n].real * hf + 1j * Hy[n].imag * hf

# Hz
@cuda.jit(cache=True)
def _dftHz_gpu(
    Hz, cHz, hf, n0,
    iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0):

    k, j, i = cuda.grid(3)
    i += iMin - 1
    j += jMin - 1
    k += kMin - 0
    if (i <= iMax) and \
       (j <= jMax) and \
       (k <= kMax):
        n = (Ni * i) + (Nj * j) + (Nk * k) + N0
        cHz[n0 + n] += Hz[n] * hf
        #cHz[n0 + n] += Hz[n].real * hf + 1j * Hz[n].imag * hf
