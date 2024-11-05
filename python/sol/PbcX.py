# -*- coding: utf-8 -*-
"""
PbcX.py
X方向周期境界条件
"""

from numba import jit

@jit(cache=True, nopython=True)
def x(Nx, Hy, Hz, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0):

    id1 = -1
    id2 = 0
    id3 = Nx - 1
    id4 = Nx

    # Hy
    for j in range(jMin - 0, jMax + 1):
        for k in range(kMin - 1, kMax + 1):
            n1 = (Ni * id1) + (Nj * j) + (Nk * k) + N0
            n2 = (Ni * id2) + (Nj * j) + (Nk * k) + N0
            n3 = (Ni * id3) + (Nj * j) + (Nk * k) + N0
            n4 = (Ni * id4) + (Nj * j) + (Nk * k) + N0
            Hy[n1] = Hy[n3]
            Hy[n4] = Hy[n2]

    # Hz
    for j in range(jMin - 1, jMax + 1):
        for k in range(kMin - 0, kMax + 1):
            n1 = (Ni * id1) + (Nj * j) + (Nk * k) + N0
            n2 = (Ni * id2) + (Nj * j) + (Nk * k) + N0
            n3 = (Ni * id3) + (Nj * j) + (Nk * k) + N0
            n4 = (Ni * id4) + (Nj * j) + (Nk * k) + N0
            Hz[n1] = Hz[n3]
            Hz[n4] = Hz[n2]
