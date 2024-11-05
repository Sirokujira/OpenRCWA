# -*- coding: utf-8 -*-
"""
PbcZ.py
Z方向周期境界条件
"""

from numba import jit

@jit(cache=True, nopython=True)
def z(Nz, Hx, Hy, iMin, iMax, jMin, jMax, Ni, Nj, Nk, N0):

    id1 = -1
    id2 = 0
    id3 = Nz - 1
    id4 = Nz

    # Hx
    for i in range(iMin - 0, iMax + 1):
        for j in range(jMin - 1, jMax + 1):
            n1 = (Ni * i) + (Nj * j) + (Nk * id1) + N0
            n2 = (Ni * i) + (Nj * j) + (Nk * id2) + N0
            n3 = (Ni * i) + (Nj * j) + (Nk * id3) + N0
            n4 = (Ni * i) + (Nj * j) + (Nk * id4) + N0
            Hx[n1] = Hx[n3]
            Hx[n4] = Hx[n2]

    # Hy
    for i in range(iMin - 1, iMax + 1):
        for j in range(jMin - 0, jMax + 1):
            n1 = (Ni * i) + (Nj * j) + (Nk * id1) + N0
            n2 = (Ni * i) + (Nj * j) + (Nk * id2) + N0
            n3 = (Ni * i) + (Nj * j) + (Nk * id3) + N0
            n4 = (Ni * i) + (Nj * j) + (Nk * id4) + N0
            Hy[n1] = Hy[n3]
            Hy[n4] = Hy[n2]
