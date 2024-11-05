# -*- coding: utf-8 -*-
"""
PbcY.py
Y方向周期境界条件
"""

from numba import jit

@jit(cache=True, nopython=True)
def y(Ny, Hz, Hx, kMin, kMax, iMin, iMax, Ni, Nj, Nk, N0):

    id1 = -1
    id2 = 0
    id3 = Ny - 1
    id4 = Ny

    # Hz
    for k in range(kMin - 0, kMax + 1):
        for i in range(iMin - 1, iMax + 1):
            n1 = (Ni * i) + (Nj * id1) + (Nk * k) + N0
            n2 = (Ni * i) + (Nj * id2) + (Nk * k) + N0
            n3 = (Ni * i) + (Nj * id3) + (Nk * k) + N0
            n4 = (Ni * i) + (Nj * id4) + (Nk * k) + N0
            Hz[n1] = Hz[n3]
            Hz[n4] = Hz[n2]

    # Hx
    for k in range(kMin - 1, kMax + 1):
        for i in range(iMin - 0, iMax + 1):
            n1 = (Ni * i) + (Nj * id1) + (Nk * k) + N0
            n2 = (Ni * i) + (Nj * id2) + (Nk * k) + N0
            n3 = (Ni * i) + (Nj * id3) + (Nk * k) + N0
            n4 = (Ni * i) + (Nj * id4) + (Nk * k) + N0
            Hx[n1] = Hx[n3]
            Hx[n4] = Hx[n2]
