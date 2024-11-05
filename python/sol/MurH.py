# -*- coding: utf-8 -*-
"""
MurH.py
"""

import math
from numba import jit, prange

# Hx/Hy/Hz共通
@jit(cache=True, nopython=True, nogil=True, parallel=True)
def calcH(h, fMurH, iMurH, Ni, Nj, Nk, N0):

    for m in prange(fMurH.shape[0]):
        i, j, k, i1, j1, k1 = iMurH[m]
        n  = (Ni * i ) + (Nj * j ) + (Nk * k ) + N0
        n1 = (Ni * i1) + (Nj * j1) + (Nk * k1) + N0
        h[n] = fMurH[m, 0] + fMurH[m, 1] * (h[n1] - h[n])
        fMurH[m, 0] = h[n1]

# Murの係数
@jit(cache=True, nopython=True)
def factor(fMaterial, d, m, cdt):

    PEC = 1
    if m != PEC:
        epsr = fMaterial[m, 0]
        amur = fMaterial[m, 2]
        vdt = cdt / math.sqrt(epsr * amur)
        return (vdt - d) / (vdt + d)
    else:
        return -1
