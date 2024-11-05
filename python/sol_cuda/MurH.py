# -*- coding: utf-8 -*-
"""
MurH.py (CUDA)
Mur-ABC (Hx/Hy/Hz共通)
"""

import math
from numba import cuda

def calcH(block1d, H, fMurH, iMurH, Ni, Nj, Nk, N0):
    numH = fMurH.shape[0]
    grid = math.ceil(numH / block1d)
    _calc_gpu[grid, block1d](H, numH, fMurH, iMurH, Ni, Nj, Nk, N0)

# (private)
@cuda.jit(cache=True)
def _calc_gpu(H, numH, fMurH, iMurH, Ni, Nj, Nk, N0):
    num = cuda.grid(1)
    if num < numH:
        i, j, k, i1, j1, k1 = iMurH[num]
        n  = (Ni * i ) + (Nj * j ) + (Nk * k ) + N0
        n1 = (Ni * i1) + (Nj * j1) + (Nk * k1) + N0
        H[n] = fMurH[num, 0] + fMurH[num, 1] * (H[n1] - H[n])
        fMurH[num, 0] = H[n1]
