# -*- coding: utf-8 -*-
"""
setupDispEz.py
分散性媒質の準備Ez
"""

import numpy as np
from numba import jit

def setEz(
    iMaterial, fMaterial, iEz,
    iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0):

    # ダミー配列を用意する(jitに必要)
    iDispEz = np.zeros((0, 3), 'i4')
    fDispEz = np.zeros((0, 4), 'f8')

    # データ数を取得する
    num = _setEz(0,
        iMaterial, fMaterial, iEz, iDispEz, fDispEz,
        iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)
    #print(num)

    # 配列を確保する
    iDispEz = np.zeros((num, 3), 'i4')  # (i, j, k)
    fDispEz = np.zeros((num, 4), 'f8')  # e, factor

    # 変数に値を代入する
    _setEz(1,
        iMaterial, fMaterial, iEz, iDispEz, fDispEz,
        iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)

    return iDispEz, fDispEz

# (private)
@jit(cache=True, nopython=True)
def _setEz(mode,
    iMaterial, fMaterial, iEz, iDispEz, fDispEz,
    iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0):

    n = 0
    for i in range(iMin, iMax + 1):
        for j in range(jMin, jMax + 1):
            for k in range(kMin, kMax + 0):
                m = iEz[(Ni * i) + (Nj * j) + (Nk * k) + N0]
                if iMaterial[m] == 2:
                    # 分散性媒質
                    if mode == 1:
                        iDispEz[n, 0] = i
                        iDispEz[n, 1] = j
                        iDispEz[n, 2] = k
                        fDispEz[n, 1] = fMaterial[m,  8]
                        fDispEz[n, 2] = fMaterial[m,  9]
                        fDispEz[n, 3] = fMaterial[m, 10]
                    n += 1
    if mode == 0:
        return n
