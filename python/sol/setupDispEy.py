# -*- coding: utf-8 -*-
"""
setupDispEy.py
分散性媒質の準備Ey
"""

import numpy as np
from numba import jit

def setEy(
    iMaterial, fMaterial, iEy,
    iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0):

    # ダミー配列を用意する(jitに必要)
    iDispEy = np.zeros((0, 3), 'i4')
    fDispEy = np.zeros((0, 4), 'f8')

    # データ数を取得する
    num = _setEy(0,
        iMaterial, fMaterial, iEy, iDispEy, fDispEy,
        iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)
    #print(num)

    # 配列を確保する
    iDispEy = np.zeros((num, 3), 'i4')  # (i, j, k)
    fDispEy = np.zeros((num, 4), 'f8')  # e, factor

    # 変数に値を代入する
    _setEy(1,
        iMaterial, fMaterial, iEy, iDispEy, fDispEy,
        iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)

    return iDispEy, fDispEy

# (private)
@jit(cache=True, nopython=True)
def _setEy(mode,
    iMaterial, fMaterial, iEy, iDispEy, fDispEy,
    iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0):

    n = 0
    for i in range(iMin, iMax + 1):
        for j in range(jMin, jMax + 0):
            for k in range(kMin, kMax + 1):
                m = iEy[(Ni * i) + (Nj * j) + (Nk * k) + N0]
                if iMaterial[m] == 2:
                    # 分散性媒質
                    if mode == 1:
                        iDispEy[n, 0] = i
                        iDispEy[n, 1] = j
                        iDispEy[n, 2] = k
                        fDispEy[n, 1] = fMaterial[m,  8]
                        fDispEy[n, 2] = fMaterial[m,  9]
                        fDispEy[n, 3] = fMaterial[m, 10]
                    n += 1
    if mode == 0:
        return n
