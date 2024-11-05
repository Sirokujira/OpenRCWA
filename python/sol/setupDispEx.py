# -*- coding: utf-8 -*-
"""
setupDispEx.py
分散性媒質の準備Ex
"""

import numpy as np
from numba import jit

def setEx(
    iMaterial, fMaterial, iEx,
    iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0):

    # ダミー配列を用意する(jitに必要)
    iDispEx = np.zeros((0, 3), 'i4')
    fDispEx = np.zeros((0, 4), 'f8')

    # データ数を取得する
    num = _setEx(0,
        iMaterial, fMaterial, iEx, iDispEx, fDispEx,
        iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)
    #print(num)

    # 配列を確保する
    iDispEx = np.zeros((num, 3), 'i4')  # (i, j, k)
    fDispEx = np.zeros((num, 4), 'f8')  # e, factor

    # 変数に値を代入する
    _setEx(1,
        iMaterial, fMaterial, iEx, iDispEx, fDispEx,
        iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)

    return iDispEx, fDispEx

# (private)
@jit(cache=True, nopython=True)
def _setEx(mode,
    iMaterial, fMaterial, iEx, iDispEx, fDispEx,
    iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0):

    n = 0
    for i in range(iMin, iMax + 0):
        for j in range(jMin, jMax + 1):
            for k in range(kMin, kMax + 1):
                m = iEx[(Ni * i) + (Nj * j) + (Nk * k) + N0]
                if iMaterial[m] == 2:
                    # 分散性媒質
                    if mode == 1:
                        iDispEx[n, 0] = i
                        iDispEx[n, 1] = j
                        iDispEx[n, 2] = k
                        fDispEx[n, 1] = fMaterial[m,  8]
                        fDispEx[n, 2] = fMaterial[m,  9]
                        fDispEx[n, 3] = fMaterial[m, 10]
                    n += 1
    if mode == 0:
        return n
