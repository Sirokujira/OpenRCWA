# -*- coding: utf-8 -*-
"""
cpmmon.py
"""

import math, cmath
import numpy as np
import sol.Fnode

# 指定した節点(i,j,k)の電界と磁界を求める
def fnode(
    ifreq, i, j, k,
    kwave, source, fPlanewave, noinc,
    Nx, Ny, Nz, Xn, Yn, Zn, cEx, cEy, cEz, cHx, cHy, cHz, Ni, Nj, Nk, N0, NN):

    # E
    cex, cey, cez = sol.Fnode.e(ifreq, i, j, k,
        Nx, Ny, Nz, cEx, cEy, cEz, Ni, Nj, Nk, N0, NN)

    # H
    chx, chy, chz = sol.Fnode.h(ifreq, i, j, k,
        cHx, cHy, cHz, Ni, Nj, Nk, N0, NN)

    # 平面波入射のとき入射界を加える
    if (source == 1) and (noinc == 0):
        ceinc, chinc = _planewave(i, j, k, fPlanewave, Xn, Yn, Zn, kwave)
        cex += ceinc[0]
        cey += ceinc[1]
        cez += ceinc[2]
        chx += chinc[0]
        chy += chinc[1]
        chz += chinc[2]

    # E/H 14成分
    f = np.zeros(14, 'f8')
    f[ 0] = math.sqrt(abs(cex)**2 + abs(cey)**2 + abs(cez)**2)
    f[ 1] = abs(cex)
    f[ 2] = abs(cey)
    f[ 3] = abs(cez)
    f[ 4] = math.degrees(cmath.phase(cex))
    f[ 5] = math.degrees(cmath.phase(cey))
    f[ 6] = math.degrees(cmath.phase(cez))
    f[ 7] = math.sqrt(abs(chx)**2 + abs(chy)**2 + abs(chz)**2)
    f[ 8] = abs(chx)
    f[ 9] = abs(chy)
    f[10] = abs(chz)
    f[11] = math.degrees(cmath.phase(chx))
    f[12] = math.degrees(cmath.phase(chy))
    f[13] = math.degrees(cmath.phase(chz))

    return f

# (private)
def _planewave(i, j, k, fPlanewave, Xn, Yn, Zn, kwave):

    x0 = (Xn[0] + Xn[-1]) / 2
    y0 = (Yn[0] + Yn[-1]) / 2
    z0 = (Zn[0] + Zn[-1]) / 2

    ei = fPlanewave[0:3]
    hi = fPlanewave[3:6]
    ri = fPlanewave[6:9]

    rri = (Xn[i] - x0) * ri[0] \
        + (Yn[j] - y0) * ri[1] \
        + (Zn[k] - z0) * ri[2]

    #phs = cmath.exp(complex(0, -kwave * rri))
    phs = math.cos(kwave * rri) - 1j * math.sin(kwave * rri)

    return ei * phs, hi * phs