# -*- coding: utf-8 -*-
"""
PmlHz.py
"""

from numba import jit, prange

@jit(cache=True, nogil=True, parallel=True, nopython=True)
def calHz(
    Nx, Ny, Hz, Ex, Ey, Hzx, Hzy, RXc, RYc,
    iPmlHz, gPmlXc, gPmlYc, rPmlH, rPml, pml_l,
    Ni, Nj, Nk, N0):

    for n in prange(iPmlHz.shape[0]):
        i, j, k, m = iPmlHz[n]

        n0   = Ni * (i    ) + Nj * (j    ) + Nk * (k    ) + N0
        ni1  = Ni * (i + 1) + Nj * (j    ) + Nk * (k    ) + N0
        nj1  = Ni * (i    ) + Nj * (j + 1) + Nk * (k    ) + N0

        ic = min(max(i, 0), Nx - 1)
        dey = Ey[ni1] - Ey[n0]
        rx = RXc[ic] * rPmlH[m]
        Hzx[n] = (Hzx[n] - (rx * dey)) / (1 + (gPmlXc[i + pml_l] * rPml[m]))

        jc = min(max(j, 0), Ny - 1)
        dex = Ex[nj1] - Ex[n0]
        ry = RYc[jc] * rPmlH[m]
        Hzy[n] = (Hzy[n] + (ry * dex)) / (1 + (gPmlYc[j + pml_l] * rPml[m]))

        Hz[n0] = Hzx[n] + Hzy[n]

"""
# 配列演算(数十倍遅い)
def calHz( \
    Nx, Ny, Hz, Ex, Ey, Hzx, Hzy, RXc, RYc,
    iPmlHz, gPmlXc, gPmlYc, rPmlH, rPml, pml_l,
    Ni, Nj, Nk, N0):
    # [:]が必要

    i = iPmlHz[:, 0]
    j = iPmlHz[:, 1]
    k = iPmlHz[:, 2]
    m = iPmlHz[:, 3]

    n0   = Ni * (i    ) + Nj * (j    ) + Nk * (k    ) + N0
    ni1  = Ni * (i + 1) + Nj * (j    ) + Nk * (k    ) + N0
    nj1  = Ni * (i    ) + Nj * (j + 1) + Nk * (k    ) + N0

    ic = np.minimum(np.maximum(i, 0), Nx - 1)
    dey = Ey[ni1] - Ey[n0]
    rx = RXc[ic] * rPmlH[m]
    Hzx[:] = (Hzx[:] - (rx * dey)) / (1 + (gPmlXc[i + pml_l] * rPml[m]))

    jc = np.minimum(np.maximum(j, 0), Ny - 1)
    dex = Ex[nj1] - Ex[n0]
    ry = RYc[jc] * rPmlH[m]
    Hzy[:] = (Hzy[:] + (ry * dex)) / (1 + (gPmlYc[j + pml_l] * rPml[m]))

    Hz[n0] = Hzx[:] + Hzy[:]
"""