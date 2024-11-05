# -*- coding: utf-8 -*-
"""
PmlHy.py
"""

from numba import jit, prange

@jit(cache=True, nogil=True, parallel=True, nopython=True)
def calHy(
    Nz, Nx, Hy, Ez, Ex, Hyz, Hyx, RZc, RXc,
    iPmlHy, gPmlZc, gPmlXc, rPmlH, rPml, pml_l,
    Ni, Nj, Nk, N0):

    for n in prange(iPmlHy.shape[0]):
        i, j, k, m = iPmlHy[n]

        n0   = Ni * (i    ) + Nj * (j    ) + Nk * (k    ) + N0
        nk1  = Ni * (i    ) + Nj * (j    ) + Nk * (k + 1) + N0
        ni1  = Ni * (i + 1) + Nj * (j    ) + Nk * (k    ) + N0

        kc = min(max(k, 0), Nz - 1)
        dex = Ex[nk1] - Ex[n0]
        rz = RZc[kc] * rPmlH[m]
        Hyz[n] = (Hyz[n] - (rz * dex)) / (1 + (gPmlZc[k + pml_l] * rPml[m]))

        ic = min(max(i, 0), Nx - 1)
        dez = Ez[ni1] - Ez[n0]
        rx = RXc[ic] * rPmlH[m]
        Hyx[n] = (Hyx[n] + (rx * dez)) / (1 + (gPmlXc[i + pml_l] * rPml[m]))

        Hy[n0] = Hyz[n] + Hyx[n]

"""
# 配列演算(数十倍遅い)
def calHy( \
    Nz, Nx, Hy, Ez, Ex, Hyz, Hyx, RZc, RXc,
    iPmlHy, gPmlZc, gPmlXc, rPmlH, rPml, pml_l,
    Ni, Nj, Nk, N0):
    # [:]が必要

    i = iPmlHy[:, 0]
    j = iPmlHy[:, 1]
    k = iPmlHy[:, 2]
    m = iPmlHy[:, 3]

    n0   = Ni * (i    ) + Nj * (j    ) + Nk * (k    ) + N0
    nk1  = Ni * (i    ) + Nj * (j    ) + Nk * (k + 1) + N0
    ni1  = Ni * (i + 1) + Nj * (j    ) + Nk * (k    ) + N0

    kc = np.minimum(np.maximum(k, 0), Nz - 1)
    dex = Ex[nk1] - Ex[n0]
    rz = RZc[kc] * rPmlH[m]
    Hyz[:] = (Hyz[:] - (rz * dex)) / (1 + (gPmlZc[k + pml_l] * rPml[m]))

    ic = np.minimum(np.maximum(i, 0), Nx - 1)
    dez = Ez[ni1] - Ez[n0]
    rx = RXc[ic] * rPmlH[m]
    Hyx[:] = (Hyx[:] + (rx * dez)) / (1 + (gPmlXc[i + pml_l] * rPml[m]))

    Hy[n0] = Hyz[:] + Hyx[:]
"""