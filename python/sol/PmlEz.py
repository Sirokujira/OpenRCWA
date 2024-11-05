# -*- coding: utf-8 -*-
"""
PmlEz.py
"""

from numba import jit, prange

@jit(cache=True, nogil=True, parallel=True, nopython=True)
def calEz(
    Nx, Ny, Ez, Hx, Hy, Ezx, Ezy, RXn, RYn,
    iPmlEz, gPmlXn, gPmlYn, rPmlE, rPml, pml_l,
    Ni, Nj, Nk, N0):

    for n in prange(iPmlEz.shape[0]):
        i, j, k, m = iPmlEz[n]

        n0   = Ni * (i    ) + Nj * (j    ) + Nk * (k    ) + N0
        ni1  = Ni * (i - 1) + Nj * (j    ) + Nk * (k    ) + N0
        nj1  = Ni * (i    ) + Nj * (j - 1) + Nk * (k    ) + N0

        ic = min(max(i, 0), Nx    )
        dhy = Hy[n0] - Hy[ni1]
        rx = RXn[ic] * rPmlE[m]
        Ezx[n] = (Ezx[n] + (rx * dhy)) / (1 + (gPmlXn[i + pml_l] * rPml[m]))

        jc = min(max(j, 0), Ny    )
        dhx = Hx[n0] - Hx[nj1]
        ry = RYn[jc] * rPmlE[m]
        Ezy[n] = (Ezy[n] - (ry * dhx)) / (1 + (gPmlYn[j + pml_l] * rPml[m]))

        Ez[n0] = Ezx[n] + Ezy[n]

"""
# 配列演算(数十倍遅い)
def calEz( \
    Nx, Ny, Ez, Hx, Hy, Ezx, Ezy, RXn, RYn,
    iPmlEz, gPmlXn, gPmlYn, rPmlE, rPml, pml_l,
    Ni, Nj, Nk, N0):
    # [:]が必要

    i = iPmlEz[:, 0]
    j = iPmlEz[:, 1]
    k = iPmlEz[:, 2]
    m = iPmlEz[:, 3]

    n0   = Ni * (i    ) + Nj * (j    ) + Nk * (k    ) + N0
    ni1  = Ni * (i - 1) + Nj * (j    ) + Nk * (k    ) + N0
    nj1  = Ni * (i    ) + Nj * (j - 1) + Nk * (k    ) + N0

    ic = np.minimum(np.maximum(i, 0), Nx    )
    dhy = Hy[n0] - Hy[ni1]
    rx = RXn[ic] * rPmlE[m]
    Ezx[:] = (Ezx[:] + (rx * dhy)) / (1 + (gPmlXn[i + pml_l] * rPml[m]))

    jc = np.minimum(np.maximum(j, 0), Ny    )
    dhx = Hx[n0] - Hx[nj1]
    ry = RYn[jc] * rPmlE[m]
    Ezy[:] = (Ezy[:] - (ry * dhx)) / (1 + (gPmlYn[j + pml_l] * rPml[m]))

    Ez[n0] = Ezx[:] + Ezy[:]
"""