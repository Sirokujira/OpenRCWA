# -*- coding: utf-8 -*-
"""
PmlEy.py
"""

from numba import jit, prange

@jit(cache=True, nogil=True, parallel=True, nopython=True)
def calEy(
    Nz, Nx, Ey, Hz, Hx, Eyz, Eyx, RZn, RXn,
    iPmlEy, gPmlZn, gPmlXn, rPmlE, rPml, pml_l,
    Ni, Nj, Nk, N0):

    for n in prange(iPmlEy.shape[0]):
        i, j, k, m = iPmlEy[n]

        n0   = Ni * (i    ) + Nj * (j    ) + Nk * (k    ) + N0
        nk1  = Ni * (i    ) + Nj * (j    ) + Nk * (k - 1) + N0
        ni1  = Ni * (i - 1) + Nj * (j    ) + Nk * (k    ) + N0

        kc = min(max(k, 0), Nz    )
        dhx = Hx[n0] - Hx[nk1]
        rz = RZn[kc] * rPmlE[m]
        Eyz[n] = (Eyz[n] + (rz * dhx)) / (1 + (gPmlZn[k + pml_l] * rPml[m]))

        ic = min(max(i, 0), Nx    )
        dhz = Hz[n0] - Hz[ni1]
        rx = RXn[ic] * rPmlE[m]
        Eyx[n] = (Eyx[n] - (rx * dhz)) / (1 + (gPmlXn[i + pml_l] * rPml[m]))

        Ey[n0] = Eyz[n] + Eyx[n]

"""
# 配列演算(数十倍遅い)
def calEy( \
    Nz, Nx, Ey, Hz, Hx, Eyz, Eyx, RZn, RXn,
    iPmlEy, gPmlZn, gPmlXn, rPmlE, rPml, pml_l,
    Ni, Nj, Nk, N0):
    # [:]が必要

    i = iPmlEy[:, 0]
    j = iPmlEy[:, 1]
    k = iPmlEy[:, 2]
    m = iPmlEy[:, 3]

    n0   = Ni * (i    ) + Nj * (j    ) + Nk * (k    ) + N0
    nk1  = Ni * (i    ) + Nj * (j    ) + Nk * (k - 1) + N0
    ni1  = Ni * (i - 1) + Nj * (j    ) + Nk * (k    ) + N0

    kc = np.minimum(np.maximum(k, 0), Nz    )
    dhx = Hx[n0] - Hx[nk1]
    rz = RZn[kc] * rPmlE[m]
    Eyz[:] = (Eyz[:] + (rz * dhx)) / (1 + (gPmlZn[k + pml_l] * rPml[m]))

    ic = np.minimum(np.maximum(i, 0), Nx    )
    dhz = Hz[n0] - Hz[ni1]
    rx = RXn[ic] * rPmlE[m]
    Eyx[:] = (Eyx[:] - (rx * dhz)) / (1 + (gPmlXn[i + pml_l] * rPml[m]))

    Ey[n0] = Eyz[:] + Eyx[:]
"""