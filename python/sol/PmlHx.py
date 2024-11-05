# -*- coding: utf-8 -*-
"""
PmlHx.py
"""

from numba import jit, prange

@jit(cache=True, nogil=True, parallel=True, nopython=True)
def calHx(
    Ny, Nz, Hx, Ey, Ez, Hxy, Hxz, RYc, RZc,
    iPmlHx, gPmlYc, gPmlZc, rPmlH, rPml, pml_l,
    Ni, Nj, Nk, N0):

    for n in prange(iPmlHx.shape[0]):
        i, j, k, m = iPmlHx[n]

        n0   = Ni * (i    ) + Nj * (j    ) + Nk * (k    ) + N0
        nj1  = Ni * (i    ) + Nj * (j + 1) + Nk * (k    ) + N0
        nk1  = Ni * (i    ) + Nj * (j    ) + Nk * (k + 1) + N0

        jc = min(max(j, 0), Ny - 1)
        dez = Ez[nj1] - Ez[n0]
        ry = RYc[jc] * rPmlH[m]
        Hxy[n] = (Hxy[n] - (ry * dez)) / (1 + (gPmlYc[j + pml_l] * rPml[m]))

        kc = min(max(k, 0), Nz - 1)
        dey = Ey[nk1] - Ey[n0]
        rz = RZc[kc] * rPmlH[m]
        Hxz[n] = (Hxz[n] + (rz * dey)) / (1 + (gPmlZc[k + pml_l] * rPml[m]))

        Hx[n0] = Hxy[n] + Hxz[n]

"""
# 配列演算(数十倍遅い)
def calHx( \
    Ny, Nz, Hx, Ey, Ez, Hxy, Hxz, RYc, RZc,
    iPmlHx, gPmlYc, gPmlZc, rPmlH, rPml, pml_l,
    Ni, Nj, Nk, N0):
    # [:]が必要

    i = iPmlHx[:, 0]
    j = iPmlHx[:, 1]
    k = iPmlHx[:, 2]
    m = iPmlHx[:, 3]

    n0   = Ni * (i    ) + Nj * (j    ) + Nk * (k    ) + N0
    nj1  = Ni * (i    ) + Nj * (j + 1) + Nk * (k    ) + N0
    nk1  = Ni * (i    ) + Nj * (j    ) + Nk * (k + 1) + N0

    jc = np.minimum(np.maximum(j, 0), Ny - 1)
    dez = Ez[nj1] - Ez[n0]
    ry = RYc[jc] * rPmlH[m]
    Hxy[:] = (Hxy[:] - (ry * dez)) / (1 + (gPmlYc[j + pml_l] * rPml[m]))

    kc = np.minimum(np.maximum(k, 0), Nz - 1)
    dey = Ey[nk1] - Ey[n0]
    rz = RZc[kc] * rPmlH[m]
    Hxz[:] = (Hxz[:] + (rz * dey)) / (1 + (gPmlZc[k + pml_l] * rPml[m]))

    Hx[n0] = Hxy[:] + Hxz[:]
"""