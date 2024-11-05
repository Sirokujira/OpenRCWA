# -*- coding: utf-8 -*-
"""
PmlEx.py
"""

from numba import jit, prange

@jit(cache=True, nogil=True, parallel=True, nopython=True)
def calEx(
    Ny, Nz, Ex, Hy, Hz, Exy, Exz, RYn, RZn,
    iPmlEx, gPmlYn, gPmlZn, rPmlE, rPml, pml_l,
    Ni, Nj, Nk, N0):

    for n in prange(iPmlEx.shape[0]):
        i, j, k, m = iPmlEx[n]

        n0   = Ni * (i    ) + Nj * (j    ) + Nk * (k    ) + N0
        nj1  = Ni * (i    ) + Nj * (j - 1) + Nk * (k    ) + N0
        nk1  = Ni * (i    ) + Nj * (j    ) + Nk * (k - 1) + N0

        jc = min(max(j, 0), Ny    )
        dhz = Hz[n0] - Hz[nj1]
        ry = RYn[jc] * rPmlE[m]
        Exy[n] = (Exy[n] + (ry * dhz)) / (1 + (gPmlYn[j + pml_l] * rPml[m]))

        kc = min(max(k, 0), Nz    )
        dhy = Hy[n0] - Hy[nk1]
        rz = RZn[kc] * rPmlE[m]
        Exz[n] = (Exz[n] - (rz * dhy)) / (1 + (gPmlZn[k + pml_l] * rPml[m]))

        Ex[n0] = Exy[n] + Exz[n]

"""
# 配列演算(数十倍遅い)
def calEx( \
    Ny, Nz, Ex, Hy, Hz, Exy, Exz, RYn, RZn,
    iPmlEx, gPmlYn, gPmlZn, rPmlE, rPml, pml_l,
    Ni, Nj, Nk, N0):
    # [:]が必要

    i = iPmlEx[:, 0]
    j = iPmlEx[:, 1]
    k = iPmlEx[:, 2]
    m = iPmlEx[:, 3]

    n0   = Ni * (i    ) + Nj * (j    ) + Nk * (k    ) + N0
    nj1  = Ni * (i    ) + Nj * (j - 1) + Nk * (k    ) + N0
    nk1  = Ni * (i    ) + Nj * (j    ) + Nk * (k - 1) + N0

    jc = np.minimum(np.maximum(j, 0), Ny    )
    dhz = Hz[n0] - Hz[nj1]
    ry = RYn[jc] * rPmlE[m]
    Exy[:] = (Exy[:] + (ry * dhz)) / (1 + (gPmlYn[j + pml_l] * rPml[m]))

    kc = np.minimum(np.maximum(k, 0), Nz    )
    dhy = Hy[n0] - Hy[nk1]
    rz = RZn[kc] * rPmlE[m]
    Exz[:] = (Exz[:] - (rz * dhy)) / (1 + (gPmlZn[k + pml_l] * rPml[m]))

    Ex[n0] = Exy[:] + Exz[:]
"""