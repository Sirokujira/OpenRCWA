# -*- coding: utf-8 -*-
"""
setupPmlHz.py
"""

from numba import jit

@jit(cache=True, nogil=True)#, parallel=True, nopython=True) # prange不可
def setHz(mode, pml_l, Nx, Ny, Nz, iHz, iPmlHz,
    iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0):

    lx = pml_l
    ly = pml_l
    lz = pml_l
    imin = iMin - lx + 0
    jmin = jMin - ly + 0
    kmin = kMin - lz + 1
    imax = iMax + lx
    jmax = jMax + ly
    kmax = kMax + lz

    num = 0
    for i in range(imin, imax):
        for j in range(jmin, jmax):
            for k in range(kmin, kmax):
                if (i < 0) or (Nx <= i) or \
                   (j < 0) or (Ny <= j) or \
                   (k < 0) or (Nz <  k):
                    if mode == 1:
                        iPmlHz[num, 0] = i
                        iPmlHz[num, 1] = j
                        iPmlHz[num, 2] = k
                        i_ = max(0, min(Nx - 1, i))
                        j_ = max(0, min(Ny - 1, j))
                        k_ = max(0, min(Nz,     k))
                        n = 0
                        if   k  <  0:
                            n = Ni * (i_    ) + Nj * (j_    ) + Nk * (0     ) + N0
                        elif Nz <  k:
                            n = Ni * (i_    ) + Nj * (j_    ) + Nk * (Nz    ) + N0
                        elif i  <  0:
                            n = Ni * (0     ) + Nj * (j_    ) + Nk * (k_    ) + N0
                        elif Nx <= i:
                            n = Ni * (Nx - 1) + Nj * (j_    ) + Nk * (k_    ) + N0
                        elif j  <  0:
                            n = Ni * (i_    ) + Nj * (0     ) + Nk * (k_    ) + N0
                        elif Ny <= j:
                            n = Ni * (i_    ) + Nj * (Ny - 1) + Nk * (k_    ) + N0
                        iPmlHz[num, 3] = iHz[n]
                    num += 1

    # 配列の大きさ
    if mode == 0:
        return num
