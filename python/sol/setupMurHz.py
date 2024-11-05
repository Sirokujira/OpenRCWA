# -*- coding: utf-8 -*-
"""
setupMurHz.py
"""

from numba import jit
import sol.MurH

@jit(cache=True, nopython=True)
def setHz(mode,
    Nx, Ny, Xn, Yn, fMaterial, iEx, iEy, fMurHz, iMurHz, cdt,
    iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0):

    imin = iMin - 1
    jmin = jMin - 1
    kmin = kMin - 0
    imax = iMax + 1
    jmax = jMax + 1
    kmax = kMax + 1

    num = 0
    for i in range(imin, imax):
        for j in range(jmin, jmax):
            for k in range(kmin, kmax):
                if ((((i < 0) or (i >= Nx)) and (j >= 0) and (j < Ny)) or \
                    (((j < 0) or (j >= Ny)) and (i >= 0) and (i < Nx))):
                    if mode == 1:
                        iMurHz[num, 0] = i
                        iMurHz[num, 1] = j
                        iMurHz[num, 2] = k
                        m = d = i1 = j1 = k1 = 0
                        n   = Ni * (i    ) + Nj * (j    ) + Nk * (k    ) + N0
                        ni1 = Ni * (i + 1) + Nj * (j    ) + Nk * (k    ) + N0
                        nj1 = Ni * (i    ) + Nj * (j + 1) + Nk * (k    ) + N0
                        if   i <   0:
                            m = iEy[ni1]
                            d = Xn[1] - Xn[0]
                            i1 = i + 1
                            j1 = j
                            k1 = k
                        elif i >= Nx:
                            m = iEy[n]
                            d = Xn[-1] - Xn[-2]
                            i1 = i - 1
                            j1 = j
                            k1 = k
                        elif j <   0:
                            m = iEx[nj1]
                            d = Yn[1] - Yn[0]
                            i1 = i
                            j1 = j + 1
                            k1 = k
                        elif j >= Ny:
                            m = iEx[n]
                            d = Yn[-1] - Yn[-2]
                            i1 = i
                            j1 = j - 1
                            k1 = k
    
                        iMurHz[num, 3] = i1
                        iMurHz[num, 4] = j1
                        iMurHz[num, 5] = k1
                        fMurHz[num, 0] = 0
                        fMurHz[num, 1] = sol.MurH.factor(fMaterial, d, m, cdt)
                    num += 1

    if mode == 0:
        return num
