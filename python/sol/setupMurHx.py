# -*- coding: utf-8 -*-
"""
setupMurHx.py
"""

from numba import jit
import sol.MurH

@jit(cache=True, nopython=True)
def setHx(mode,
    Ny, Nz, Yn, Zn, fMaterial, iEy, iEz, fMurHx, iMurHx, cdt,
    iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0):

    imin = iMin - 0
    jmin = jMin - 1
    kmin = kMin - 1
    imax = iMax + 1
    jmax = jMax + 1
    kmax = kMax + 1

    num = 0
    for i in range(imin, imax):
        for j in range(jmin, jmax):
            for k in range(kmin, kmax):
                if ((((j < 0) or (j >= Ny)) and (k >= 0) and (k < Nz)) or \
                    (((k < 0) or (k >= Nz)) and (j >= 0) and (j < Ny))):
                    if mode == 1:
                        iMurHx[num, 0] = i
                        iMurHx[num, 1] = j
                        iMurHx[num, 2] = k
                        m = d = i1 = j1 = k1 = 0
                        n   = Ni * (i    ) + Nj * (j    ) + Nk * (k    ) + N0
                        nj1 = Ni * (i    ) + Nj * (j + 1) + Nk * (k    ) + N0
                        nk1 = Ni * (i    ) + Nj * (j    ) + Nk * (k + 1) + N0
                        if   j <   0:
                            m = iEz[nj1]
                            d = Yn[1] - Yn[0]
                            i1 = i
                            j1 = j + 1
                            k1 = k
                        elif j >= Ny:
                            m = iEz[n]
                            d = Yn[-1] - Yn[-2]
                            i1 = i
                            j1 = j - 1
                            k1 = k
                        elif k <   0:
                            m = iEy[nk1]
                            d = Zn[1] - Zn[0]
                            i1 = i
                            j1 = j
                            k1 = k + 1
                        elif k >= Nz:
                            m = iEy[n]
                            d = Zn[-1] - Zn[-2]
                            i1 = i
                            j1 = j
                            k1 = k - 1
    
                        iMurHx[num, 3] = i1
                        iMurHx[num, 4] = j1
                        iMurHx[num, 5] = k1
                        fMurHx[num, 0] = 0
                        fMurHx[num, 1] = sol.MurH.factor(fMaterial, d, m, cdt)
                    num += 1

    if mode == 0:
        return num
