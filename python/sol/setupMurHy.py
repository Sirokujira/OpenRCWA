# -*- coding: utf-8 -*-
"""
setupMurHy.py
"""

from numba import jit
import sol.MurH

@jit(cache=True, nopython=True)
def setHy(mode,
    Nx, Nz, Xn, Zn, fMaterial, iEz, iEx, fMurHy, iMurHy, cdt,
    iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0):

    imin = iMin - 1
    jmin = jMin - 0
    kmin = kMin - 1
    imax = iMax + 1
    jmax = jMax + 1
    kmax = kMax + 1

    num = 0
    for i in range(imin, imax):
        for j in range(jmin, jmax):
            for k in range(kmin, kmax):
                if ((((k < 0) or (k >= Nz)) and (i >= 0) and (i < Nx)) or \
                    (((i < 0) or (i >= Nx)) and (k >= 0) and (k < Nz))):
                    if mode == 1:
                        iMurHy[num, 0] = i
                        iMurHy[num, 1] = j
                        iMurHy[num, 2] = k
                        m = d = i1 = j1 = k1 = 0
                        n   = Ni * (i    ) + Nj * (j    ) + Nk * (k    ) + N0
                        nk1 = Ni * (i    ) + Nj * (j    ) + Nk * (k + 1) + N0
                        ni1 = Ni * (i + 1) + Nj * (j    ) + Nk * (k    ) + N0
                        if   k <   0:
                            m = iEx[nk1]
                            d = Zn[1] - Zn[0]
                            i1 = i
                            j1 = j
                            k1 = k + 1
                        elif k >= Nz:
                            m = iEx[n]
                            d = Zn[-1] - Zn[-2]
                            i1 = i
                            j1 = j
                            k1 = k - 1
                        elif i <   0:
                            m = iEz[ni1]
                            d = Xn[1] - Xn[0]
                            i1 = i + 1
                            j1 = j
                            k1 = k
                        elif i >= Nx:
                            m = iEz[n]
                            d = Xn[-1] - Xn[-2]
                            i1 = i - 1
                            j1 = j
                            k1 = k
    
                        iMurHy[num, 3] = i1
                        iMurHy[num, 4] = j1
                        iMurHy[num, 5] = k1
                        fMurHy[num, 0] = 0
                        fMurHy[num, 1] = sol.MurH.factor(fMaterial, d, m, cdt)
                    num += 1

    if mode == 0:
        return num
