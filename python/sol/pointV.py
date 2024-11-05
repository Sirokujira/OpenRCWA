# -*- coding: utf-8 -*-
"""
pointV.py
"""

import sol.planewave

# 指定した時刻の観測点の電圧を代入する
def v(
    itime, Parm, fPlanewave, Xn, Yn, Zn, Xc, Yc, Zc, iPoint, fPoint, VPoint,
    Ex, Ey, Ez,
    iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0):

    npoint = iPoint.shape[0]
    #print(npoint, itime)

    if npoint == 0:
        return

    source = Parm['source']
    t = (itime + 1) * Parm['dt']

    for ipoint in range(npoint):
        idir = iPoint[ipoint, 0]
        i    = iPoint[ipoint, 1]
        j    = iPoint[ipoint, 2]
        k    = iPoint[ipoint, 3]
        ds   = fPoint[ipoint, 3 + idir]  # = dx/dy/dz (X/Y/Z)

        n = (Ni * i) + (Nj * j) + (Nk * k) + N0

        e = 0
        if   (idir == 0) and \
            (iMin <= i) and (i <  iMax) and \
            (jMin <= j) and (j <= jMax) and \
            (kMin <= k) and (k <= kMax):  # MPI
            e = Ex[n]
            if source == 1:
                fi, _ = sol.planewave.f(Xc[i], Yn[j], Zn[k], t, fPlanewave[0], fPlanewave)
                e += fi
        elif (idir == 1) and \
            (iMin <= i) and (i <= iMax) and \
            (jMin <= j) and (j <  jMax) and \
            (kMin <= k) and (k <= kMax):  # MPI
            e = Ey[n]
            if source == 1:
                fi, _ = sol.planewave.f(Xn[i], Yc[j], Zn[k], t, fPlanewave[1], fPlanewave)
                e += fi
        elif (idir == 2) and \
             (iMin <= i) and (i <= iMax) and \
             (jMin <= j) and (j <= jMax) and \
             (kMin <= k) and (k <  kMax):  # MPI
            e = Ez[n]
            if source == 1:
                fi, _ = sol.planewave.f(Xn[i], Yn[j], Zc[k], t, fPlanewave[2], fPlanewave)
                e += fi

        VPoint[ipoint, itime] = e * (-ds)
