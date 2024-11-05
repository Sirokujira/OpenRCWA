# -*- coding: utf-8 -*-
"""
inductor.py
"""

# inductor関係の変数を更新する
def calcL(
    Parm, iInductor, fInductor,
    Ex, Ey, Ez, Hx, Hy, Hz,
    iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0):

    ninductor = iInductor.shape[0]

    if ninductor <= 0:
        return

    cdt = Parm['C'] * Parm['dt']

    for m in range(ninductor):

        idir       = iInductor[m, 0]
        i, j, k    = iInductor[m, 1:4]
        dx, dy, dz = fInductor[m, 3:6]
        ldat       = fInductor[m, 7:10]  # factor, e-new, e-sum (作業配列のpointer)

        n   = Ni * (i    ) + Nj * (j    ) + Nk * (k    ) + N0
        ni1 = Ni * (i - 1) + Nj * (j    ) + Nk * (k    ) + N0
        nj1 = Ni * (i    ) + Nj * (j - 1) + Nk * (k    ) + N0
        nk1 = Ni * (i    ) + Nj * (j    ) + Nk * (k - 1) + N0

        if   (idir == 0) and \
            (iMin <= i) and (i <  iMax) and \
            (jMin <= j) and (j <= jMax) and \
            (kMin <= k) and (k <= kMax):  # MPI
            # X方向
            roth = (Hz[n] - Hz[nj1]) / dy \
                 - (Hy[n] - Hy[nk1]) / dz
            Ex[n] = ldat[1] + (cdt * roth) - (ldat[0] * cdt**2 * ldat[2])
            ldat[1] = Ex[n]
            ldat[2] += ldat[1]
        elif (idir == 1) and \
            (iMin <= i) and (i <= iMax) and \
            (jMin <= j) and (j <  jMax) and \
            (kMin <= k) and (k <= kMax):  # MPI
            # Y方向
            roth = (Hx[n] - Hx[nk1]) / dz \
                 - (Hz[n] - Hz[ni1]) / dx
            Ey[n] = ldat[1] + (cdt * roth) - (ldat[0] * cdt**2 * ldat[2])
            ldat[1] = Ey[n]
            ldat[2] += ldat[1]
        elif (idir == 2) and \
            (iMin <= i) and (i <= iMax) and \
            (jMin <= j) and (j <= jMax) and \
            (kMin <= k) and (k <  kMax):  # MPI
            # Z方向
            roth = (Hy[n] - Hy[ni1]) / dx \
                 - (Hx[n] - Hx[nj1]) / dy
            Ez[n] = ldat[1] + (cdt * roth) - (ldat[0] * cdt**2 * ldat[2])
            ldat[1] = Ez[n]
            ldat[2] += ldat[1]
