# -*- coding: utf-8 -*-
"""
feed.py (CUDA)
給電点
"""

import math
from numba import cuda

def evi(block1d,
    itime, Parm, iFeed, fFeed, VFeed, IFeed,
    Ex, Ey, Ez, Hx, Hy, Hz, iEx, iEy, iEz,
    iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0):

    nfeed = iFeed.shape[0]

    t = (itime + 1) * Parm['dt']

    block = min(block1d, nfeed)
    grid = math.ceil(nfeed / block)

    _evi_gpu[grid, block](
        itime, t,
        Parm['tw'], Parm['rfeed'], Parm['PEC'], Parm['ETA0'],
        nfeed, iFeed, fFeed, VFeed, IFeed,
        Ex, Ey, Ez, Hx, Hy, Hz, iEx, iEy, iEz,
        iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)

# (private)
@cuda.jit(cache=True)
def _evi_gpu(
    itime, t,
    tw, rfeed, PEC, ETA0,
    nfeed, iFeed, fFeed, VFeed, IFeed,
    Ex, Ey, Ez, Hx, Hy, Hz, iEx, iEy, iEz,
    iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0):

    eps = 1e-6

    ifeed = cuda.grid(1)

    if ifeed < nfeed:

        idir  = iFeed[ifeed, 0]
        i     = iFeed[ifeed, 1]
        j     = iFeed[ifeed, 2]
        k     = iFeed[ifeed, 3]

        v     = fFeed[ifeed, 3]
        delay = fFeed[ifeed, 4]
        dx    = fFeed[ifeed, 6]
        dy    = fFeed[ifeed, 7]
        dz    = fFeed[ifeed, 8]

        n = (Ni * i) + (Nj * j) + (Nk * k) + N0

        # 給電電圧
        v0 = _volt(t, tw, delay)
        v *= v0

        # 給電点のV/Iを求め、給電点にEを代入する
        c = 0
        if   (idir == 0) and \
            (iMin <= i) and (i <  iMax) and \
            (jMin <= j) and (j <= jMax) and \
            (kMin <= k) and (k <= kMax):  # MPI
            # X方向
            c = dz * (Hz[n] - Hz[n - Nj]) \
              - dy * (Hy[n] - Hy[n - Nk])
            c /= ETA0
            v -= rfeed * c
            if (iEx[n] == PEC) or (abs(v0) > eps):
                Ex[n] = -v / dx
        elif (idir == 1) and \
            (iMin <= i) and (i <= iMax) and \
            (jMin <= j) and (j <  jMax) and \
            (kMin <= k) and (k <= kMax):  # MPI
            # Y方向
            c = dx * (Hx[n] - Hx[n - Nk]) \
              - dz * (Hz[n] - Hz[n - Ni])
            c /= ETA0
            v -= rfeed * c
            if (iEy[n] == PEC) or (abs(v0) > eps):
                Ey[n] = -v / dy
        elif (idir == 2) and \
            (iMin <= i) and (i <= iMax) and \
            (jMin <= j) and (j <= jMax) and \
            (kMin <= k) and (k <  kMax):  # MPI
            # Z方向
            c = dy * (Hy[n] - Hy[n - Ni]) \
              - dx * (Hx[n] - Hx[n - Nj])
            c /= ETA0
            v -= rfeed * c
            if (iEz[n] == PEC) or (abs(v0) > eps):
                Ez[n] = -v / dz

        # 給電点のV/I時間波形を保存する
        #print(ifeed, itime, v, c)
        VFeed[ifeed, itime] = v
        IFeed[ifeed, itime] = c

# 指定した時刻の電圧 (device関数)
@cuda.jit(device=True, cache=True)
def _volt(t, tw, td):
    arg = (t - tw - td) / (tw / 4.0)

    return math.sqrt(2.0) * math.exp(0.5) * arg * math.exp(-arg**2)
