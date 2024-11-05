# -*- coding: utf-8 -*-
"""
feed.py
給電点
"""

import math
import numpy as np
import sol.dft

# 指定した時刻の電圧
def volt(t, tw, td):
    arg = (t - tw - td) / (tw / 4.0)

    return math.sqrt(2.0) * math.exp(0.5) * arg * math.exp(-arg**2)

# V/I波形を保存する
def evi(
    itime, Parm, iFeed, fFeed, VFeed, IFeed,
    Ex, Ey, Ez, Hx, Hy, Hz, iEx, iEy, iEz,
    iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0):

    if iFeed.shape[0] <= 0:
        return

    eps   = Parm['EPS'] #1e-6
    dt    = Parm['dt']
    tw    = Parm['tw']
    rfeed = Parm['rfeed']
    PEC   = Parm['PEC']
    ETA0  = Parm['ETA0']

    t = (itime + 1) * dt

    for ifeed in range(iFeed.shape[0]):

        idir  = iFeed[ifeed, 0]
        i     = iFeed[ifeed, 1]
        j     = iFeed[ifeed, 2]
        k     = iFeed[ifeed, 3]

        v     = fFeed[ifeed, 3]
        delay = fFeed[ifeed, 4]
        dx    = fFeed[ifeed, 6]
        dy    = fFeed[ifeed, 7]
        dz    = fFeed[ifeed, 8]

        n   = (Ni * i) + (Nj * j) + (Nk * k) + N0
        ni1 = n - Ni
        nj1 = n - Nj
        nk1 = n - Nk

        # V
        v0 = volt(t, tw, delay)
        v *= v0

        # E, V, I
        c = 0
        if   (idir == 0) and \
            (iMin <= i) and (i <  iMax) and \
            (jMin <= j) and (j <= jMax) and \
            (kMin <= k) and (k <= kMax):  # MPI
            # X方向
            c = dz * (Hz[n] - Hz[nj1]) \
              - dy * (Hy[n] - Hy[nk1])
            c /= ETA0
            v -= rfeed * c
            if (iEx[n] == PEC) or (abs(v0) > eps):
                Ex[n] = -v / dx
        elif (idir == 1) and \
            (iMin <= i) and (i <= iMax) and \
            (jMin <= j) and (j <  jMax) and \
            (kMin <= k) and (k <= kMax):  # MPI
            # Y方向
            c = dx * (Hx[n] - Hx[nk1]) \
              - dz * (Hz[n] - Hz[ni1])
            c /= ETA0
            v -= rfeed * c
            if (iEy[n] == PEC) or (abs(v0) > eps):
                Ey[n] = -v / dy
        elif (idir == 2) and \
            (iMin <= i) and (i <= iMax) and \
            (jMin <= j) and (j <= jMax) and \
            (kMin <= k) and (k <  kMax):  # MPI
            # Z方向
            c = dy * (Hy[n] - Hy[ni1]) \
              - dx * (Hx[n] - Hx[nj1])
            c /= ETA0
            v -= rfeed * c
            if (iEz[n] == PEC) or (abs(v0) > eps):
                Ez[n] = -v / dz

        # V/I波形
        VFeed[ifeed, itime] = v
        IFeed[ifeed, itime] = c

# 入力インピーダンスと反射係数を計算する(第1周波数)
def calcZin(Parm, fFeed, Freq1, VFeed, IFeed, Ntime):

    dt = Parm['dt']
    nfeed = fFeed.shape[0]
    nfreq = len(Freq1)

    Zin = np.zeros((nfeed, nfreq), 'c16')
    Ref = np.zeros((nfeed, nfreq), 'f8')

    for ifeed in range(nfeed):
        fv = VFeed[ifeed, :]
        fi = IFeed[ifeed, :]
        z0 = fFeed[ifeed, 5]
        for ifreq, freq in enumerate(Freq1):
            # 入力インピーダンス Zin = Vin / Iin
            vin = sol.dft.calc(Ntime, fv, freq, dt, 0)
            iin = sol.dft.calc(Ntime, fi, freq, dt, -0.5)
            Zin[ifeed, ifreq] = vin / iin

            # 反射係数 Ref = (Zin - Z0) / (Zin + Z0) (振幅のみ)
            Ref[ifeed, ifreq] = abs((Zin[ifeed, ifreq] - z0) \
                                  / (Zin[ifeed, ifreq] + z0))

    return Zin, Ref

# 給電電力を計算する(第2周波数、ポスト処理で使用される)
def calcPin(Parm, fFeed, Freq2, VFeed, IFeed, Ntime, cFdft):

    dt = Parm['dt']
    nfeed = fFeed.shape[0]
    nfreq = len(Freq2)

    Pin = np.zeros((nfreq, nfeed, 2), 'f8')

    for ifeed in range(nfeed):
        fv = VFeed[ifeed, :]
        fi = IFeed[ifeed, :]
        z0 = fFeed[ifeed, 5]
        for ifreq, freq in enumerate(Freq2):
            vin = sol.dft.calc(Ntime, fv, freq, dt, 0)    / cFdft[ifreq]
            iin = sol.dft.calc(Ntime, fi, freq, dt, -0.5) / cFdft[ifreq]
            zin = vin / iin
            ref = abs((zin - z0) / (zin + z0))
            Pin[ifreq, ifeed, 0] = (vin * iin.conjugate()).real
            Pin[ifreq, ifeed, 1] = Pin[ifreq, ifeed, 0] / max(1 - ref**2, Parm['EPS'])

    return Pin
