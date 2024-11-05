# -*- coding: utf-8 -*-
"""
dft.py
DFT関係
"""

import math
import numpy as np
import sol.feed, sol.planewave

# DFT用の係数を準備する
def setup(Parm, fPlanewave, Xn, Yn, Zn, Freq2):

    maxiter = Parm['solver'][0]
    dt      = Parm['dt']

    # 正規化因子
    dsize = len(Freq2)
    cFdft = np.zeros(dsize, 'c16')

    # 定数
    if   Parm['source'] == 0:
        # 給電点
        tw = Parm['tw']
    elif Parm['source'] == 1:
        # 平面波入射
        x0 = (Xn[0] + Xn[-1]) / 2
        y0 = (Yn[0] + Yn[-1]) / 2
        z0 = (Zn[0] + Zn[-1]) / 2

    # 入射界のDFT
    for ifreq in range(len(Freq2)):
        omega = 2 * math.pi * Freq2[ifreq]
        csum = 0
        for itime in range(maxiter + 1):
            t = itime * dt
            fi = 0
            if   Parm['source'] == 0:
                # 給電点
                fi = sol.feed.volt(t, tw, 0)
            elif Parm['source'] == 1:
                # 平面波入射
                fi, _ = sol.planewave.f(x0, y0, z0, t, 1.0, fPlanewave)
            phase = omega * t
            csum += fi * (math.cos(phase) - 1j * math.sin(phase))
        cFdft[ifreq] = csum

    # DFT係数
    dsize = len(Freq2) * (maxiter + 1)
    cEdft = np.zeros(dsize, 'c16')
    cHdft = np.zeros(dsize, 'c16')
    for itime in range(maxiter + 1):
        for ifreq in range(len(Freq2)):
            omega = 2 * math.pi * Freq2[ifreq]
            hphase = omega * (itime + 0.5) * dt
            ephase = omega * (itime + 1.0) * dt
            idx = (itime * len(Freq2)) + ifreq
            cHdft[idx] = (math.cos(hphase) - 1j * math.sin(hphase)) / cFdft[ifreq]
            cEdft[idx] = (math.cos(ephase) - 1j * math.sin(ephase)) / cFdft[ifreq]
    #print(cEdft)
    #print(cHdft)

    return cEdft, cHdft, cFdft

# 関数fのDFTを計算する (汎用関数)
def calc(ntime, f, freq, dt, shift):
    """
    omega = 2 * math.pi * freq
    csum = complex(0, 0)
    for n in range(ntime):
        ot = omega * (n + shift) * dt
        csum += cmath.exp(complex(0, - ot)) * f[n]

    return csum
    """
    omega = 2 * math.pi * freq
    n = np.arange(ntime)
    ot = omega * (n + shift) * dt
    cexp = np.cos(ot) - 1j * np.sin(ot)
    return np.sum(cexp * f[0:ntime])