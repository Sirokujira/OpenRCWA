# -*- coding: utf-8 -*-
"""
point.py
観測点
"""

import cmath
import numpy as np
import sol.dft

# Sパラメーターの周波数特性を計算する
def calcSpara(Parm, iPoint, Freq1, VPoint, Ntime):

    npoint = iPoint.shape[0]
    nfreq = len(Freq1)

    if (npoint <= 2) or (nfreq <= 0):
        return

    dt  = Parm['dt']
    eps = 1e-10#Parm['EPS']
    
    # alloc
    Spara = np.zeros((npoint - 2, nfreq), 'c16')  # 1+/1-は除く
    cv    = np.zeros((npoint, nfreq), 'c16')

    # DFT
    for ipoint in range(npoint):
        for ifreq in range(nfreq):
            freq  = Freq1[ifreq]
            cv[ipoint, ifreq] = sol.dft.calc(Ntime, VPoint[ipoint, :], freq, dt, 0)

    # S-parameters
    for ifreq in range(nfreq):
        cv0 = cv[0,          ifreq]  # V1
        cvp = cv[npoint - 2, ifreq]  # V1+
        cvm = cv[npoint - 1, ifreq]  # V1-
        if (abs(cv0) < eps) or (abs(cvp) < eps) or (abs(cvm) < eps):
            continue  # 0割回避
        c1 = (cvp + cvm) / cv0
        c2 = cmath.sqrt(c1**2 - 4)
        c3 = c1 + c2
        if c3.imag < 0:
            c3 = c1 - c2    # Im > 0
        c3 = 0.5 * c3       # exp(+gd)
        c4 = 1.0 / c3       # exp(-gd)
        c5 = c4**2 - c3**2  # exp(-2gd) - exp(2gd)
        c6 = ((cvp * c4) - (cvm * c3)) / c5  # V+
        c7 = ((cvm * c4) - (cvp * c3)) / c5  # V-
        # S11 = E- / E+
        Spara[0, ifreq] = c7 / c6
        # Sn1 (n > 1)
        for ipoint in range(1, npoint - 2):
            Spara[ipoint, ifreq] = cv[ipoint, ifreq] / c6  # Sn1 = Vn / V+

    return Spara

# 結合度の周波数特性を計算する
def calcCoupling(Parm, iFeed, iPoint, Freq1, VFeed, VPoint, Ntime):

    npoint = iPoint.shape[0] - 2
    nfeed = iFeed.shape[0]
    nfreq = len(Freq1)

    if (npoint <= 0) or (nfeed <= 0) or (nfreq <= 0):
        return

    dt = Parm['dt']

    # DFT
    cvf = np.zeros((nfeed,  nfreq), 'c16')
    cvp = np.zeros((npoint, nfreq), 'c16')
    for ifeed in range(nfeed):
        for ifreq, freq in enumerate(Freq1):
            cvf[ifeed, ifreq] = sol.dft.calc(Ntime, VFeed[ifeed, :], freq, dt, 0)
    for ipoint in range(npoint):
        for ifreq, freq in enumerate(Freq1):
            cvp[ipoint, ifreq] = sol.dft.calc(Ntime, VPoint[ipoint, :], freq, dt, 0)

    # 結合度
    Coupling = np.zeros((nfeed, npoint, nfreq), 'c16')
    for ifeed in range(nfeed):
        for ipoint in range(npoint):
            for ifreq in range(nfreq):
                Coupling[ifeed, ipoint, ifreq] = cvp[ipoint, ifreq] / cvf[ifeed, ifreq]

    return Coupling
