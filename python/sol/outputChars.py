# -*- coding: utf-8 -*-
"""
outputChars.py
"""

import sys, math, cmath
import numpy as np
import sol.feed, sol.point, sol.farfield

# 計算結果の一部を出力する(計算の確認)
def out(
    fp, Parm, Nx, Ny, Nz, Xn, Yn, Zn, Xc, Yc, Zc,
    iFeed, fFeed, iPoint, Freq1, Freq2,
    cEx, cEy, cEz, cHx, cHy, cHz,
    VFeed, IFeed, VPoint, Ntime, cFdft,
    Ni, Nj, Nk, N0, NN):

    # 空配列
    Zin      = np.zeros((0, 0),    'c16')
    Ref      = np.zeros((0, 0),    'f8')
    Pin      = np.zeros((0, 0, 0), 'f8')
    Spara    = np.zeros((0, 0),    'c16')
    Coupling = np.zeros((0, 0, 0), 'c16')

    # 遠方界計算の準備(ポスト処理でも使う)
    fSurface, cSurface \
    = sol.farfield.setup(
        Nx, Ny, Nz, Xn, Yn, Zn, Xc, Yc, Zc, Freq2, Parm['pbc'],
        cEx, cEy, cEz, cHx, cHy, cHz,
        Ni, Nj, Nk, N0, NN)

    # 入力インピーダンス、反射係数
    if iFeed.shape[0] > 0:
        # 計算
        Zin, Ref = sol.feed.calcZin(Parm, fFeed, Freq1, VFeed, IFeed, Ntime)
        Pin = sol.feed.calcPin(Parm, fFeed, Freq2, VFeed, IFeed, Ntime, cFdft)  # ポスト処理用
        # 出力
        _outputZin(fp,         iFeed, fFeed, Freq1, Zin, Ref)
        _outputZin(sys.stdout, iFeed, fFeed, Freq1, Zin, Ref)

    # Sパラメーター
    if iPoint.shape[0] > 0:
        # 計算
        Spara = sol.point.calcSpara(Parm, iPoint, Freq1, VPoint, Ntime)
        # 出力
        _outputSpara(fp,         iPoint, Freq1, Spara)
        _outputSpara(sys.stdout, iPoint, Freq1, Spara)

    # 結合度
    if (iFeed.shape[0] > 0) and (iPoint.shape[0] > 2):
        # 計算
        Coupling = sol.point.calcCoupling(Parm, iFeed, iPoint, Freq1, VFeed, VPoint, Ntime)
        # 出力
        _outputCoupling(fp,         iFeed, iPoint, Freq1, Coupling)
        _outputCoupling(sys.stdout, iFeed, iPoint, Freq1, Coupling)

    # 散乱断面積
    if Parm['source'] == 1:
        # 計算
        bcs, fcs = _cs(Parm, iFeed, Freq2, Pin, fSurface, cSurface)
        # 出力
        _outputCross(fp,         Freq2, bcs, fcs)
        _outputCross(sys.stdout, Freq2, bcs, fcs)

    return Zin, Ref, Pin, Spara, Coupling, fSurface, cSurface

# (private) 散乱断面積を計算する(第2周波数)
def _cs(Parm, iFeed, Freq2, Pin, fSurface, cSurface):

    nfeed = iFeed.shape[0]
    nfreq = len(Freq2)

    # alloc
    bcs = np.zeros(nfreq, 'f8')
    fcs = np.zeros(nfreq, 'f8')

    theta = Parm['planewave'][0]
    phi   = Parm['planewave'][1]
    C     = Parm['C']
    ETA0  = Parm['ETA0']

    for ifreq, freq in enumerate(Freq2):
        # 波数と遠方界係数
        k     = (2 * math.pi * freq) / C
        ffctr = sol.farfield.factor(ifreq, k, Pin, nfeed, 0, ETA0)

        # BCS
        etheta, ephi = sol.farfield.calc(ifreq,       theta,       phi, \
            ffctr, k, fSurface, cSurface)
        bcs[ifreq] = abs(etheta)**2 + abs(ephi)**2

        # FCS
        etheta, ephi = sol.farfield.calc(ifreq, 180 - theta, 180 + phi, \
            ffctr, k, fSurface, cSurface)
        fcs[ifreq] = abs(etheta)**2 + abs(ephi)**2

    return bcs, fcs

# (private) 入力インピーダンスと反射係数を出力する
def _outputZin(fp, iFeed, fFeed, Freq1, Zin, Ref):

    nfeed = iFeed.shape[0]

    fp.write("=== input impedance ===\n")

    for ifeed in range(nfeed):
        fp.write("feed #%d (Z0[ohm] = %.2f)\n" % (ifeed + 1, fFeed[ifeed, 5]))
        fp.write("  frequency[Hz] Rin[ohm]   Xin[ohm]    Gin[mS]    Bin[mS]    Ref[dB]       VSWR\n")
        for ifreq, freq in enumerate(Freq1):
            zin = Zin[ifeed, ifreq]
            yin = 1e3 / zin
            ref = Ref[ifeed, ifreq]
            refdb = 20 * math.log10(max(abs(ref), 1e-10))
            vswr = (1 + ref) / (1 - ref) if (abs(1 - ref) > 1e-6) else 1000
            fp.write("%13.5e%11.3f%11.3f%11.3f%11.3f%11.3f%11.3f\n" % \
                (freq, zin.real, zin.imag, yin.real, yin.imag, refdb, vswr))

    fp.flush()

# (private) Sパラメーターを出力する
def _outputSpara(fp, iPoint, Freq1, Spara):

    npoint = iPoint.shape[0] - 2

    fp.write("=== S-parameter ===\n")

    fp.write("  frequency[Hz]")
    for ipoint in range(npoint):
        fp.write("  S%d1[dB] S%d1[deg]" % (ipoint + 1, ipoint + 1))
    fp.write("\n")

    for ifreq in range(len(Freq1)):
        fp.write("%13.5e  " % Freq1[ifreq])
        for ipoint in range(npoint):
            s = Spara[ipoint, ifreq]
            amp = abs(s)
            ampdb = 20 * math.log10(max(amp, 1e-10))
            deg = math.degrees(cmath.phase(s))
            fp.write("%9.3f%9.3f" % (ampdb, deg))
        fp.write("\n")

    fp.flush()

# (private) 結合度を出力する
def _outputCoupling(fp, iFeed, iPoint, Freq1, Coupling):

    nfeed = iFeed.shape[0]
    npoint = iPoint.shape[0] - 2

    fp.write("=== coupling ===\n")

    fp.write("  frequency[Hz]")
    for ifeed in range(nfeed):
        for ipoint in range(npoint):
            fp.write("  C%d%d[dB] C%d%d[deg]" % (ipoint + 1, ifeed + 1, ipoint + 1, ifeed + 1))
    fp.write("\n")

    for ifreq, freq in enumerate(Freq1):
        fp.write("%13.5e  " % freq)
        for ifeed in range(nfeed):
            for ipoint in range(npoint):
                c = Coupling[ifeed, ipoint, ifreq]
                amp = abs(c)
                ampdb = 20 * math.log10(max(amp, 1e-10))
                deg = math.degrees(cmath.phase(c))
                fp.write("%9.3f%9.3f" % (ampdb, deg))
        fp.write("\n")

# (private) 散乱断面積を出力する
def _outputCross(fp, Freq2, bcs, fcs):

    fp.write("=== cross section ===\n")
    fp.write("  frequency[Hz] backward[m*m] forward[m*m]\n")

    for ifreq in range(len(Freq2)):
        fp.write("%13.5e  %12.4e  %12.4e\n" % (Freq2[ifreq], bcs[ifreq], fcs[ifreq]))

    fp.flush()
