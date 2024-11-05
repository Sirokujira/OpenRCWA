# -*- coding: utf-8 -*-
"""
farfield.py
"""

import math, cmath
import numpy as np
from numba import jit, prange
import sol.Fnode

# 遠方界計算の準備
# PBCのときは面積=0として飛ばす
# 並列化は時間が増えるので不可
#@jit(cache=True, nogil=True, parallel=True, nopython=True)
@jit(cache=True, nopython=True)
def setup(
    Nx, Ny, Nz, Xn, Yn, Zn, Xc, Yc, Zc, Freq2, pbc,
    cEx, cEy, cEz, cHx, cHy, cHz,
    Ni, Nj, Nk, N0, NN):

    nfreq = len(Freq2)
    num = 2 * ((Nx * Ny) + (Ny * Nz) + (Nz * Nx))

    fSurface = np.zeros((7, num), 'f8')
    cSurface = np.zeros((6, num, nfreq), 'c8')

    # 作業配列
    cex = np.zeros((2, 2), 'c8')
    cey = np.zeros((2, 2), 'c8')
    cez = np.zeros((2, 2), 'c8')
    chx = np.zeros((2, 2), 'c8')
    chy = np.zeros((2, 2), 'c8')
    chz = np.zeros((2, 2), 'c8')

    # counter
    n = 0

    # X surface
    for side in range(2):
        i = [0, Nx][side]
        for j in range(Ny):
            for k in range(Nz):
                for ifreq in range(len(Freq2)):
                    for jm in range(2):
                        for km in range(2):
                            _, cey[jm][km], cez[jm][km] = \
                                sol.Fnode.e(ifreq, i, j + jm, k + km, \
                                    Nx, Ny, Nz, cEx, cEy, cEz, Ni, Nj, Nk, N0, NN)#NEx, NEy, NEz)
                            _, chy[jm][km], chz[jm][km] = \
                                sol.Fnode.h(ifreq, i, j + jm, k + km, \
                                    cHx, cHy, cHz, Ni, Nj, Nk, N0, NN)#NHx, NHy, NHz)
                    cSurface[0, n, ifreq] = complex(0, 0)
                    cSurface[1, n, ifreq] = np.sum(cey) / 4
                    cSurface[2, n, ifreq] = np.sum(cez) / 4
                    cSurface[3, n, ifreq] = complex(0, 0)
                    cSurface[4, n, ifreq] = np.sum(chy) / 4
                    cSurface[5, n, ifreq] = np.sum(chz) / 4
                # 周波数によらない形状データ
                fSurface[0, n] = [-1, +1][side]
                fSurface[1, n] = 0
                fSurface[2, n] = 0
                fSurface[3, n] = Xn[i]
                fSurface[4, n] = Yc[j]
                fSurface[5, n] = Zc[k]
                fSurface[6, n] = (Yn[j + 1] - Yn[j]) * (Zn[k + 1] - Zn[k]) * (1 - pbc[0])
                n += 1

    # Y surface
    for side in range(2):
        j = [0, Ny][side]
        for k in range(Nz):
            for i in range(Nx):
                for ifreq in range(len(Freq2)):
                    for km in range(2):
                        for im in range(2):
                            cex[km][im], _, cez[km][im] = \
                                sol.Fnode.e(ifreq, i + im, j, k + km, \
                                    Nx, Ny, Nz, cEx, cEy, cEz, Ni, Nj, Nk, N0, NN)#NEx, NEy, NEz)
                            chx[km][im], _, chz[km][im] = \
                                sol.Fnode.h(ifreq, i + im, j, k + km, \
                                    cHx, cHy, cHz, Ni, Nj, Nk, N0, NN)#NHx, NHy, NHz)
                    cSurface[0, n, ifreq] = np.sum(cex) / 4
                    cSurface[1, n, ifreq] = complex(0, 0)
                    cSurface[2, n, ifreq] = np.sum(cez) / 4
                    cSurface[3, n, ifreq] = np.sum(chx) / 4
                    cSurface[4, n, ifreq] = complex(0, 0)
                    cSurface[5, n, ifreq] = np.sum(chz) / 4
                # 周波数によらない形状データ
                fSurface[0, n] = 0
                fSurface[1, n] = [-1, +1][side]
                fSurface[2, n] = 0
                fSurface[3, n] = Xc[i]
                fSurface[4, n] = Yn[j]
                fSurface[5, n] = Zc[k]
                fSurface[6, n] = (Zn[k + 1] - Zn[k]) * (Xn[i + 1] - Xn[i]) * (1 - pbc[1])
                n += 1

    # Z surface
    for side in range(2):
        k = [0, Nz][side]
        for i in range(Nx):
            for j in range(Ny):
                for ifreq in range(len(Freq2)):
                    for im in range(2):
                        for jm in range(2):
                            cex[im][jm], cey[im][jm], _ = \
                                sol.Fnode.e(ifreq, i + im, j + jm, k, \
                                    Nx, Ny, Nz, cEx, cEy, cEz, Ni, Nj, Nk, N0, NN)#NEx, NEy, NEz)
                            chx[im][jm], chy[im][jm], _ = \
                                sol.Fnode.h(ifreq, i + im, j + jm, k, \
                                    cHx, cHy, cHz, Ni, Nj, Nk, N0, NN)#NHx, NHy, NHz)
                    cSurface[0, n, ifreq] = np.sum(cex) / 4
                    cSurface[1, n, ifreq] = np.sum(cey) / 4
                    cSurface[2, n, ifreq] = complex(0, 0)
                    cSurface[3, n, ifreq] = np.sum(chx) / 4
                    cSurface[4, n, ifreq] = np.sum(chy) / 4
                    cSurface[5, n, ifreq] = complex(0, 0)
                # 周波数によらない形状データ
                fSurface[0, n] = 0
                fSurface[1, n] = 0
                fSurface[2, n] = [-1, +1][side]
                fSurface[3, n] = Xc[i]
                fSurface[4, n] = Yc[j]
                fSurface[5, n] = Zn[k]
                fSurface[6, n] = (Xn[i + 1] - Xn[i]) * (Yn[j + 1] - Yn[j]) * (1 - pbc[2])
                n += 1

    assert(n == num)
    #print(np.sum(fSurface))
    #print(np.sum(cSurface))

    return fSurface, cSurface

@jit(cache=True, nogil=True, parallel=True, nopython=True)
#@jit(cache=True, nopython=True)
def calc(ifreq, theta, phi, ffctr, k, fSurface, cSurface):

    # sin, cos
    cost = math.cos(math.radians(theta))
    sint = math.sin(math.radians(theta))
    cosp = math.cos(math.radians(phi))
    sinp = math.sin(math.radians(phi))

    # (r1, t1, p1) : (r, theta, phi)方向の単位ベクトル
    r1 = np.array([+ sint * cosp, + sint * sinp, + cost])
    t1 = np.array([+ cost * cosp, + cost * sinp, - sint])
    p1 = np.array([- sinp,        + cosp,          0   ])

    # ポテンシャルL/N初期化(複素数ベクトル)
    pl = np.zeros(3, 'c8')
    pn = np.zeros(3, 'c8')

    # 境界要素に関するループ
    for n in prange(fSurface.shape[1]):
        # 境界要素
        nv = fSurface[0:3, n]
        pv = fSurface[3:6, n]
        ds = fSurface[6,   n]

        # 電磁界ベクトル
        ev = cSurface[0:3, n, ifreq]
        hv = cSurface[3:6, n, ifreq]

        # Z0 * J = n X (Z0 * H)
        cj = np.cross(nv, hv)

        # M = -n X E
        cm = - np.cross(nv, ev)

        # exp(jkr * r) * dS
        #expds = ds * cmath.exp(complex(0, k * np.dot(r1, pv)))
        rr = (r1[0] * pv[0]) + (r1[1] * pv[1]) + (r1[2] * pv[2])
        expds = ds * cmath.exp(complex(0, k * rr))

        # L += M * exp(jkr * r) * dS
        pl += cm * expds

        # Z0 * N += (Z0 * J) * exp(jkr * r) * dS
        pn += cj * expds

    # Z0 * N-theta, Z0 * N-phi (np.dotはjit不可)
    #pnt = np.dot(t1, pn)
    #pnp = np.dot(p1, pn)
    pnt = (t1[0] * pn[0]) + (t1[1] * pn[1]) + (t1[2] * pn[2])
    pnp = (p1[0] * pn[0]) + (p1[1] * pn[1]) + (p1[2] * pn[2])

    # L-theta, L-phi
    #plt = np.dot(t1, pl)
    #plp = np.dot(p1, pl)
    plt = (t1[0] * pl[0]) + (t1[1] * pl[1]) + (t1[2] * pl[2])
    plp = (p1[0] * pl[0]) + (p1[1] * pl[1]) + (p1[2] * pl[2])

    # F-theta, F-phi
    etheta = ffctr * (pnt + plp)
    ephi   = ffctr * (pnp - plt)

    return etheta, ephi

# 遠方界因子
def factor(ifreq, k, Pin, nfeed, matchingloss, ETA0):

    ffctr = 0
    if nfeed > 0:
        # 給電点 (post only)
        s = 0
        for ifeed in range(nfeed):
            s += 0.5 * Pin[ifreq, ifeed, matchingloss]
        ffctr = k / math.sqrt(8 * math.pi * ETA0 * s)
    else:
        # 平面波入射 (solver + post)
        einc = 1
        ffctr = k / (einc * math.sqrt(4 * math.pi))

    return ffctr

# 遠方界成分
# etheta, ephi : 複素数
def farComponent(etheta, ephi):

    e = np.zeros(7, 'f8')
    
    # abs
    e[0] = math.sqrt(abs(etheta)**2 + abs(ephi)**2)

    # theta/phi
    e[1] = abs(etheta)
    e[2] = abs(ephi)

    # major/minor
    tmp1 = abs(etheta**2) + abs(ephi**2)
    tmp2 = abs(etheta**2 + ephi**2)
    e[3] = math.sqrt((tmp1 + tmp2) / 2)
    e[4] = math.sqrt(max(tmp1 - tmp2, 0) / 2)

    # RHCP/LHCP
    e[5] = abs(etheta + 1j * ephi) / math.sqrt(2)
    e[6] = abs(etheta - 1j * ephi) / math.sqrt(2)

    return e
