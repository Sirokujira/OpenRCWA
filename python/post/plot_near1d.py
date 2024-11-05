# -*- coding: utf-8 -*-
"""
plot_near1d.py
"""

import math
import numpy as np
import matplotlib.pyplot as plt
import post.common

# 近傍界線上分布図
def plot(Post, fPlanewave, Xn, Yn, Zn, Freq2, cEx, cEy, cEz, cHx, cHy, cHz, Ni, Nj, Nk, N0, NN):

    Ni = Post['Ni']
    Nj = Post['Nj']
    Nk = Post['Nk']
    N0 = Post['N0']
    NN = Post['NN']
    Nx = len(Xn) - 1
    Ny = len(Yn) - 1
    Nz = len(Zn) - 1

    nfreq  = len(Freq2)
    nline  = len(Post['n1ddir'])
    source = Post['source']
    noinc  = Post['n1dnoinc']
    #print(Ni, Nj, Nk, N0, NN, Nx, Ny, Nz)
    #print(nline, nfreq)

    if (nline < 1) or (nfreq < 1):
        return

    # log
    fname = 'near1d.log'
    fp = open(fname, 'wt', encoding='utf-8')

    # plot
    nfig = 0
    for iline in range(nline):
        compo = Post['n1dcompo'][iline]
        idir = Post['n1ddir'][iline]
        [pos1, pos2] = Post['n1dpos'][iline]

        # 単位
        if Post['n1ddb'] == 1:
            strunit = '[dBV/m]' if compo.startswith('E') else '[dBA/m]'
        else:
            strunit = '[V/m]' if compo.startswith('E') else '[A/m]'

        # 座標
        i = j = k = 0
        x = []  # 横軸
        xlabel = ''
        if   idir == 0:
            j = np.argmin(abs(Yn - pos1))
            k = np.argmin(abs(Zn - pos2))
            x = Xn
            xlabel = 'X [m]  (Y = %.4gm, Z = %.4gm)' % (Yn[j], Zn[k])
        elif idir == 1:
            k = np.argmin(abs(Zn - pos1))
            i = np.argmin(abs(Xn - pos2))
            x = Yn
            xlabel = 'Y [m]  (Z = %.4gm, X = %.4gm)' % (Zn[k], Xn[i])
        elif idir == 2:
            i = np.argmin(abs(Xn - pos1))
            j = np.argmin(abs(Yn - pos2))
            x = Zn
            xlabel = 'Z [m]  (X = %.4gm, Y = %.4gm)' % (Xn[i], Yn[j])
        #print(iline, i, j, k)

        # 縦軸(14成分)
        # E, Ex, Ey, Ez, deg(Ex), deg(Ey), deg(Ez)
        # H, Hx, Hy, Hz, deg(Hx), deg(Hy), deg(Hz)
        y = np.zeros((len(x), 14), 'f8')

        for ifreq, freq in enumerate(Freq2):
            # 波数k
            kwave = (2 * math.pi * freq) / Post['C']

            # E, H
            if   idir == 0:
                for i in range(Nx + 1):
                    y[i, :] = post.common.fnode( \
                        ifreq, i, j, k, \
                        kwave, source, fPlanewave, noinc, \
                        Nx, Ny, Nz, Xn, Yn, Zn, cEx, cEy, cEz, cHx, cHy, cHz, Ni, Nj, Nk, N0, NN)
            elif idir == 1:
                for j in range(Ny + 1):
                    y[j, :] = post.common.fnode( \
                        ifreq, i, j, k, \
                        kwave, source, fPlanewave, noinc, \
                        Nx, Ny, Nz, Xn, Yn, Zn, cEx, cEy, cEz, cHx, cHy, cHz, Ni, Nj, Nk, N0, NN)
            elif idir == 2:
                for k in range(Nz + 1):
                    y[k, :] = post.common.fnode( \
                        ifreq, i, j, k, \
                        kwave, source, fPlanewave, noinc, \
                        Nx, Ny, Nz, Xn, Yn, Zn, cEx, cEy, cEz, cHx, cHy, cHz, Ni, Nj, Nk, N0, NN)
            #print(y[:, 0])

            # log
            _log_n1d(fp, iline, freq, compo, idir, x, pos1, pos2, y)

            # dBに変換する
            if Post['n1ddb'] == 1:
                for m in [0, 1, 2, 3, 7, 8, 9, 10]:
                    y[:, m] = 20 * np.log10(np.maximum(y[:, m], 1e-10))
            
            # 描画する成分のリスト
            icmp = []
            scol = []
            scmp = []
            if   compo == 'E':
                icmp = [0, 1, 2, 3]
                scol = ['k', 'r', 'g', 'b']
                scmp = ['E', 'Ex', 'Ey', 'Ez']
            elif compo == 'Ex':
                icmp = [1, 4]
            elif compo == 'Ey':
                icmp = [2, 5]
            elif compo == 'Ez':
                icmp = [3, 6]
            elif compo == 'H':
                icmp = [7, 8, 9, 10]
                scol = ['k', 'r', 'g', 'b']
                scmp = ['H', 'Hx', 'Hy', 'Hz']
            elif compo == 'Hx':
                icmp = [8, 11]
            elif compo == 'Hy':
                icmp = [9, 12]
            elif compo == 'Hz':
                icmp = [10, 13]

            # 最大/最小
            dmax = np.max(y[:, icmp[0]], axis=0)
            if Post['n1dscale'][0] == 0:
                # 自動スケール
                ymax = dmax
                if Post['n1ddb'] == 1:
                    ymin = ymax - 50
                else:
                    ymin = 0
                ydiv = 10
            else:
                # 指定スケール
                ymin = Post['n1dscale'][1]
                ymax = Post['n1dscale'][2]
                ydiv = Post['n1dscale'][3]
            #print(ymin, ymax, ydiv)

            # figure
            nfig += 1
            strfig = 'OpenRCWA - near field 1d (%d/%d)' % (nfig, nline * nfreq)
            fig = plt.figure(strfig, figsize=(Post['w2d'][0], Post['w2d'][1]))
            ax = fig.add_subplot()

            # plot
            if (compo == 'E') or (compo == 'H'):
                # E/Ex/Ey/Ez or H/Hx/Hy/Hz : amplitude only
                for m in range(len(icmp)):
                    ax.plot(x, y[:, icmp[m]], color=scol[m], label=scmp[m])
                ax.set_xlim(x[0], x[-1])
                ax.set_ylabel(compo + ' ' + strunit)
                ax.set_ylim(ymin, ymax)
                if Post['n1dscale'][0] == 1:
                    ax.set_yticks(np.linspace(ymin, ymax, ydiv + 1))
                ax.legend(loc='best')
                ax.grid(True)
            else:
                # Ex/Ey/Ez/Hx/Hy/Hz : amplitude and phase
                # amplitude
                ax.plot(x, y[:, icmp[0]], color='k')
                ax.set_xlim(x[0], x[-1])
                ax.set_ylabel('amplitude ' + strunit)
                ax.set_ylim(ymin, ymax)
                if Post['n1dscale'][0] == 1:
                    ax.set_yticks(np.linspace(ymin, ymax, ydiv + 1))
                ax.grid(True)

                # phase
                ax2 = ax.twinx()
                plt.plot(x, y[:, icmp[1]], color='r', linestyle='--')
                ax2.set_xlim(x[0], x[-1])
                ax2.set_ylabel('phase [deg]', color='r')
                ax2.set_ylim(-180, 180)
                ax2.set_yticks([-180, -135, -90, -45, 0, 45, 90, 135, 180])
                ax2.grid(False)

            # x-label
            ax.set_xlabel(xlabel)
            
            # title
            ax.set_title('%s\n%s, f = %.3f%s, max = %.4g%s' %
                (Post['title'], compo, freq * Post['fscale'], Post['funit'], dmax, strunit))

# (private) near1d.log
def _log_n1d(fp, iline, freq, compo, idir, x, p1, p2, y):

    # header
    fp.write('#%d : frequency[Hz] = %.3e\n' % (iline, freq))
    fp.write(' No.     X[m]        Y[m]        Z[m]       ')
    if compo.startswith('E'):
        fp.write('E[V/m]      Ex[V/m]   Ex[deg]    Ey[V/m]   Ey[deg]    Ez[V/m]   Ez[deg]\n')
    else:
        fp.write('H[A/m]      Hx[A/m]   Hx[deg]    Hy[A/m]   Hy[deg]    Hz[A/m]   Hz[deg]\n')

    # body
    p = [0] * 3
    fmt1 = '%4d%12.3e%12.3e%12.3e'
    fmt2 = '%12.4e%12.4e%9.3f%12.4e%9.3f%12.4e%9.3f\n'
    for n0, p0 in enumerate(x):
        if   idir == 0:
            p = [p0, p1, p2]
        elif idir == 1:
            p = [p2, p0, p1]
        elif idir == 2:
            p = [p1, p2, p0]
        f = y[n0, :]
        fp.write(fmt1 % (n0, p[0], p[1], p[2]))
        if compo.startswith('E'):
            fp.write(fmt2 % (f[0], f[1], f[4], f[2], f[5], f[3], f[6]))
        else:
            fp.write(fmt2 % (f[7], f[8], f[11], f[9], f[12], f[10], f[13]))
