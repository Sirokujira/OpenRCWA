# -*- coding: utf-8 -*-
"""
plot_near2d.py
"""

import math
import numpy as np
import matplotlib.pyplot as plt
import post.common

# 近傍界面上分布図 (2D only)
# TODO : 3D, animation
def plot(Post, fPlanewave, Xn, Yn, Zn, Freq2, cEx, cEy, cEz, cHx, cHy, cHz, Ni, Nj, Nk, N0, NN, gline):

    Ni = Post['Ni']
    Nj = Post['Nj']
    Nk = Post['Nk']
    N0 = Post['N0']
    NN = Post['NN']
    Nx = len(Xn) - 1
    Ny = len(Yn) - 1
    Nz = len(Zn) - 1

    nfreq  = len(Freq2)
    nplane = len(Post['n2ddir'])
    source = Post['source']
    noinc  = Post['n2dnoinc']
    zoom   = Post['n2dzoom']
    #print(nplane, nfreq)
    #print(gline)
    #print(gline.shape)
    #print(zoom)

    if (nplane < 1) or (nfreq < 1):
        return
    if (Post['n2dcontour'] < 0) or (Post['n2dcontour'] > 3):
        return

    # log
    fname = 'near2d.log'
    fp = open(fname, 'wt', encoding='utf-8')

    # ページ数
    nfig = 0
    for iplane in range(nplane):
        compo = Post['n2dcompo'][iplane]
        if (compo == 'E') or (compo == 'H'):
            nfig += 1 * nfreq  # 振幅:1ページ
        else:
            nfig += 2 * nfreq  # 振幅と位相:2ページ

    ifig = 0
    for iplane in range(nplane):
        # setup
        idir  = Post['n2ddir'][iplane]  # =0/1/2:X/Y/Z
        compo = Post['n2dcompo'][iplane]
        pos   = Post['n2dpos'][iplane]
        #print(iplane, idir, compo, pos)

        # 単位
        strunit = ['', '[deg]']
        if Post['n2ddb'] == 1:
            strunit[0] = '[dBV/m]' if compo.startswith('E') else '[dBA/m]'
        else:
            strunit[0] = '[V/m]' if compo.startswith('E') else '[A/m]'

        # 座標
        i = j = k = 0
        i0 = i1 = j0 = j1 = k0 = k1 = 0
        x = y = []  # 横軸、縦軸
        xlabel = ylabel = ''
        if   idir == 0:
            # X面
            i = np.argmin(abs(Xn - pos))
            j0 = 0  if zoom[0] == 0 else np.argmin(abs(Yn - zoom[1]))
            j1 = Ny if zoom[0] == 0 else np.argmin(abs(Yn - zoom[2]))
            k0 = 0  if zoom[0] == 0 else np.argmin(abs(Zn - zoom[3]))
            k1 = Nz if zoom[0] == 0 else np.argmin(abs(Zn - zoom[4]))
            x = Yn[j0: j1 + 1]
            y = Zn[k0: k1 + 1]
            clabel = 'X'
            xlabel = 'Y'
            ylabel = 'Z'
        elif idir == 1:
            # Y面
            j = np.argmin(abs(Yn - pos))
            k0 = 0  if zoom[0] == 0 else np.argmin(abs(Zn - zoom[1]))
            k1 = Nz if zoom[0] == 0 else np.argmin(abs(Zn - zoom[2]))
            i0 = 0  if zoom[0] == 0 else np.argmin(abs(Xn - zoom[3]))
            i1 = Nx if zoom[0] == 0 else np.argmin(abs(Xn - zoom[4]))
            x = Zn[k0: k1 + 1]
            y = Xn[i0: i1 + 1]
            clabel = 'Y'
            xlabel = 'Z'
            ylabel = 'X'
        elif idir == 2:
            # Z面
            k = np.argmin(abs(Zn - pos))
            i0 = 0  if zoom[0] == 0 else np.argmin(abs(Xn - zoom[1]))
            i1 = Nx if zoom[0] == 0 else np.argmin(abs(Xn - zoom[2]))
            j0 = 0  if zoom[0] == 0 else np.argmin(abs(Yn - zoom[3]))
            j1 = Ny if zoom[0] == 0 else np.argmin(abs(Yn - zoom[4]))
            x = Xn[i0: i1 + 1]
            y = Yn[j0: j1 + 1]
            clabel = 'Z'
            xlabel = 'X'
            ylabel = 'Y'
        #print(iplane, i, j, k)
        #print(len(x), len(y))
        #print(x)
        #print(y)
        if (len(x) < 2) or (len(y) < 2):
            print("*** no area, plane #%d" % (iplane + 1))
            continue

        # 縦軸(14成分)
        # E, Ex, Ey, Ez, deg(Ex), deg(Ey), deg(Ez)
        # H, Hx, Hy, Hz, deg(Hx), deg(Hy), deg(Hz)
        # (Y, X)の順に注意
        f = np.zeros((len(y), len(x), 14), 'f8')
 
        # 図示する成分番号(1個/2個)
        icmp = []
        if   compo == 'E':
            icmp = [0]
        elif compo == 'Ex':
            icmp = [1, 4]
        elif compo == 'Ey':
            icmp = [2, 5]
        elif compo == 'Ez':
            icmp = [3, 6]
        elif compo == 'H':
            icmp = [7]
        elif compo == 'Hx':
            icmp = [8, 11]
        elif compo == 'Hy':
            icmp = [9, 12]
        elif compo == 'Hz':
            icmp = [10, 13]
        
        for ifreq, freq in enumerate(Freq2):
            # 波数k
            kwave = (2 * math.pi * freq) / Post['C']

            # E, H
            if   idir == 0:
                for j in range(len(x)):
                    for k in range(len(y)):
                        f[k, j, :] = post.common.fnode( \
                            ifreq, i0 + i, j0 + j, k0 + k, \
                            kwave, source, fPlanewave, noinc, \
                            Nx, Ny, Nz, Xn, Yn, Zn, cEx, cEy, cEz, cHx, cHy, cHz, Ni, Nj, Nk, N0, NN)
            elif idir == 1:
                for k in range(len(x)):
                    for i in range(len(y)):
                        f[i, k, :] = post.common.fnode( \
                            ifreq, i0 + i, j0 + j, k0 + k, \
                            kwave, source, fPlanewave, noinc, \
                            Nx, Ny, Nz, Xn, Yn, Zn, cEx, cEy, cEz, cHx, cHy, cHz, Ni, Nj, Nk, N0, NN)
            elif idir == 2:
                for i in range(len(x)):
                    for j in range(len(y)):
                        f[j, i, :] = post.common.fnode( \
                            ifreq, i0 + i, j0 + j, k0 + k, \
                            kwave, source, fPlanewave, noinc, \
                            Nx, Ny, Nz, Xn, Yn, Zn, cEx, cEy, cEz, cHx, cHy, cHz, Ni, Nj, Nk, N0, NN)

            # log
            _log_n2d(fp, iplane, freq, compo, idir, pos, x, y, f)

            # dBに変換する
            if Post['n2ddb'] == 1:
                for m in [0, 1, 2, 3, 7, 8, 9, 10]:
                    f[:, :, m] = 20 * np.log10(np.maximum(f[:, :, m], 1e-10))

            # 最小/最大
            dmax = np.max(f[:, :, icmp[0]])
            if Post['n2dscale'][0] == 0:
                # 自動スケール
                fmax = dmax
                if Post['n2ddb'] == 1:
                    fmin = fmax - 30
                else:
                    fmin = 0
            else:
                # 指定スケール
                fmin = Post['n2dscale'][1]
                fmax = Post['n2dscale'][2]

            # 値が一定の時は描かない
            if abs(fmin - fmax) < 1e-10:
                print('*** constant data : %s = %e\n' % (compo, fmin))
                continue

            # 下限=最小値
            f = np.maximum(f, fmin)

            # plot, m = 0/1 : amplitude/phase
            for m in range(len(icmp)):
                # figure
                ifig += 1
                strfig = 'OpenRCWA - near field 2d (%d/%d)' % (ifig, nfig)
                fig = plt.figure(strfig, figsize=(Post['w2d'][0], Post['w2d'][1]))
                ax = fig.add_subplot()

                # 等高線のレベル
                if m == 0:
                    levels = np.linspace(fmin, fmax, Post['n2dscale'][3] + 1)
                else:
                    levels = np.linspace(-180, 180, Post['n2dscale'][3] + 1)
 
                # 等高線図
                if   (Post['n2dcontour'] == 0) or (Post['n2dcontour'] == 2):
                    # カラー塗りつぶし
                    CS = ax.contourf(x, y, f[:, :, icmp[m]], levels, cmap='rainbow')
                elif (Post['n2dcontour'] == 1) or (Post['n2dcontour'] == 3):
                    #　カラー等高線
                    CS = ax.contour(x, y, f[:, :, icmp[m]], levels, cmap='rainbow')
                ax.set_aspect('equal')

                # カラーバー
                if m == 0:
                    cbar = fig.colorbar(CS)
                else:
                    cbar = fig.colorbar(CS, ticks=[-180, -120, -60, 0, 60, 120, 180])
                cbar.set_label(compo + ' ' + strunit[m])

                # 物体形状を上書きする
                if Post['n2dobj'] == 1:
                    _n2dobj(ax, idir, gline, [x[0], x[-1]], [y[0], y[-1]])

                # X,Yラベル
                ax.set_xlabel(xlabel + ' [m]')
                ax.set_ylabel(ylabel + ' [m]')

                # タイトル
                strmax = ', max = %.4g%s' % (dmax, strunit[0]) if m == 0 else ''
                ax.set_title('%s\n%s = %g[m], f = %.3f%s%s' % \
                    (Post['title'], clabel, pos, freq * Post['fscale'], Post['funit'], strmax))
    
# (private) 物体形状を上書きする
def _n2dobj(ax, idir, gline, xlim, ylim):

    eps = (abs(xlim[0] - xlim[1]) + abs(ylim[0] - ylim[1]))**2 * 1e-12
    #print(eps)

    x1 = y1 = x2 = y2 = []
    if   idir == 0:
        x1 = gline[:, 0, 1]
        y1 = gline[:, 0, 2]
        x2 = gline[:, 1, 1]
        y2 = gline[:, 1, 2]
    elif idir == 1:
        x1 = gline[:, 0, 2]
        y1 = gline[:, 0, 0]
        x2 = gline[:, 1, 2]
        y2 = gline[:, 1, 0]
    elif idir == 2:
        x1 = gline[:, 0, 0]
        y1 = gline[:, 0, 1]
        x2 = gline[:, 1, 0]
        y2 = gline[:, 1, 1]

    ias = \
        ((x1 - xlim[0]) * (x1 - xlim[1]) <= eps) & \
        ((x2 - xlim[0]) * (x2 - xlim[1]) <= eps) & \
        ((y1 - ylim[0]) * (y1 - ylim[1]) <= eps) & \
        ((y2 - ylim[0]) * (y2 - ylim[1]) <= eps)
    #print(ias)
    #print(gline)

    ax.plot([x1[ias], x2[ias]], [y1[ias], y2[ias]], color='k', lw=1)

# (private) near2d.log
def _log_n2d(fp, iplane, freq, compo, idir, p0, x, y, fc):
    # header
    fp.write('#%d : frequency[Hz] = %.3e\n' % (iplane + 1, freq))
    fp.write('  No.  No.     X[m]        Y[m]        Z[m]       ')
    if   compo.startswith('E'):
        fp.write('E[V/m]      Ex[V/m]   Ex[deg]    Ey[V/m]   Ey[deg]    Ez[V/m]   Ez[deg]\n')
    elif compo.startswith('H'):
        fp.write('H[A/m]      Hx[A/m]   Hx[deg]    Hy[A/m]   Hy[deg]    Hz[A/m]   Hz[deg]\n')

    # body
    p = [0] * 3
    fmt1 = '%5d%5d%12.3e%12.3e%12.3e'
    fmt2 = '%12.4e%12.4e%9.3f%12.4e%9.3f%12.4e%9.3f\n'
    for n1, p1 in enumerate(x):
        for n2, p2 in enumerate(y):
            if   idir == 0:
                p = [p0, p1, p2]
            elif idir == 1:
                p = [p2, p0, p1]
            elif idir == 2:
                p = [p1, p2, p0]
            f = fc[n2, n1, :]
            fp.write(fmt1 % (n1, n2, p[0], p[1], p[2]))
            if   compo.startswith('E'):
                fp.write(fmt2 % (f[0], f[1], f[4], f[2], f[5], f[3], f[6]))
            elif compo.startswith('H'):
                fp.write(fmt2 % (f[7], f[8], f[11], f[9], f[12], f[10], f[13]))
