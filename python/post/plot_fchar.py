# -*- coding: utf-8 -*-
"""
plot_fchar.py
"""

import math
import numpy as np
import matplotlib.pyplot as plt
import sol.farfield

EPS = 1e-12

# 周波数特性
def plot(Post, Freq1, Freq2, Zin, Ref, Pin, Spara, Coupling, fSurface, cSurface):

    if Post['smith'] > 0:
        _smith(Post, Freq1, Zin)

    if Post['zin'][0] > 0:
        _fchar(1, Post, Freq1, Freq2, Zin, Ref, Pin, Spara, Coupling, fSurface, cSurface)

    if (Post['yin'][0] > 0):
        _fchar(2, Post, Freq1, Freq2, Zin, Ref, Pin, Spara, Coupling, fSurface, cSurface)

    if (Post['ref'][0] > 0):
        _fchar(3, Post, Freq1, Freq2, Zin, Ref, Pin, Spara, Coupling, fSurface, cSurface)

    if (Post['spara'][0] > 0):
        _fchar(4, Post, Freq1, Freq2, Zin, Ref, Pin, Spara, Coupling, fSurface, cSurface)

    if (Post['coupling'][0] > 0):
        _fchar(5, Post, Freq1, Freq2, Zin, Ref, Pin, Spara, Coupling, fSurface, cSurface)

    if (Post['f0d'][2] > 0):
        _fchar(6, Post, Freq1, Freq2, Zin, Ref, Pin, Spara, Coupling, fSurface, cSurface)

# (private)　周波数特性
def _fchar(itype, Post, Freq1, Freq2, Zin, Ref, Pin, Spara, Coupling, fSurface, cSurface):

    nfreq = len(Freq1) if ((itype >= 1) and (itype <= 5)) else len(Freq2)
    nfeed = Zin.shape[0]
    npoint = Spara.shape[0]
    #print(itype, nfreq, nfeed, npoint)

    if nfreq < 2:
        print('*** single frequency.')
        return
    if (itype < 1) or (itype > 6):
        return

    name = ['input impedance', 'input admittance', 'reflection', 'S-parameter', 'coupling', 'far field']
    #comp = [['Rin', 'Xin'], ['Gin', 'Bin'], ['', ''], ['', ''], ['', ''], ['', '']]
    unit = ['[ohm]', '[mS]', '[dB]', '[dB]', '[dB]', '']

    # open far0d.log
    if itype == 6:
        fname = 'far0d.log'
        fp = open(fname, 'wt', encoding='utf-8')

    # ページ数
    npage = 0
    if   (itype == 1) or (itype == 2) or (itype == 3):
        npage = nfeed
    elif (itype == 4):
        npage = 1 if (npoint > 0) else 0
    elif (itype == 5):
        npage = nfeed
    elif (itype == 6):
        npage = 1

    x = y = []
    color = []
    label = []
    #print(itype, nfeed, npoint, nfreq, npage)

    for ipage in range(npage):
        # 描画データ
        if   itype == 1:
            # 入力インピーダンス(実部虚部の2本)
            ifeed = ipage
            x = Freq1 * Post['fscale']
            y = np.zeros((2, nfreq), 'f8')
            y[0, :] = Zin[ifeed, :].real
            y[1, :] = Zin[ifeed, :].imag
            label = ['Rin', 'Xin']
            color = ['r', 'b']
        elif itype == 2:
            # 入力アドミッタンス(実部虚部の2本)
            ifeed = ipage
            x = Freq1 * Post['fscale']
            y = np.zeros((2, nfreq), 'f8')
            y[0, :] = (1e3 / Zin[ifeed, :]).real
            y[1, :] = (1e3 / Zin[ifeed, :]).imag
            label = ['Gin', 'Bin']
            color = ['r', 'b']
        elif itype == 3:
            # 反射係数(1本)
            ifeed = ipage
            x = Freq1 * Post['fscale']
            y = np.zeros((1, nfreq), 'f8')
            y[0, :] = 20 * np.log10(np.maximum(np.abs(Ref[ifeed, :]), EPS))
            label = ['']
            color = ['k']
        elif itype == 4:
            # Sパラメーター(観測点の数)
            x = Freq1 * Post['fscale']
            y = np.zeros((npoint, nfreq), 'f8')
            for ipoint in range(npoint):
                y[ipoint, :] = 20 * np.log10(np.maximum(np.abs(Spara[ipoint, :]), EPS))
            label = [''] * npoint
            color = ['k'] * npoint
            if npoint > 0:
                color[0] = 'r'
            if npoint > 1:
                color[1] = 'b'
            if npoint > 2:
                color[2] = 'g'
        elif itype == 5:
            # 結合度(観測点の数)
            ifeed = ipage
            x = Freq1 * Post['fscale']
            y = np.zeros((npoint, nfreq), 'f8')
            for ipoint in range(npoint):
                y[ipoint, :] = 20 * np.log10(np.maximum(np.abs(Coupling[ifeed, ipoint, :]), EPS))
            label = [''] * npoint
            color = ['k'] * npoint
            if npoint > 0:
                color[0] = 'r'
            if npoint > 1:
                color[1] = 'b'
            if npoint > 2:
                color[2] = 'g'
        elif itype == 6:
            # 遠方界(1本)
            x = Freq2 * Post['fscale']
            y = np.zeros((1, nfreq), 'f8')
            theta = Post['f0d'][0]
            phi   = Post['f0d'][1]
            C     = Post['C']
            ETA0  = Post['ETA0']

            # log (header)
            _log_f0d_1(fp, Post)

            for ifreq, freq in enumerate(Freq2):
                # 波数と遠方界係数
                kwave = (2 * math.pi * freq) / C
                ffctr = sol.farfield.factor(ifreq, kwave, Pin, nfeed, Post['mloss'], ETA0)
                # 遠方界
                etheta, ephi = sol.farfield.calc(ifreq, theta, phi, \
                    ffctr, kwave, fSurface, cSurface)
                efar = sol.farfield.farComponent(etheta, ephi)
                y[0, ifreq] = 20 * np.log10(max(efar[0], EPS))
                # log (body)
                _log_f0d_2(fp, ifreq, freq, efar)

            label = ['']
            color = ['k']

        # figure
        strfig = 'OpenRCWA - %s (%d/%d)' % (name[itype - 1], ipage + 1, npage)
        fig = plt.figure(strfig, figsize=(Post['w2d'][0], Post['w2d'][1]))
        ax = fig.add_subplot()

        # plot
        for idata in range(y.shape[0]):
            ax.plot(x, y[idata, :], color=color[idata], label=label[idata])
        ax.grid()
        
        # X軸
        ax.set_xlim(x[0], x[-1])
        if nfreq > 1:
            fdiv = np.linspace(x[0], x[-1], Post['freqdiv'] + 1)
        else:
            fdiv = x[0]
        ax.set_xticks(fdiv)

        # X軸ラベル
        ax.set_xlabel('frequency ' + Post['funit'])

        # Y軸スケール(指定スケール)
        ymin = ymax = ydiv = 0
        if   (itype == 1) and (Post['zin'][0] == 2):
            [ymin, ymax, ydiv] = Post['zin'][1:4]
        elif (itype == 2) and (Post['yin'][0] == 2):
            [ymin, ymax, ydiv] = Post['yin'][1:4]
        elif (itype == 3) and (Post['ref'][0] == 2):
            [ymin, ymax, ydiv] = Post['ref'][1:4]
        elif (itype == 4) and (Post['spara'][0] == 2):
            [ymin, ymax, ydiv] = Post['spara'][1:4]
        elif (itype == 5) and (Post['coupling'][0] == 2):
            [ymin, ymax, ydiv] = Post['coupling'][1:4]
        elif (itype == 6) and (Post['f0d'][2] == 2):
            [ymin, ymax, ydiv] = Post['f0d'][3:6]
        if (ymin < ymax):
            ax.set_ylim(ymin, ymax)
            ax.set_yticks(np.linspace(ymin, ymax, ydiv + 1))

        # Y軸ラベル
        if   (itype == 1) or (itype == 2) or (itype == 3) or (itype == 4) or (itype == 5):
            strunit = unit[itype - 1]
        elif (itype == 6):
            name[itype - 1] = Post['farname']
            strunit = Post['f0dunit']
        ax.set_ylabel(name[itype - 1] + ' ' + strunit)

        # XY範囲を取得
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()

        # Y=0線を描く
        if (ylim[0] < 0) and (ylim[1] > 0):
            ax.plot([x[0], x[-1]], [0, 0], color='gray')
        
        # 見本
        if (itype == 1) or (itype == 2):
            ax.legend(loc='best')

        # タイトル
        title = Post['title']
        if   (itype == 1) or (itype == 2):
            strtitle = '%s\n(feed# = %d)' % (title, ifeed + 1)
        elif (itype == 3):
            strtitle = '%s\nfeed# = %d, Z0 = %gohm, min = %.3fdB' % (title, ifeed + 1, Post['z0'][ifeed], np.min(y))
        elif (itype == 4):
            strtitle = '%s\nmax = %.3fdB, min = %.3fdB' % (title, np.max(y), np.min(y))
        elif (itype == 5):
            strtitle = '%s\nmax = %.3fdB, min = %.3fdB' % (title, np.max(y), np.min(y))
        elif (itype == 6):
            strtitle = '%s\n(theta, phi) = (%g, %g)deg, max = %.3f%s' % (title, Post['f0d'][0], Post['f0d'][1], np.max(y), strunit)
        ax.set_title(strtitle)

        # 右Y軸表示
        if itype == 3:
            # VSWR
            vs = np.array([1.1, 1.2, 1.5, 2.0, 3.0, 5.0, 10.0])
            vsdb = 20 * np.log10((vs - 1) / (vs + 1))
            lx = xlim[1] - xlim[0]
            ly = ylim[1] - ylim[0]
            ax.text(xlim[1], ylim[1] + 0.01 * ly, 'VSWR', color='r')
            for ivs, db in enumerate(vsdb):
                if (db >= ylim[0]) and (db <= ylim[1]):
                    xvs = [xlim[1], xlim[1] - 0.015 * lx]
                    yvs = [db, db]
                    ax.plot(xvs, yvs, color='r')
                    ax.text(xlim[1], db + 0.01 * ly, str(vs[ivs]), color='r')
        elif itype == 4:
            # Sパラメーター
            for ipoint in range(npoint):
                ax.text(x[-1], y[ipoint, -1], 'S%d1' % (ipoint + 1), color=color[ipoint])
        elif itype == 5:
            # 結合度
            for ipoint in range(npoint):
                ax.text(x[-1], y[ipoint, -1], 'C%d%d' % (ipoint + 1, ifeed + 1), color=color[ipoint])

# (private) smith chart
def _smith(Post, Freq1, Zin):

    nfeed = Zin.shape[0]
    nfreq = Zin.shape[1]

    if (nfeed < 1) or (nfreq < 1):
        return

    # Zin
    for ifeed in range(nfeed):
        # figure
        strfig = 'OpenRCWA - SmithChart (%d/%d)' % (ifeed + 1, nfeed)
        fig = plt.figure(strfig, figsize=(Post['w2d'][0], Post['w2d'][1]))
        ax = fig.add_subplot()

        z0 = Post['z0'][ifeed]

        # Smith Chart
        _smithchart(ax, z0)

        ref = (Zin[ifeed, :] - z0) / (Zin[ifeed, :] + z0)
        ax.plot(ref.real, ref.imag, color='k', marker='o', markersize=3)

        # #n : start frequency
        #ax.text(ref[0].real, ref[0].imag, '#' + str(ifeed + 1))
        ax.text(ref[0].real, ref[0].imag, 'f0', color='r')
        ax.text(ref[-1].real, ref[-1].imag, 'f1', color='b')

        # title
        ax.set_title('%s\nf%s = %g - %g, feed#%d, Z0 = %gohm' %
            (Post['title'], Post['funit'], Freq1[0] * Post['fscale'], Freq1[-1] * Post['fscale'], ifeed + 1, z0))

        # layout
        ax.set_aspect('equal')
        ax.axis('off')

# (private) Smith Chart
def _smithchart(ax, z0):

    dh = 0.06 # string shift
    
    # R=0
    div = 50
    angle = np.linspace(0, 2 * np.pi, div + 1)
    ax.plot(np.cos(angle), np.sin(angle), color='gray')

    # X=0
    ax.plot([-1, +1], [0, 0], color='gray', linestyle='--')
    
    # DB
    const = [0.2, 0.5, 1, 2, 5]
    zmin  = [0,    1.1, 2.2, 3.5,  5,  20,  200]
    zmax  = [1,    2,   3,   4,   10, 100, 1000]
    zdiv  = [20,   9,   4,   1,    5,   8,    4]
    
    # R-const
    for ir in range(len(const)):
        xlist = np.array([])
        for ix in range(len(zmin)):
            xlist = np.hstack((xlist, np.linspace(zmin[ix], zmax[ix], zdiv[ix] + 1)))

        z = const[ir] + (1j * xlist)
        g = (z - 1) / (z + 1)
        ax.plot(g.real, +g.imag, color='gray', linestyle='--')
        ax.plot(g.real, -g.imag, color='gray', linestyle='--')
        strr = '%d' % (z0 * const[ir])
        ax.text(g[0].real - dh, 0, strr)

    # X-const
    for ix in range(len(const)):
        rlist = np.array([])
        for ir in range(len(zmin)):
            rlist = np.hstack((rlist, np.linspace(zmin[ir], zmax[ir], zdiv[ir] + 1)))

        z = rlist + (1j * const[ix])
        g = (z - 1) / (z + 1)
        ax.plot(g.real, +g.imag, color='gray', linestyle='--')
        ax.plot(g.real, -g.imag, color='gray', linestyle='--')
        strx = '%dj' % (z0 * const[ix])
        ax.text(g[0].real - dh, +g[0].imag, '+' + strx)
        ax.text(g[0].real - dh, -g[0].imag, '-' + strx)

# (private) far0d.log (header)
def _log_f0d_1(fp, Post):

    fp.write('theta=%.3f[deg] phi=%.3f[deg]\n' % (Post['f0d'][0], Post['f0d'][1]))
    fp.write('  No. frequency[Hz]    E-abs[dB]  E-theta[dB]    E-phi[dB]  E-major[dB]  E-minor[dB]   E-RHCP[dB]   E-LHCP[dB] AxialRatio[dB]\n')

# (private) far0d.log (body)
def _log_f0d_2(fp, ifreq, freq, efar):

    fmt = '%4d%15.5e%13.4f%13.4f%13.4f%13.4f%13.4f%13.4f%13.4f%13.4f\n'
    fdb = 20 * np.log10(np.maximum(efar, EPS))
    fp.write(fmt % (ifreq, freq, fdb[0], fdb[1], fdb[2], fdb[3], fdb[4], fdb[5], fdb[6], fdb[3] - fdb[4]))
