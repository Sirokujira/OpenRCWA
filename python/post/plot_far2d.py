# -*- coding: utf-8 -*-
"""
plot_far2d.py
"""

import math
import numpy as np
import matplotlib.pyplot as plt
import sol.farfield

EPS = 1e-20

# far field (3D)
def plot(Post, Freq2, Pin, fSurface, cSurface, gline):
    
    if (Post['f2d'][0] < 1) or (Post['f2d'][1] < 2) or (Post['f2d'][2] < 2):
        return

    nfreq = len(Freq2)
    nfeed = Pin.shape[1]
    #print(nfreq, nfeed)
 
    # log
    fname = 'far2d.log'
    fp = open(fname, 'wt', encoding='utf-8')

    # far field factor
    kwave = np.zeros(nfreq, 'f8')
    ffctr = np.zeros(nfreq, 'f8')
    for ifreq, freq in enumerate(Freq2):
        kwave[ifreq] = (2 * math.pi * freq) / Post['C']
        ffctr[ifreq] = sol.farfield.factor(ifreq, kwave[ifreq], Pin, nfeed, Post['mloss'], Post['ETA0'])
 
    # alloc
    nth = Post['f2d'][1] + 1
    nph = Post['f2d'][2] + 1
    pfar = np.zeros((nth, nph, 7), 'f8')
    th = np.linspace(0, 180, nth)
    ph = np.linspace(0, 360, nph)
    pfar = np.zeros((nth, nph, 7), 'f8')

    # 物体形状の大きさ
    gxmin = np.min(gline[:, :, 0])
    gxmax = np.max(gline[:, :, 0])
    gymin = np.min(gline[:, :, 1])
    gymax = np.max(gline[:, :, 1])
    gzmin = np.min(gline[:, :, 2])
    gzmax = np.max(gline[:, :, 2])
    gx0 = (gxmin + gxmax) / 2
    gy0 = (gymin + gymax) / 2
    gz0 = (gzmin + gzmax) / 2
    gfctr = Post['f2dobj'] / max(gxmax - gxmin, gymax - gymin, gzmax - gzmin)
    #print(gxmin, gxmax, gymin, gymax, gzmin, gzmax, gx0, gy0, gz0, gfctr)

    # plot
    nfig = 0
    for ifreq, freq in enumerate(Freq2):
        # 遠方界電力
        for ith in range(nth):
            for iph in range(nph):
                etheta, ephi = sol.farfield.calc(ifreq, th[ith], ph[iph], \
                    ffctr[ifreq], kwave[ifreq], fSurface, cSurface)
                efar = sol.farfield.farComponent(etheta, ephi)
                pfar[ith, iph, 0:7] = efar**2

        # 値が0である成分は飛ばす
        for m in range(7):
            if Post['f2dcompo'][m] == 1:
                if (np.max(np.max(pfar[:, :, m])) < EPS):
                    print('*** %s : max = 0' % Post['farcomp'][m])
                    Post['f2dcompo'][m] = 0

        # 統計量
        stat = _statistics(Post, nfeed, nth, nph, pfar)

        # log
        _log_f2d(fp, freq, nth, nph, th, ph, pfar)

        # dBに変換
        if Post['f2ddb'] == 1:
            pfar = 10 * np.log10(np.maximum(pfar, EPS))

        # 1成分=1ページ
        for m in range(7):
            if Post['f2dcompo'][m] == 0:
                continue

            # 最大値とその方向
            rmax = -float('inf')
            thmax = 0
            phmax = 0
            for ith in range(nth):
                for iph in range(nph):
                    if (pfar[ith, iph, m] > rmax):
                        rmax = pfar[ith, iph, m]
                        thmax = th[ith]
                        phmax = ph[iph]
            #print(thmax, phmax)

            # 最小値
            if Post['f2ddb'] == 1:
                rmin = rmax - abs(Post['f2dscale'][2] - Post['f2dscale'][1])
            else:
                rmin = 0

            # figure
            nfig += 1
            strfig = 'OpenRCWA - far field (3D) (%d/%d)' % (nfig, nfreq * sum(Post['f2dcompo']))
            fig = plt.figure(strfig, figsize=(Post['w3d'][0], Post['w3d'][1]))
            ax = fig.add_subplot(projection='3d')

            # plot
            ph2d, th2d = np.meshgrid(np.deg2rad(ph), np.deg2rad(th))
            r = np.maximum((pfar[:, :, m] - rmin) / (rmax - rmin), 0)
            x = r * np.cos(ph2d) * np.sin(th2d)
            y = r * np.sin(ph2d) * np.sin(th2d)
            z = r                * np.cos(th2d)
            ax.plot_surface(x, y, z)#, cmap='rainbow')

            # 物体形状の上書き
            if (Post['f2dobj'] > 0):
                for n in range(gline.shape[0]):
                    ax.plot(gfctr * (gline[n, :, 0] - gx0), \
                            gfctr * (gline[n, :, 1] - gy0), \
                            gfctr * (gline[n, :, 2] - gz0), color='gray')

            # title
            str2 = '%s, f = %.3g%s' % (Post['farcomp'][m], freq * Post['fscale'], Post['funit'])
            str4 = 'max = %.4g%s @ (theta, phi) = (%.1f, %.1f)deg' % (rmax, Post['f2dunit'], thmax, phmax)
            ax.set_title('%s\n%s, %s\n%s' % (Post['title'], str2, stat, str4))

            # label
            ax.set_xlabel('X')
            ax.set_ylabel('Y')
            ax.set_zlabel('Z')

            # 視点
            ax.view_init(elev = 90 - Post['w3d'][2], azim = Post['w3d'][3], roll = 0)

            # アスペクト比=1
            ax.set_aspect('equal')

# (private) statistics
def _statistics(Post, nfeed, nth, nph, pfar):
    stat = ''

    pmax = 0
    psum = 0
    fctr = (math.pi / (nth - 1)) * (2 * math.pi / (nph - 1)) / (4 * math.pi)
    for ith in range(nth):
        for iph in range(nph - 1):
            th = np.pi * ith / (nth - 1)
            pabs = pfar[ith, iph, 0]
            psum += fctr * math.sin(th) * pabs
            pmax = max(pabs, pmax)

    if nfeed > 0:
        # 給電点
        pmax = pmax / psum
        if Post['f2ddb'] == 1:
            pmaxdb = 10 * math.log10(max(pmax, EPS))
            stat = 'directive gain = %.3fdBi' % pmaxdb
        else:
            stat = 'directive gain = %.3f' % pmax
        #stat2 = 'efficiency = %.3f[%%]' % (psum * 100)
    else:
        # 平面波入射
        if Post['f2ddb'] == 1:
            psumdb = 10 * math.log10(max(psum, EPS))
            stat = 'total cross section = %.3fdBm^2' % psumdb
        else:
            stat = 'total cross section = %.3em^2' % psum

    return stat

# (private) far2d.log
def _log_f2d(fp, freq, nth, nph, th, ph, pfar):
    # header
    fp.write('frequency[Hz] = %.3e\n' % freq)
    fp.write(' No. No. theta[deg] phi[deg]   E-abs[dB]  E-theta[dB]    E-phi[dB]  E-major[dB]  E-minor[dB]   E-RHCP[dB]   E-LHCP[dB] AxialRatio[dB]\n')

    # body
    fmt = '%4d%4d %9.1f%9.1f%13.4f%13.4f%13.4f%13.4f%13.4f%13.4f%13.4f%13.4f\n'
    for ith in range(nth):
        for iph in range(nph):
            pdb = 10 * np.log10(np.maximum(pfar[ith, iph, :], EPS))
            fp.write(fmt % \
                (ith, iph, th[ith], ph[iph], pdb[0], pdb[1], pdb[2], pdb[3], pdb[4], pdb[5], pdb[6], pdb[3] - pdb[4]))
