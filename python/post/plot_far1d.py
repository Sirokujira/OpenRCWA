# -*- coding: utf-8 -*-
"""
plot_far1d.py
"""

import math
import numpy as np
import matplotlib.pyplot as plt
import sol.farfield

EPS = 1e-20

def plot(Post, Freq2, Pin, fSurface, cSurface):

    nplane = len(Post['f1ddiv'])
    nfreq = len(Freq2)
    nfeed = Pin.shape[1]
    #print(nplane, nfreq, nfeed)
    
    if (nplane < 1) or (nfreq < 1):
        return

    # log
    fname = 'far1d.log'
    fp = open(fname, 'wt', encoding='utf-8')

    # far field factor
    kwave = np.zeros(nfreq, 'f8')
    ffctr = np.zeros(nfreq, 'f8')
    for ifreq, freq in enumerate(Freq2):
        kwave[ifreq] = (2 * math.pi * freq) / Post['C']
        ffctr[ifreq] = sol.farfield.factor(ifreq, kwave[ifreq], Pin, nfeed, Post['mloss'], Post['ETA0'])

    # color
    col = ['k', 'r', 'b', 'r', 'b', 'r', 'b']

    # label
    strlabel = ['E-abs', 'E-theta', 'E-phi', 'E-major', 'E-minor', 'E-RHCP', 'E-LHCP']

    # plot
    nfig = 0
    for ip in range(nplane):
        adir = Post['f1ddir'][ip]
        adiv = Post['f1ddiv'][ip]
        pfar = np.zeros((adiv + 1, 7), 'f8')
        angle = np.linspace(0, 360, adiv + 1)

        for ifreq, freq in enumerate(Freq2):

            # set data
            for ia in range(adiv + 1):
                # direction
                theta = 0
                phi = 0
                strplane = ''
                if   adir == 'X':
                    theta = angle[ia]
                    phi   = 90
                    strplane = 'phi = 90deg'
                    strx = 'theta'
                    theta_direction = 'clockwise'
                    theta_zero_location = 'N'
                elif adir == 'Y':
                    theta = angle[ia]
                    phi   = 0
                    strplane = 'phi = 0deg'
                    strx = 'theta'
                    theta_direction = 'clockwise'
                    theta_zero_location = 'N'
                elif adir == 'Z':
                    theta = 90
                    phi   = angle[ia]
                    strplane = 'theta = 90deg'
                    strx = 'phi'
                    theta_direction = 'counterclockwise'
                    theta_zero_location = 'E'
                elif adir == 'V':
                    theta = angle[ia]
                    phi   = Post['f1dangle'][ip]
                    strplane = 'phi=%gdeg' % phi
                    strx = 'theta'
                    theta_direction = 'clockwise'
                    theta_zero_location = 'N'
                elif adir == 'H':
                    theta = Post['f1dangle'][ip]
                    phi   = angle[ia]
                    strplane = 'theta=%gdeg' % theta
                    strx = 'phi'
                    theta_direction = 'counterclockwise'
                    theta_zero_location = 'E'

                # far field
                etheta, ephi = sol.farfield.calc(ifreq, theta, phi, \
                    ffctr[ifreq], kwave[ifreq], fSurface, cSurface)
                efar = sol.farfield.farComponent(etheta, ephi)
                pfar[ia, 0:7] = efar**2

            # log
            _log_f1d(fp, ip, adir, adiv, freq, Post['f1dangle'][ip], pfar)

            # normalization
            if Post['f1dnorm'] == 1:
                pmax = np.max(pfar[:, 0])
                pfar[:, 0:7] /= pmax

            # to dB
            if Post['f1ddb'] == 1:
                pfar = 10 * np.log10(np.maximum(pfar, EPS))

            # scale
            dmax = np.max(pfar)
            if Post['f1dscale'][0] == 0:
                # auto scale
                if Post['f1ddb'] == 1:
                    ymax = dmax
                    ymin = dmax - 40
                else:
                    ymax = dmax
                    ymin = 0
            else:
                # user scale
                ymin = Post['f1dscale'][1]
                ymax = Post['f1dscale'][2]

            # plot
            for icomp in range(3):
                if Post['f1dcompo'][icomp] == 0:
                    continue

                # component
                if   icomp == 0:
                    ic = [1, 1, 1, 0, 0, 0, 0]
                elif icomp == 1:
                    ic = [1, 0, 0, 1, 1, 0, 0]
                elif icomp == 2:
                    ic = [1, 0, 0, 0, 0, 1, 1]

                # figure
                nfig += 1
                strfig = 'OpenRCWA - far field (2D) (%d/%d)' % (nfig, nplane * nfreq * np.sum(Post['f1dcompo']))
                fig = plt.figure(strfig, figsize=(Post['w2d'][0], Post['w2d'][1]))
                if Post['f1dstyle'] == 0:
                    ax = fig.add_subplot(projection='polar')
                else:
                    ax = fig.add_subplot()

                # plot
                for m in range(7):
                    if ic[m] == 0:
                        continue
                    if Post['f1dstyle'] == 0:
                        ax.plot(np.deg2rad(angle), np.maximum(pfar[:, m], ymin), label=strlabel[m], color=col[m])
                    else:
                        ax.plot(angle, pfar[:, m], label=strlabel[m], color=col[m])
                        ax.grid()
                ax.legend(loc='best')
                #ax.legend(loc='upper right')

                # axis
                if Post['f1dstyle'] == 0:
                    # polar plot
                    ax.set_rlim(ymin, ymax)
                    if Post['f1dscale'][0] == 1:
                        rdiv = np.linspace(ymin, ymax, Post['f1dscale'][3] + 1)
                        ax.set_rticks(rdiv)
                    ax.set_theta_direction(theta_direction)
                    ax.set_theta_zero_location(theta_zero_location)
                else:
                    # XY plot
                    # X-axis
                    ax.set_xlim(0, 360)
                    xdiv = np.linspace(0, 360, 7)
                    ax.set_xticks(xdiv)
                    ax.set_xlabel(strx + ' [deg]')
                    # Y-axis
                    ax.set_ylim(ymin, ymax)
                    if Post['f1dscale'][0] == 1:
                        ydiv = np.linspace(ymin, ymax, Post['f1dscale'][3] + 1)
                        ax.set_yticks(ydiv)
                    ax.set_ylabel(Post['farname'] + ' ' + Post['f1dunit'])

                # title
                iamax = np.argmax(pfar[:, m])
                ax.set_title('%s\n%s, f = %.3g%s, max = %.4g%s @ %gdeg' % \
                    (Post['title'], strplane, freq * Post['fscale'], Post['funit'], dmax, Post['f1dunit'], angle[iamax]))

# (private) far1d.log
def _log_f1d(fp, ip, adir, adiv, freq, angle0, pfar):

    # header
    fp.write('#%d : %s-plane' % (ip + 1, adir))
    if   adir == 'V':
        fp.write(' (phi = %.2fdeg)' % angle0)
    elif adir == 'H':
        fp.write(' (theta = %.2fdeg)' % angle0)
    fp.write(', frequency[Hz] = %.3e\n' % freq)
    fp.write('  No.   deg    E-abs[dB]  E-theta[dB]    E-phi[dB]  E-major[dB]  E-minor[dB]   E-RHCP[dB]   E-LHCP[dB] AxialRatio[dB]\n')

    # body
    fmt = '%4d%7.1f%13.4f%13.4f%13.4f%13.4f%13.4f%13.4f%13.4f%13.4f\n'
    for ia in range(adiv + 1):
        angle = 360 * ia / adiv
        fdb = 10 * np.log10(np.maximum(pfar[ia], EPS))
        fp.write(fmt % \
            (ia, angle, fdb[0], fdb[1], fdb[2], fdb[3], fdb[4], fdb[5], fdb[6], fdb[3] - fdb[4]))
