# -*- coding: utf-8 -*-
"""
plot_iter.py
"""

import numpy as np
import matplotlib.pyplot as plt
import sol.dft

# 収束状況、時間波形。スペクトル
def plot(Post, Freq1, Eiter, Hiter, Iiter, Niter, VFeed, IFeed, VPoint, Ntime):

    if Post['iter'] > 0:
        _average(Post, Eiter, Hiter, Iiter, Niter)

    if Post['feed'] > 0:
        _feed_waveform(Post, VFeed, IFeed, Ntime)
        _feed_spectrum(Post, Freq1, VFeed, IFeed, Ntime)

    if Post['point'] > 0:
        _point_waveform(Post, VPoint, Ntime)
        _point_spectrum(Post, Freq1, VPoint, Ntime)

# (private) 収束状況:E/H平均電磁界 (2D)
def _average(Post, Eiter, Hiter, Iiter, Niter):

    # figure
    strfig = 'OpenRCWA - iteration (2D)'
    fig = plt.figure(strfig, figsize=(Post['w2d'][0], Post['w2d'][1]))
    ax = fig.add_subplot()

    # データ
    n = Niter
    x  = Iiter[0: n]
    y1 = Eiter[0: n]
    y2 = Hiter[0: n]

    # 最大・最小
    ymax = max(max(y1), max(y2))
    ymin = ymax * 1e-6
    y1 = np.maximum(y1, ymin)
    y2 = np.maximum(y2, ymin)

    # plot
    ax.semilogy(x, y1, color='r', label='<E>')
    ax.semilogy(x, y2, color='b', label='<H>')
    ax.legend(loc='best')

    # X軸
    ax.set_xlim(x[0], x[-1])
    ax.set_xlabel('time step')

    # Y軸
    ax.set_ylim(ymin, ymax)
    ax.set_ylabel('E/H average')

    # タイトル
    ax.set_title(Post['title'])

    # 罫線ON
    ax.grid(True)

# (private) 給電点の時間波形 (2D)
def _feed_waveform(Post, VFeed, IFeed, Ntime):

    dt = Post['dt']

    # X軸 : 時間
    x = np.linspace(0, (Ntime - 1) * dt, Ntime)

    for ifeed in range(VFeed.shape[0]):
        # figure
        strfig = 'OpenRCWA - waveform feed #%d (2D)' % (ifeed + 1)
        fig = plt.figure(strfig, figsize=(Post['w2d'][0], Post['w2d'][1]))
        ax1 = fig.add_subplot()

        # V　左目盛り
        y1 = abs(VFeed[ifeed, 0: Ntime])
        ax1.semilogy(x, y1, color='r')
        y1max = max(y1)
        ax1.set_ylim([y1max * 1e-6, y1max])
        ax1.set_ylabel('V(t) [V]', color='r')

        # I 右目盛り
        ax2 = ax1.twinx()
        y2 = abs(IFeed[ifeed, 0: Ntime])
        ax2.semilogy(x, y2, color='b')
        y2max = max(y2)
        ax2.set_ylim([y2max * 1e-6, y2max])
        ax2.set_ylabel('I(t) [A]', color='b')

        ax1.grid(True)

        # X軸
        ax1.set_xlabel('time [sec]')
        ax1.set_xlim(x[0], x[-1])

        # タイトル
        strtitle = "%s\nwaveform feed #%d" % (Post['title'], ifeed + 1)
        ax1.set_title(strtitle)

# (private) 給電点のスペクトル (2D)
def _feed_spectrum(Post, Freq1, VFeed, IFeed, Ntime):

    if len(Freq1) <= 1:
        print("*** frequency1 : single frequency")
        return

    dt = Post['dt']

    # X軸 : 周波数
    x = Freq1 * Post['fscale']

    for ifeed in range(VFeed.shape[0]):
        # figure
        strfig = 'OpenRCWA - spectrum feed #%d (2D)' % (ifeed + 1)
        fig = plt.figure(strfig, figsize=(Post['w2d'][0], Post['w2d'][1]))
        ax = fig.add_subplot()

        # V/I
        y1 = np.zeros(len(Freq1), 'f8')
        y2 = np.zeros(len(Freq1), 'f8')
        for ifreq in range(len(Freq1)):
            y1[ifreq] = abs(sol.dft.calc(Ntime, VFeed[ifeed, 0: Ntime], Freq1[ifreq], dt, 0))
            y2[ifreq] = abs(sol.dft.calc(Ntime, IFeed[ifeed, 0: Ntime], Freq1[ifreq], dt, 0))
        # 最大値で正規化
        y1 /= max(y1)
        y2 /= max(y2)

        # plot
        ax.plot(x, y1, color='r', label='V(f)')
        ax.plot(x, y2, color='b', label='I(f)')
        ax.legend(loc='best')
        ax.grid(True)

        # Y軸
        ax.set_ylim([0, 1])

        # X軸
        ax.set_xlabel('frequency %s' % Post['funit'])
        ax.set_xlim(x[0], x[-1])
        ax.set_xticks(np.linspace(x[0], x[-1], Post['freqdiv'] + 1))

        # タイトル
        strtitle = "%s\nspectrum feed #%d" % (Post['title'], ifeed + 1)
        ax.set_title(strtitle)

# (private) 観測点の時間波形 (2D)
def _point_waveform(Post, VPoint, Ntime):

    dt = Post['dt']
    npoint = VPoint.shape[0] - (2 if VPoint.shape[0] > 2 else 0)

    # X軸 : 時間
    x = np.linspace(0, (Ntime - 1) * dt, Ntime)

    for ipoint in range(npoint):
        # figure
        strfig = 'OpenRCWA - waveform point #%d (2D)' % (ipoint + 1)
        fig = plt.figure(strfig, figsize=(Post['w2d'][0], Post['w2d'][1]))
        ax = fig.add_subplot()

        # V
        y = abs(VPoint[ipoint, 0: Ntime])

        # plot
        ax.semilogy(x, y)
        ax.grid(True)

        # Y軸
        ymax = max(y)
        ax.set_ylim([ymax * 1e-6, ymax])
        ax.set_ylabel('V(t) [V]')

        # X軸
        ax.set_xlabel('time [sec]')
        ax.set_xlim(x[0], x[-1])

        # タイトル
        strtitle = "%s\nwaveform point #%d" % (Post['title'], ipoint + 1)
        ax.set_title(strtitle)

# (private) 観測点のスペクトル (2D)
def _point_spectrum(Post, Freq1, VPoint, Ntime):

    if len(Freq1) <= 1:
        print("*** frequency1 : single frequency")
        return

    dt = Post['dt']
    npoint = VPoint.shape[0] - (2 if VPoint.shape[0] > 2 else 0)

    # X軸 : 周波数
    x = Freq1 * Post['fscale']

    for ipoint in range(npoint):
        # figure
        strfig = 'OpenRCWA - spectrum point #%d (2D)' % (ipoint + 1)
        fig = plt.figure(strfig, figsize=(Post['w2d'][0], Post['w2d'][1]))
        ax = fig.add_subplot()

        # スペクトル
        y = np.zeros(len(Freq1), 'f8')
        for ifreq, freq in enumerate(Freq1):
            y[ifreq] = abs(sol.dft.calc(Ntime, VPoint[ipoint, 0: Ntime], freq, dt, 0))
            #print(ifreq, y[ifreq])
        # 最大値で正規化
        y /= max(y)

        # plot
        ax.plot(x, y)
        ax.grid(True)

        # X軸
        ax.set_xlabel('frequency %s' % Post['funit'])
        ax.set_xlim(x[0], x[-1])
        ax.set_xticks(np.linspace(x[0], x[-1], Post['freqdiv'] + 1))

        # Y軸
        ax.set_ylim([0, 1])
        ax.set_ylabel('V(f)')

        # タイトル
        strtitle = "%s\nspectrum point #%d" % (Post['title'], ipoint + 1)
        ax.set_title(strtitle)
