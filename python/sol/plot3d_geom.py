# -*- coding: utf-8 -*-
"""
plot3d_geom.py
入力データを3D図形表示する
"""

import datetime
import matplotlib.pyplot as plt

# (private) メッシュその他を描く
# 戻り値 : ax
#@jit(cache=True)
def _plot(Parm, Nx, Ny, Nz, Xn, Yn, Zn):

    # figure
    strfig = 'OpenRCWA - geometry (3D) - %s' % datetime.datetime.now().ctime()
    fig = plt.figure(strfig, figsize=(5, 5))  # 5 inches
    ax = fig.add_subplot(projection='3d')

    # X-const mesh
    for i in range(Nx + 1):
        x = [Xn[i]] * 3
        y = [Yn[0], Yn[0], Yn[-1]]
        z = [Zn[-1], Zn[0], Zn[0]]
        ax.plot(x, y, z, color='lightgray')

    # Y-const mesh
    for j in range(Ny + 1):
        y = [Yn[j]] * 3
        z = [Zn[0], Zn[0], Zn[-1]]
        x = [Xn[-1], Xn[0], Xn[0]]
        ax.plot(x, y, z, color='lightgray')

    # Z-const mesh
    for k in range(Nz + 1):
        z = [Zn[k]] * 3
        x = [Xn[0], Xn[0], Xn[-1]]
        y = [Yn[-1], Yn[0], Yn[0]]
        ax.plot(x, y, z, color='lightgray')

    # layout
    ax.grid(False)
    #ax.axis(False)
    ax.set_aspect('equal')
    ax.view_init(elev = 30, azim = 30, roll = 0)

    # axis
    ax.set_xlim(Xn[0], Xn[-1])
    ax.set_ylim(Yn[0], Yn[-1])
    ax.set_zlim(Zn[0], Zn[-1])
    ax.set_xticks([Xn[0], Xn[-1]])
    ax.set_yticks([Yn[0], Yn[-1]])
    ax.set_zticks([Zn[0], Zn[-1]])
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    # title
    ax.set_title('%s\ncells = %d x %d x %d' % (Parm['title'], Nx, Ny, Nz))

    return ax

# 物体形状の外枠を描く
def shape(Parm, Nx, Ny, Nz, Xn, Yn, Zn, iFeed, gline, mline):

    ax = _plot(Parm, Nx, Ny, Nz, Xn, Yn, Zn)

    #print(gline)
    #print(mline)
    # 物体形状
    for n in range(len(mline)):
        color = 'k' if mline[n] == 1 else 'm'
        ax.plot([gline[n, 0, 0], gline[n, 1, 0]], \
                [gline[n, 0, 1], gline[n, 1, 1]], \
                [gline[n, 0, 2], gline[n, 1, 2]], color=color)

    # 給電点
    for n in range(iFeed.shape[0]):
        idir, i, j, k = iFeed[n, 0:4]
        if   idir == 0:
            ax.plot([Xn[i], Xn[i + 1]], [Yn[j]] * 2, [Zn[k]] * 2, 'r')
        elif idir == 1:
            ax.plot([Xn[i]] * 2, [Yn[j], Yn[j + 1]], [Zn[k]] * 2, 'r')
        elif idir == 2:
            ax.plot([Xn[i]] * 2, [Yn[j]] * 2, [Zn[k], Zn[k + 1]], 'r')

# Yee格子の電界点を描く
#@jit(cache=True)
def cell(Parm, Nx, Ny, Nz, Xn, Yn, Zn, iEx, iEy, iEz, Ni, Nj, Nk, N0):

    ax = _plot(Parm, Nx, Ny, Nz, Xn, Yn, Zn)

    # Ex
    for i in range(Nx + 0):
        for j in range(Ny + 1):
            for k in range(Nz + 1):
                n = (Ni * i) + (Nj * j) + (Nk * k) + N0
                m = iEx[n]
                if m > 0:
                    color = ['k', 'm', 'c', 'y'][min(m - 1, 3)]
                    ax.plot([Xn[i], Xn[i + 1]], [Yn[j]]*2, [Zn[k]]*2, color=color)

    # Ey
    for i in range(Nx + 1):
        for j in range(Ny + 0):
            for k in range(Nz + 1):
                n = (Ni * i) + (Nj * j) + (Nk * k) + N0
                m = iEy[n]
                if m > 0:
                    color = ['k', 'm', 'c', 'y'][min(m - 1, 3)]
                    ax.plot([Xn[i]]*2, [Yn[j], Yn[j + 1]], [Zn[k]]*2, color=color)

    # Ez
    for i in range(Nx + 1):
        for j in range(Ny + 1):
            for k in range(Nz + 0):
                n = (Ni * i) + (Nj * j) + (Nk * k) + N0
                m = iEz[n]
                if m > 0:
                    color = ['k', 'm', 'c', 'y'][min(m - 1, 3)]
                    ax.plot([Xn[i]]*2, [Yn[j]]*2, [Zn[k], Zn[k + 1]], color=color)
