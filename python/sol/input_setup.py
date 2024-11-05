# -*- coding: utf-8 -*-
"""
input_setup.py
"""

import math
import numpy as np

# セル中心
def cell_center(
    Xn, Yn, Zn):

    Xc = (Xn[0:-1] + Xn[1:]) / 2
    Yc = (Yn[0:-1] + Yn[1:]) / 2
    Zc = (Zn[0:-1] + Zn[1:]) / 2
    #print(Xc, Yc, Zc)

    return Xc, Yc, Zc

# 平面波入射
def planewave(Parm, Xn, Yn, Zn, Freq2, fPlanewave):
    # 入力データ
    theta = Parm['planewave'][0]
    phi   = Parm['planewave'][1]
    pol   = Parm['planewave'][2]

    # sin, cos
    cost = math.cos(math.radians(theta))
    sint = math.sin(math.radians(theta))
    cosp = math.cos(math.radians(phi))
    sinp = math.sin(math.radians(phi))

    # (r1, t1, p1) : (r, theta, phi)方向の単位ベクトル
    r1 = np.array([+ sint * cosp, + sint * sinp, + cost])
    t1 = np.array([+ cost * cosp, + cost * sinp, - sint])
    p1 = np.array([- sinp,        + cosp,          0   ])

    # ri : 進行方向単位ベクトル
    ri = -r1

    # ei : 電界単位ベクトル
    if   pol == 1:
        # V-pol
        ei = -t1
    elif pol == 2:
        # H-pol
        ei = +p1

    # hi = ei X r : 磁界単位ベクトル
    hi = np.cross(ei, r1)

    # 初期位置
    f0 = (Freq2[0] + Freq2[-1]) / 2
    r = math.sqrt((Xn[0] - Xn[-1])**2 + \
                  (Yn[0] - Yn[-1])**2 + \
                  (Zn[0] - Zn[-1])**2) / 2 + (0.5 * Parm['C'] / f0)
    r0 = np.array([Xn[0] + Xn[-1], \
                   Yn[0] + Yn[-1], \
                   Zn[0] + Zn[-1]]) / 2 - (r * ri)

    # waveform parameter
    ai = 4 / (1.27 / f0)

    # 変数に代入
    #fPlanewave = np.zeros(15, 'f8')
    fPlanewave[ 0:  3] = ei
    fPlanewave[ 3:  6] = hi
    fPlanewave[ 6:  9] = ri
    fPlanewave[ 9: 12] = r0
    fPlanewave[12]     = ai
    fPlanewave[13]     = Parm['dt']; assert(fPlanewave[13] > 0)  # コピー
    fPlanewave[14]     = Parm['C'];  assert(fPlanewave[14] > 0)  # コピー

    #return fPlanewave

# 入力データに必要なデータを追加する
def setup(
    Parm, Xn, Yn, Zn, Xc, Yc, Zc,
    iMaterial, fMaterial, iGeometry, fGeometry, iFeed, fFeed, 
    iPoint, fPoint, iInductor, fInductor, Freq1, Freq2):

    # タイムステップ (既定値)
    if Parm['dt'] < 1e-20:
        Parm['dt'] = _timestep(Parm, Xn, Yn, Zn)

    # パルス幅 (既定値)
    if Parm['tw'] < 1e-20:
        Parm['tw'] = _pulsewidth(Parm, Freq1, Freq2)

    # 給電点
    if Parm['source'] == 0:
        _feed(iFeed, fFeed, Xn, Yn, Zn, Xc, Yc, Zc)
    #print(iFeed, fFeed)

    # 観測点
    if iPoint.shape[0] > 0:
        _point(Parm, Xn, Yn, Zn, Xc, Yc, Zc, iPoint, fPoint)

    # R/C
    _load(Parm, Xn, Yn, Zn, Xc, Yc, Zc, iMaterial, fMaterial, iGeometry, fGeometry)

    # inductor
    if iInductor.shape[0] > 0:
        _inductor(Parm, Xn, Yn, Zn, Xc, Yc, Zc, iInductor, fInductor)

    # 厚さのない板を節点に寄せる
    if iGeometry.shape[0] > 0:
        _fit_geometry(Parm, Xn, Yn, Zn, iGeometry, fGeometry)

# === 以下はprivate関数 ===

# タイムステップ 
def _timestep(Parm, Xn, Yn, Zn):
    dxmin = np.min(Xn[1:] - Xn[:-1])
    dymin = np.min(Yn[1:] - Yn[:-1])
    dzmin = np.min(Zn[1:] - Zn[:-1])
    dt = 1 / math.sqrt(1 / dxmin**2 + 1 / dymin**2 + 1 / dzmin**2) / Parm['C']
    #print(dxmin, dymin, dzmin, dt)

    return dt

# パルス幅
def _pulsewidth(Parm, Freq1, Freq2):
    f0 = 0
    if   len(Freq1) > 0:
        f0 = (Freq1[0] + Freq1[-1]) / 2
    elif len(Freq2) > 0:
        f0 = (Freq2[0] + Freq2[-1]) / 2
    else:
        f0 = 1 / (20 * Parm['dt'])

    tw = 1.27 / f0
    #print("f0=%e Tw=%e" % (f0, tw))

    return tw

# 給電点
def _feed(iFeed, fFeed, Xn, Yn, Zn, Xc, Yc, Zc):
    for n in range(iFeed.shape[0]):
        # (i, j, k) と (dx, dy, dz)
        iFeed[n, 1:4], fFeed[n, 6:9] = _index_length(iFeed[n, 0], fFeed[n, 0:3], Xn, Yn, Zn, Xc, Yc, Zc)
    #print(iFeed, fFeed)

# 観測点
def _point(Parm, Xn, Yn, Zn, Xc, Yc, Zc, iPoint, fPoint):
    # 1+, 1- 観測点を加える
    npoint = iPoint.shape[0] - 2
    iPoint[npoint + 0, :] = iPoint[0, :]
    fPoint[npoint + 0, :] = fPoint[0, :]
    iPoint[npoint + 1, :] = iPoint[0, :]
    fPoint[npoint + 1, :] = fPoint[0, :]

    # (i, j, k), (dx, dy, dz) を計算する
    for n in range(npoint + 2):
        iPoint[n, 1:4], fPoint[n, 3:6] = _index_length(iPoint[n, 0], fPoint[n, 0:3], Xn, Yn, Zn, Xc, Yc, Zc)

    # prop = 観測点1の伝搬方向 = +X/-X/+Y/-Y/+Z/-Z = 0...5
    prop = Parm['prop']
    idir = prop // 2  # = X/Y/Z = 0/1/2
    sign = +1 if (prop % 2) == 0 else -1  # = +1/-1

    # 1+/1- のindexを調整する
    idx = 1 + idir  # = X/Y/Z = i/j/k = 1/2/3
    iPoint[npoint + 0, idx] = iPoint[0, idx] + sign
    iPoint[npoint + 1, idx] = iPoint[0, idx] - sign
    #print(iPoint, fPoint)

# R/C
def _load(Parm, Xn, Yn, Zn, Xc, Yc, Zc, iMaterial, fMaterial, iGeometry, fGeometry):

    for n in range(iGeometry.shape[0]):
        shape = iGeometry[n, 1]
        # RとCのみ
        if ((90 <= shape) and (shape <= 95)):
            # 物体形状
            rorc = 0 if shape <= 92 else 1  # 0/1=R/C
            idir = shape % 3  # 0/1/2=X/Y/Z方向
            gpos = fGeometry[n, 0:3]
            load = fGeometry[n, 3]
            [i, j, k], [dx, dy, dz] = _index_length(idir, gpos, Xn, Yn, Zn, Xc, Yc, Zc)
            lpos = None
            if   idir == 0:
                lpos = [Xn[i], Xn[i + 1], Yn[j], Yn[j], Zn[k], Zn[k]]
                dlds = dx / (dy * dz)
            elif idir == 1:
                lpos = [Xn[i], Xn[i], Yn[j], Yn[j + 1], Zn[k], Zn[k]]
                dlds = dy / (dz * dx)
            elif idir == 2:
                lpos = [Xn[i], Xn[i], Yn[j], Yn[j], Zn[k], Zn[k + 1]]
                dlds = dz / (dx * dy)
            iGeometry[n, 1] = 1
            fGeometry[n, 0:6] = lpos
            # 物性値
            nmaterial = iGeometry[n, 0]  # 物性値番号
            iMaterial[nmaterial] = 1  # 通常媒質
            fMaterial[nmaterial, 0] = 1 if (rorc == 0) else (load / Parm['EPS0']) * dlds
            fMaterial[nmaterial, 1] = (1 / load) * dlds if (rorc == 0) else 0
            fMaterial[nmaterial, 2] = 1
            fMaterial[nmaterial, 3] = 0
            #print(g)
            #print(nmaterial)
            #print(iMaterial[nmaterial])
            #print(fMaterial[nmaterial, :])

# inductor (L)
def _inductor(Parm, Xn, Yn, Zn, Xc, Yc, Zc, iInductor, fInductor):

    for n in range(iInductor.shape[0]):
        idir = iInductor[n, 0]  # = X/Y/Z = 0/1/2
        # (i, j, k) と (dx, dy, dz)
        iInductor[n, 1:4], fInductor[n, 3:6] = _index_length(idir, fInductor[n, 0:3], Xn, Yn, Zn, Xc, Yc, Zc)

        # その他の変数
        dx, dy, dz = fInductor[n, 3:6]
        dlds = [dx / (dy * dz), dy / (dz * dx), dz / (dx * dy)][idir]
        fInductor[n, 7] = Parm['MU0'] * dlds / fInductor[n, 6]
        fInductor[n, 8:10] = 0
    #print(iInductor, fInductor)

# 厚さのない物体形状を節点に寄せる
def _fit_geometry(Parm, Xn, Yn, Zn, iGeometry, fGeometry):
    # 微小量
    d0 = Parm['EPS'] * ( \
        abs(Xn[0] - Xn[-1]) + \
        abs(Yn[0] - Yn[-1]) + \
        abs(Zn[0] - Zn[-1]))

    for n in range(iGeometry.shape[0]):
        shape = iGeometry[n, 1]
        p = fGeometry[n, :]
        x1, x2, y1, y2, z1, z2 = p[0:6]

        if   (shape == 1) or (shape == 2):
            # rectangle, ellipsoid
            if abs(x1 - x2) < d0:
                i = np.argmin(abs(Xn - x1))
                p[0:2] = Xn[i]
            if abs(y1 - y2) < d0:
                j = np.argmin(abs(Yn - y1))
                p[2:4] = Yn[j]
            if abs(z1 - z2) < d0:
                k = np.argmin(abs(Zn - z1))
                p[4:6] = Zn[k]
        elif shape == 11:
            # X cylinder
            if abs(x1 - x2) < d0:
                i = np.argmin(abs(Xn - x1))
                p[0:2] = Xn[i]
        elif shape == 12:
            # Y cylinder
            if abs(y1 - y2) < d0:
                j = np.argmin(abs(Yn - y1))
                p[2:4] = Yn[j]
        elif shape == 13:
            # Z cylinder
            if abs(z1 - z2) < d0:
                k = np.argmin(abs(Zn - z1))
                p[4:6] = Zn[k]

# (private) indexとメッシュ長
# 指定した向き(idir=0/1/2)と場所(pos[0:3])に一番近いYee格子点の(i,j,k)と(dx,dy.dz)を求める (private)
def _index_length(idir, pos, Xn, Yn, Zn, Xc, Yc, Zc):
    x = pos[0]
    y = pos[1]
    z = pos[2]

    Nx = len(Xc)
    Ny = len(Yc)
    Nz = len(Zc)

    if   idir == 0:
        # X direction
        i = np.argmin(abs(Xc - x))
        j = np.argmin(abs(Yn - y))
        k = np.argmin(abs(Zn - z))
        j = max(1, min(Ny - 1, j))
        k = max(1, min(Nz - 1, k))
        dx = Xn[i + 1] - Xn[i]
        dy = (Yc[j] - Yc[j - 1]) if (Ny > 1) else (Yn[1] - Yn[0])
        dz = (Zc[k] - Zc[k - 1]) if (Nz > 1) else (Zn[1] - Zn[0])
    elif idir == 1:
        # Y direction
        j = np.argmin(abs(Yc - y))
        k = np.argmin(abs(Zn - z))
        i = np.argmin(abs(Xn - x))
        k = max(1, min(Nz - 1, k))
        i = max(1, min(Nx - 1, i))
        dy = Yn[j + 1] - Yn[j]
        dz = (Zc[k] - Zc[k - 1]) if (Nz > 1) else (Zn[1] - Zn[0])
        dx = (Xc[i] - Xc[i - 1]) if (Nx > 1) else (Xn[1] - Xn[0])
    elif idir == 2:
        # Z direction
        k = np.argmin(abs(Zc - z))
        i = np.argmin(abs(Xn - x))
        j = np.argmin(abs(Yn - y))
        i = max(1, min(Nx - 1, i))
        j = max(1, min(Ny - 1, j))
        dz = Zn[k + 1] - Zn[k]
        dx = (Xc[i] - Xc[i - 1]) if (Nx > 1) else (Xn[1] - Xn[0])
        dy = (Yc[j] - Yc[j - 1]) if (Ny > 1) else (Yn[1] - Yn[0])

    return [i, j, k], [dx, dy, dz]
