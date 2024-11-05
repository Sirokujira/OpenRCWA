# -*- coding: utf-8 -*-
"""
geometry.py
物体形状関係の処理
OpenRCWA/OpenSTF共通
"""

import math
import numpy as np
from numba import jit

# 物体形状の3D線分データを取得する
def lines(iGeometry, fGeometry):

    eps = 1e-6

    # 円の分割数
    rdiv = 36

    # 線分数を取得する
    nline = _nlines(iGeometry, rdiv)
    #print(nline)

    # 線分配列作成
    gline = np.empty((nline, 2, 3), 'f8')
    mline = np.empty(nline, 'i4')

    # 線分データに代入する
    nline = 0
    for n in range(iGeometry.shape[0]):
        m = iGeometry[n, 0]
        s = iGeometry[n, 1]
        p = fGeometry[n, :]

        x1, x2, y1, y2, z1, z2 = p[0:6]
        x0, y0, z0 = np.array([x1 + x2, y1 + y2, z1 + z2]) / 2

        if   s == 1:
            # 直方体
            dx = abs(x2 - x1)
            dy = abs(y2 - y1)
            dz = abs(z2 - z1)

            # 線分
            if   (dy < eps) and (dz < eps):
                nline = _line(nline, gline, mline, m, x1, y1, z1, x2, y1, z1)
            elif (dz < eps) and (dx < eps):
                nline = _line(nline, gline, mline, m, x1, y1, z1, x1, y2, z1)
            elif (dx < eps) and (dy < eps):
                nline = _line(nline, gline, mline, m, x1, y1, z1, x1, y1, z2)
            # 平面
            elif (dx < eps):
                nline = _rectangle(nline, gline, mline, m, 'X', x1, y1, z1, y2, z2)
            elif (dy < eps):
                nline = _rectangle(nline, gline, mline, m, 'Y', y1, z1, x1, z2, x2)
            elif (dz < eps):
                nline = _rectangle(nline, gline, mline, m, 'Z', z1, x1, y1, x2, y2)
            # 直方体
            else:
                nline = _rectangle(nline, gline, mline, m, 'Z', z1, x1, y1, x2, y2)
                nline = _rectangle(nline, gline, mline, m, 'Z', z2, x1, y1, x2, y2)
                nline = _line(nline, gline, mline, m, x1, y1, z1, x1, y1, z2)
                nline = _line(nline, gline, mline, m, x2, y1, z1, x2, y1, z2)
                nline = _line(nline, gline, mline, m, x2, y2, z1, x2, y2, z2)
                nline = _line(nline, gline, mline, m, x1, y2, z1, x1, y2, z2)
        elif s == 2:
            # 楕円体
            nline = _ellipse(nline, gline, mline, m, 'X', x0, y1, z1, y2, z2, rdiv)
            nline = _ellipse(nline, gline, mline, m, 'Y', y0, z1, x1, z2, x2, rdiv)
            nline = _ellipse(nline, gline, mline, m, 'Z', z0, x1, y1, x2, y2, rdiv)
        elif s == 11:
            # X方向楕円柱
            nline = _ellipse(nline, gline, mline, m, 'X', x1, y1, z1, y2, z2, rdiv)
            nline = _ellipse(nline, gline, mline, m, 'X', x2, y1, z1, y2, z2, rdiv)
            nline = _line(nline, gline, mline, m, x1, y1, z0, x2, y1, z0)
            nline = _line(nline, gline, mline, m, x1, y2, z0, x2, y2, z0)
            nline = _line(nline, gline, mline, m, x1, y0, z1, x2, y0, z1)
            nline = _line(nline, gline, mline, m, x1, y0, z2, x2, y0, z2)
        elif s == 12:
            # Y方向楕円柱
            nline = _ellipse(nline, gline, mline, m, 'Y', y1, z1, x1, z2, x2, rdiv)
            nline = _ellipse(nline, gline, mline, m, 'Y', y2, z1, x1, z2, x2, rdiv)
            nline = _line(nline, gline, mline, m, x0, y1, z1, x0, y2, z1)
            nline = _line(nline, gline, mline, m, x0, y1, z2, x0, y2, z2)
            nline = _line(nline, gline, mline, m, x1, y1, z0, x1, y2, z0)
            nline = _line(nline, gline, mline, m, x2, y1, z0, x2, y2, z0)
        elif s == 13:
            # Z方向楕円柱
            nline = _ellipse(nline, gline, mline, m, 'Z', z1, x1, y1, x2, y2, rdiv)
            nline = _ellipse(nline, gline, mline, m, 'Z', z2, x1, y1, x2, y2, rdiv)
            nline = _line(nline, gline, mline, m, x1, y0, z1, x1, y0, z2)
            nline = _line(nline, gline, mline, m, x2, y0, z1, x2, y0, z2)
            nline = _line(nline, gline, mline, m, x0, y1, z1, x0, y1, z2)
            nline = _line(nline, gline, mline, m, x0, y2, z1, x0, y2, z2)
        elif s == 31:
            # X方向三角柱
            px = p[0:2]
            py = p[2:5]
            pz = p[5:8]
            for i in range(2):
                nline = _line(nline, gline, mline, m, px[i], py[0], pz[0], px[i], py[1], pz[1])
                nline = _line(nline, gline, mline, m, px[i], py[1], pz[1], px[i], py[2], pz[2])
                nline = _line(nline, gline, mline, m, px[i], py[2], pz[2], px[i], py[0], pz[0])
            for i in range(3):
                nline = _line(nline, gline, mline, m, px[0], py[i], pz[i], px[1], py[i], pz[i])
        elif s == 32:
            # Y方向三角柱
            py = p[0:2]
            pz = p[2:5]
            px = p[5:8]
            for i in range(2):
                nline = _line(nline, gline, mline, m, px[0], py[i], pz[0], px[1], py[i], pz[1])
                nline = _line(nline, gline, mline, m, px[1], py[i], pz[1], px[2], py[i], pz[2])
                nline = _line(nline, gline, mline, m, px[2], py[i], pz[2], px[0], py[i], pz[0])
            for i in range(3):
                nline = _line(nline, gline, mline, m, px[i], py[0], pz[i], px[i], py[1], pz[i])
        elif s == 33:
            # Z方向三角柱
            pz = p[0:2]
            px = p[2:5]
            py = p[5:8]
            for i in range(2):
                nline = _line(nline, gline, mline, m, px[0], py[0], pz[i], px[1], py[1], pz[i])
                nline = _line(nline, gline, mline, m, px[1], py[1], pz[i], px[2], py[2], pz[i])
                nline = _line(nline, gline, mline, m, px[2], py[2], pz[i], px[0], py[0], pz[i])
            for i in range(3):
                nline = _line(nline, gline, mline, m, px[i], py[i], pz[0], px[i], py[i], pz[1])
        elif s == 41:
            # X方向角錐台
            x1, x2  = p[0:2]
            y0, z0  = p[2:4]
            h1y, h1z = p[4:6] / 2
            h2y, h2z = p[6:8] / 2
            nline = _rectangle(nline, gline, mline, m, 'X', x1, y0 - h1y, z0 - h1z, y0 + h1y, z0 + h1z)
            nline = _rectangle(nline, gline, mline, m, 'X', x2, y0 - h2y, z0 - h2z, y0 + h2y, z0 + h2z)
            nline = _line(nline, gline, mline, m, x1, y0 - h1y, z0 - h1z, x2, y0 - h2y, z0 - h2z)
            nline = _line(nline, gline, mline, m, x1, y0 - h1y, z0 + h1z, x2, y0 - h2y, z0 + h2z)
            nline = _line(nline, gline, mline, m, x1, y0 + h1y, z0 - h1z, x2, y0 + h2y, z0 - h2z)
            nline = _line(nline, gline, mline, m, x1, y0 + h1y, z0 + h1z, x2, y0 + h2y, z0 + h2z)
        elif s == 42:
            # Y方向角錐台
            y1, y2  = p[0:2]
            z0, x0  = p[2:4]
            h1z, h1x = p[4:6] / 2
            h2z, h2x = p[6:8] / 2
            nline = _rectangle(nline, gline, mline, m, 'Y', y1, z0 - h1z, x0 - h1x, z0 + h1z, x0 + h1x)
            nline = _rectangle(nline, gline, mline, m, 'Y', y2, z0 - h2z, x0 - h2x, z0 + h2z, x0 + h2x)
            nline = _line(nline, gline, mline, m, x0 - h1x, y1, z0 - h1z, x0 - h2x, y2, z0 - h2z)
            nline = _line(nline, gline, mline, m, x0 - h1x, y1, z0 + h1z, x0 - h2x, y2, z0 + h2z)
            nline = _line(nline, gline, mline, m, x0 + h1x, y1, z0 - h1z, x0 + h2x, y2, z0 - h2z)
            nline = _line(nline, gline, mline, m, x0 + h1x, y1, z0 + h1z, x0 + h2x, y2, z0 + h2z)
        elif s == 43:
                # Z方向角錐台
            z1, z2  = p[0:2]
            x0, y0  = p[2:4]
            h1x, h1y = p[4:6] / 2
            h2x, h2y = p[6:8] / 2
            nline = _rectangle(nline, gline, mline, m, 'Z', z1, x0 - h1x, y0 - h1y, x0 + h1x, y0 + h1y)
            nline = _rectangle(nline, gline, mline, m, 'Z', z2, x0 - h2x, y0 - h2y, x0 + h2x, y0 + h2y)
            nline = _line(nline, gline, mline, m, x0 - h1x, y0 - h1y, z1, x0 - h2x, y0 - h2y, z2)
            nline = _line(nline, gline, mline, m, x0 - h1x, y0 + h1y, z1, x0 - h2x, y0 + h2y, z2)
            nline = _line(nline, gline, mline, m, x0 + h1x, y0 - h1y, z1, x0 + h2x, y0 - h2y, z2)
            nline = _line(nline, gline, mline, m, x0 + h1x, y0 + h1y, z1, x0 + h2x, y0 + h2y, z2)
        elif s == 51:
            # X方向円錐台
            x1, x2  = p[0:2]
            y0, z0  = p[2:4]
            r1y, r1z = p[4:6] / 2
            r2y, r2z = p[6:8] / 2
            nline = _ellipse(nline, gline, mline, m, 'X', x1, y0 - r1y, z0 - r1z, y0 + r1y, z0 + r1z, rdiv)
            nline = _ellipse(nline, gline, mline, m, 'X', x2, y0 - r2y, z0 - r2z, y0 + r2y, z0 + r2z, rdiv)
            nline = _line(nline, gline, mline, m, x1, y0,       z0 + r1z, x2, y0,       z0 + r2z)
            nline = _line(nline, gline, mline, m, x1, y0,       z0 - r1z, x2, y0,       z0 - r2z)
            nline = _line(nline, gline, mline, m, x1, y0 + r1y, z0,       x2, y0 + r2y, z0      )
            nline = _line(nline, gline, mline, m, x1, y0 - r1y, z0,       x2, y0 - r2y, z0      )
        elif s == 52:
            # Y方向円錐台
            y1, y2  = p[0:2]
            z0, x0  = p[2:4]
            r1z, r1x = p[4:6] / 2
            r2z, r2x = p[6:8] / 2
            nline = _ellipse(nline, gline, mline, m, 'Y', y1, z0 - r1z, x0 - r1x, z0 + r1z, x0 + r1x, rdiv)
            nline = _ellipse(nline, gline, mline, m, 'Y', y2, z0 - r2z, x0 - r2x, z0 + r2z, x0 + r2x, rdiv)
            nline = _line(nline, gline, mline, m, x0,       y1, z0 + r1z, x0,       y2, z0 + r2z)
            nline = _line(nline, gline, mline, m, x0,       y1, z0 - r1z, x0,       y2, z0 - r2z)
            nline = _line(nline, gline, mline, m, x0 + r1x, y1, z0,       x0 + r2x, y2, z0      )
            nline = _line(nline, gline, mline, m, x0 - r1x, y1, z0,       x0 - r2x, y2, z0      )
        elif s == 53:
            # Z方向円錐台
            z1, z2  = p[0:2]
            x0, y0  = p[2:4]
            r1x, r1y = p[4:6] / 2
            r2x, r2y = p[6:8] / 2
            nline = _ellipse(nline, gline, mline, m, 'Z', z1, x0 - r1x, y0 - r1y, x0 + r1x, y0 + r1y, rdiv)
            nline = _ellipse(nline, gline, mline, m, 'Z', z2, x0 - r2x, y0 - r2y, x0 + r2x, y0 + r2y, rdiv)
            nline = _line(nline, gline, mline, m, x0,       y0 + r1y, z1, x0,       y0 + r2y, z2)
            nline = _line(nline, gline, mline, m, x0,       y0 - r1y, z1, x0,       y0 - r2y, z2)
            nline = _line(nline, gline, mline, m, x0 + r1x, y0,       z1, x0 + r2x, y0,       z2)
            nline = _line(nline, gline, mline, m, x0 - r1x, y0,       z1, x0 - r2x, y0,       z2)

    # 配列の大きさを正味に修正する
    #print(nline)
    gline = np.resize(gline, (nline, 2, 3))
    mline = np.resize(mline, nline)

    return gline, mline

# (private) 線分データ追加
def _line(nline, gline, mline, m, x1, y1, z1, x2, y2, z2):

    gline[nline] = [[x1, y1, z1], [x2, y2, z2]]
    mline[nline] = m

    return nline + 1

# (private) 長方形データ追加
def _rectangle(nline, gline, mline, m, cdir, c0, p1, q1, p2, q2):

    if   cdir == 'X':
        nline = _line(nline, gline, mline, m, c0, p1, q1, c0, p2, q1)
        nline = _line(nline, gline, mline, m, c0, p2, q1, c0, p2, q2)
        nline = _line(nline, gline, mline, m, c0, p2, q2, c0, p1, q2)
        nline = _line(nline, gline, mline, m, c0, p1, q2, c0, p1, q1)
    elif cdir == 'Y':
        nline = _line(nline, gline, mline, m, q1, c0, p1, q1, c0, p2)
        nline = _line(nline, gline, mline, m, q1, c0, p2, q2, c0, p2)
        nline = _line(nline, gline, mline, m, q2, c0, p2, q2, c0, p1)
        nline = _line(nline, gline, mline, m, q2, c0, p1, q1, c0, p1)
    elif cdir == 'Z':
        nline = _line(nline, gline, mline, m, p1, q1, c0, p2, q1, c0)
        nline = _line(nline, gline, mline, m, p2, q1, c0, p2, q2, c0)
        nline = _line(nline, gline, mline, m, p2, q2, c0, p1, q2, c0)
        nline = _line(nline, gline, mline, m, p1, q2, c0, p1, q1, c0)

    return nline

# (private) 楕円データ追加
def _ellipse(nline, gline, mline,m, cdir, c0, p1, q1, p2, q2, rdiv):

    p0 = (p1 + p2) / 2
    q0 = (q1 + q2) / 2
    pr = abs(p1 - p2) / 2
    qr = abs(q1 - q2) / 2

    for i in range(rdiv):
        a1 = (i + 0) * (2 * math.pi) / rdiv
        a2 = (i + 1) * (2 * math.pi) / rdiv
        x1 = y1 = z1 = x2 = y2 = z2 = 0

        if   cdir == 'X':
            x1 = c0
            y1 = p0 + pr * math.cos(a1)
            z1 = q0 + qr * math.sin(a1)
            x2 = c0
            y2 = p0 + pr * math.cos(a2)
            z2 = q0 + qr * math.sin(a2)
        elif cdir == 'Y':
            y1 = c0
            z1 = p0 + pr * math.cos(a1)
            x1 = q0 + qr * math.sin(a1)
            y2 = c0
            z2 = p0 + pr * math.cos(a2)
            x2 = q0 + qr * math.sin(a2)
        elif cdir == 'Z':
            z1 = c0
            x1 = p0 + pr * math.cos(a1)
            y1 = q0 + qr * math.sin(a1)
            z2 = c0
            x2 = p0 + pr * math.cos(a2)
            y2 = q0 + qr * math.sin(a2)

        nline = _line(nline, gline, mline, m, x1, y1, z1, x2, y2, z2)

    return nline

# (private) 線分数を取得する
# 注意 : 厚さなしのときは多めに取っている
def _nlines(iGeometry, rdiv):

    nline = 0

    for n in range(iGeometry.shape[0]):
        s = iGeometry[n, 1]
        if   s in [1, 41, 42, 43]:
            nline += 12
        elif s in [2]:
            nline += 3 * rdiv
        elif s in [11, 12, 13, 51, 52, 53]:
            nline += (2 * rdiv) + 4
        elif s in [31, 32, 33]:
            nline += 9

    return nline

# 物体形状を囲む範囲を求める
@jit(cache=True, nopython=True)
def boundingbox(s, p):

    x1 = x2 = y1 = y2 = z1 = z2 = 0

    # cube / sphere / cylinder
    if (s == 1) or (s == 2) or \
       (s == 11) or (s == 12) or (s == 13):
        x1 = min(p[0], p[1])
        x2 = max(p[0], p[1])
        y1 = min(p[2], p[3])
        y2 = max(p[2], p[3])
        z1 = min(p[4], p[5])
        z2 = max(p[4], p[5])
    # triangle pillar
    elif (s == 31):
        x1 = min(p[0], p[1])
        x2 = max(p[0], p[1])
        y1 = min(p[2], p[3], p[4])
        y2 = max(p[2], p[3], p[4])
        z1 = min(p[5], p[6], p[7])
        z2 = max(p[5], p[6], p[7])
    elif (s == 32):
        y1 = min(p[0], p[1])
        y2 = max(p[0], p[1])
        z1 = min(p[2], p[3], p[4])
        z2 = max(p[2], p[3], p[4])
        x1 = min(p[5], p[6], p[7])
        x2 = max(p[5], p[6], p[7])
    elif (s == 33):
        z1 = min(p[0], p[1])
        z2 = max(p[0], p[1])
        x1 = min(p[2], p[3], p[4])
        x2 = max(p[2], p[3], p[4])
        y1 = min(p[5], p[6], p[7])
        y2 = max(p[5], p[6], p[7])
    # pyramid / cone
    elif (s == 41) or (s == 51):
        hy = max(p[4], p[6]) / 2
        hz = max(p[5], p[7]) / 2
        x1 = min(p[0], p[1])
        x2 = max(p[0], p[1])
        y1 = p[2] - hy
        y2 = p[2] + hy
        z1 = p[3] - hz
        z2 = p[3] + hz
    elif (s == 42) or (s == 52):
        hz = max(p[4], p[6]) / 2
        hx = max(p[5], p[7]) / 2
        y1 = min(p[0], p[1])
        y2 = max(p[0], p[1])
        z1 = p[2] - hz
        z2 = p[2] + hz
        x1 = p[3] - hx
        x2 = p[3] + hx
    elif (s == 43) or (s == 53):
        hx = max(p[4], p[6]) / 2
        hy = max(p[5], p[7]) / 2
        z1 = min(p[0], p[1])
        z2 = max(p[0], p[1])
        x1 = p[2] - hx
        x2 = p[2] + hx
        y1 = p[3] - hy
        y2 = p[3] + hy

    return x1, x2, y1, y2, z1, z2

# 点(x,y,z)は物体形状(s,p)の内部または境界上にあるか
# eps : 計算対象の大きさより十分距離
@jit(cache=True, nopython=True)
def inside(x, y, z, s, p, eps):

    zero = 1e-6
    eps2 = eps**2

    # rectangle
    if   s == 1:
        if ((x - p[0]) * (x - p[1]) <= eps2) and \
           ((y - p[2]) * (y - p[3]) <= eps2) and \
           ((z - p[4]) * (z - p[5]) <= eps2):
            return True
    # ellipsoid
    elif s == 2:
        x0 = (p[0] + p[1]) / 2
        y0 = (p[2] + p[3]) / 2
        z0 = (p[4] + p[5]) / 2
        xr = abs(p[0] - p[1]) / 2
        yr = abs(p[2] - p[3]) / 2
        zr = abs(p[4] - p[5]) / 2
        if (x - x0) * (x - x0) / (xr * xr) \
         + (y - y0) * (y - y0) / (yr * yr) \
         + (z - z0) * (z - z0) / (zr * zr) < 1 + zero:
            return True
    # cylinder
    elif s == 11:
        # X cylinder
        y0 = (p[2] + p[3]) / 2
        z0 = (p[4] + p[5]) / 2
        yr = abs(p[2] - p[3]) / 2
        zr = abs(p[4] - p[5]) / 2
        if ((x - p[0]) * (x - p[1]) <= eps2) and \
            ((y - y0) * (y - y0) / (yr * yr) \
           + (z - z0) * (z - z0) / (zr * zr) < 1 + zero):
            return True
    elif s == 12:
        # Y clinder
        x0 = (p[0] + p[1]) / 2
        z0 = (p[4] + p[5]) / 2
        xr = abs(p[0] - p[1]) / 2
        zr = abs(p[4] - p[5]) / 2
        if ((y - p[2]) * (y - p[3]) <= eps2) and \
            ((z - z0) * (z - z0) / (zr * zr) \
           + (x - x0) * (x - x0) / (xr * xr) < 1 + zero):
            return True
    elif s == 13:
        # Z cylinder
        x0 = (p[0] + p[1]) / 2
        y0 = (p[2] + p[3]) / 2
        xr = abs(p[0] - p[1]) / 2
        yr = abs(p[2] - p[3]) / 2
        if ((z - p[4]) * (z - p[5]) <= eps2) and \
            ((x - x0) * (x - x0) / (xr * xr) \
           + (y - y0) * (y - y0) / (yr * yr) < 1 + zero):
            return True
    # pillar
    elif s == 31:
        # X-pillar
        if (x - p[0]) * (x - p[1]) > eps2:
            return False
        else:
            return _inout3(y, z, p[2:5], p[5:8], zero)
    elif s == 32:
        # Y-pillar
        if (y - p[0]) * (y - p[1]) > eps2:
            return False
        else:
            return _inout3(z, x, p[2:5], p[5:8], zero)
    elif s == 33:
        # Z-pillar
        if (z - p[0]) * (z - p[1]) > eps2:
            return False
        else:
            return _inout3(x, y, p[2:5], p[5:8], zero)
    # pyramid
    elif s == 41:
        # X-pyramid
        x1, x2  = p[0:2]
        y0, z0  = p[2:4]
        h1y, h1z = p[4:6] / 2
        h2y, h2z = p[6:8] / 2
        f = (x - x1) / (x2 - x1) if (abs(x1 - x2) > eps) else 0
        hy = h1y + f * (h2y - h1y)
        hz = h1z + f * (h2z - h1z)
        if (x - x1) * (x - x2) < eps2 and \
            abs(y - y0) < hy + eps and \
            abs(z - z0) < hz + eps:
            return True
    elif s == 42:
        # Y-pyramid
        y1, y2  = p[0:2]
        z0, x0  = p[2:4]
        h1z, h1x = p[4:6] / 2
        h2z, h2x = p[6:8] / 2
        f = (y - y1) / (y2 - y1) if (abs(y1 - y2) > eps) else 0
        hz = h1z + f * (h2z - h1z)
        hx = h1x + f * (h2x - h1x)
        if (y - y1) * (y - y2) < eps2 and \
            abs(z - z0) < hz + eps and \
            abs(x - x0) < hx + eps:
            return True
    elif s == 43:
        # Z-pyramid
        z1, z2  = p[0:2]
        x0, y0  = p[2:4]
        h1x, h1y = p[4:6] / 2
        h2x, h2y = p[6:8] / 2
        f = (z - z1) / (z2 - z1) if (abs(z1 - z2) > eps) else 0
        hx = h1x + f * (h2x - h1x)
        hy = h1y + f * (h2y - h1y)
        if (z - z1) * (z - z2) < eps2 and \
            abs(x - x0) < hx + eps and \
            abs(y - y0) < hy + eps:
            return True

    # cone
    elif s == 51:
        # X-cone
        x1, x2  = p[0:2]
        y0, z0  = p[2:4]
        r1y, r1z = p[4:6] / 2
        r2y, r2z = p[6:8] / 2
        f = (x - x1) / (x2 - x1) if (abs(x1 - x2) > eps) else 0
        ry = r1y + f * (r2y - r1y)
        rz = r1z + f * (r2z - r1z)
        if (x - x1) * (x - x2) < eps2 and \
            ((y - y0) / ry)**2 + ((z - z0) / rz)**2 < 1 + zero:
            return True
    elif s == 52:
        # Y-cone
        y1, y2  = p[0:2]
        z0, x0  = p[2:4]
        r1z, r1x = p[4:6] / 2
        r2z, r2x = p[6:8] / 2
        f = (y - y1) / (y2 - y1) if (abs(y1 - y2) > eps) else 0
        rz = r1z + f * (r2z - r1z)
        rx = r1x + f * (r2x - r1x)
        if (y - y1) * (y - y2) < eps2 and \
            ((z - z0) / rz)**2 + ((x - x0) / rx)**2 < 1 + zero:
            return True
    elif s == 53:
        # Z-cone
        z1, z2  = p[0:2]
        x0, y0  = p[2:4]
        r1x, r1y = p[4:6] / 2
        r2x, r2y = p[6:8] / 2
        f = (z - z1) / (z2 - z1) if (abs(z1 - z2) > eps) else 0
        rx = r1x + f * (r2x - r1x)
        ry = r1y + f * (r2y - r1y)
        if (z - z1) * (z - z2) < eps2 and \
            ((x - x0) / rx)**2 + ((y - y0) / ry)**2 < 1 + zero:
            return True

    return False

# (private) 点(x,y)は三角形の内部か境界上にあるか
# zero : 無次元の微小量
# tx[3], ty[3] : 三角形のXY座標
@jit(cache=True, nopython=True)
def _inout3(x, y, tx, ty, zero):

    eps = zero * ((max(tx) - min(tx)) + (max(ty) - min(ty)))

    x1, x2, x3 = tx
    y1, y2, y3 = ty

    det = ((x2 - x1) * (y3 - y1)) - ((x3 - x1) * (y2 - y1))

    if abs(det) > eps**2:
        a = + ((x - x1) * (y3 - y1) - (y - y1) * (x3 - x1)) / det
        b = - ((x - x1) * (y2 - y1) - (y - y1) * (x2 - x1)) / det
        if (a > -zero) and (b > -zero) and (a + b < 1 + zero):
            return True

    return False
