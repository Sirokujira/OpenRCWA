# -*- coding: utf-8 -*-
"""
material.py
"""

import math
import numpy as np
from numba import jit
import sol.geometry

# Yee格子の電界点と磁界点の物性値番号を計算する
@jit(cache=True, nopython=True)
def setup(
    Nx, Ny, Nz, Xn, Yn, Zn, Xc, Yc, Zc, iGeometry, fGeometry,
    iEx, iEy, iEz, iHx, iHy, iHz,
    iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0):

    # 微小量
    eps = 1e-6 * math.sqrt( \
        (Xn[0] - Xn[-1])**2 + \
        (Yn[0] - Yn[-1])**2 + \
        (Zn[0] - Zn[-1])**2)
    #print(eps)

    for n in range(iGeometry.shape[0]):
        m = iGeometry[n, 0]
        s = iGeometry[n, 1]
        p = fGeometry[n, :]
        #print(s, m, p)

        # bounding box
        x1, x2, y1, y2, z1, z2 = sol.geometry.boundingbox(s, p)
        #print(x1, x2, y1, y2, z1, z2)

        i1 = i2 = j1 = j2 = k1 = k2 = 0

        # Ex
        if s == 1:
            i1, i2 = _getSpan(Xc, Nx,     iMin, iMax - 1, x1, x2, eps)
            j1, j2 = _getSpan(Yn, Ny + 1, jMin, jMax,     y1, y2, eps)
            k1, k2 = _getSpan(Zn, Nz + 1, kMin, kMax,     z1, z2, eps)
        else:
            i1 = iMin
            i2 = iMax - 1
            j1 = jMin
            j2 = jMax
            k1 = kMin
            k2 = kMax

        for i in range(i1, i2 + 1):
            for j in range(j1, j2 + 1):
                for k in range(k1, k2 + 1):
                    if sol.geometry.inside(Xc[i], Yn[j], Zn[k], s, p, eps):
                        n = (Ni * i) + (Nj * j) + (Nk * k) + N0
                        iEx[n] = m

        # Ey
        if s == 1:
            j1, j2 = _getSpan(Yc, Ny,     jMin, jMax - 1, y1, y2, eps)
            k1, k2 = _getSpan(Zn, Nz + 1, kMin, kMax,     z1, z2, eps)
            i1, i2 = _getSpan(Xn, Nx + 1, iMin, iMax,     x1, x2, eps)
        else:
            j1 = jMin
            j2 = jMax - 1
            k1 = kMin
            k2 = kMax
            i1 = iMin
            i2 = iMax

        for i in range(i1, i2 + 1):
            for j in range(j1, j2 + 1):
                for k in range(k1, k2 + 1):
                    if sol.geometry.inside(Xn[i], Yc[j], Zn[k], s, p, eps):
                        n = (Ni * i) + (Nj * j) + (Nk * k) + N0
                        iEy[n] = m

        # Ez
        if s == 1:
            k1, k2 = _getSpan(Zc, Nz,     kMin, kMax - 1, z1, z2, eps)
            i1, i2 = _getSpan(Xn, Nx + 1, iMin, iMax,     x1, x2, eps)
            j1, j2 = _getSpan(Yn, Ny + 1, jMin, jMax,     y1, y2, eps)
        else:
            k1 = kMin
            k2 = kMax - 1
            i1 = iMin
            i2 = iMax
            j1 = jMin
            j2 = jMax

        for i in range(i1, i2 + 1):
            for j in range(j1, j2 + 1):
                for k in range(k1, k2 + 1):
                    if sol.geometry.inside(Xn[i], Yn[j], Zc[k], s, p, eps):
                        n = (Ni * i) + (Nj * j) + (Nk * k) + N0
                        iEz[n] = m

        # Hx
        if s == 1:
            i1, i2 = _getSpan(Xn, Nx + 1, iMin,     iMax, x1, x2, eps)
            j1, j2 = _getSpan(Yc, Ny,     jMin - 1, jMax, y1, y2, eps)
            k1, k2 = _getSpan(Zc, Nz,     kMin - 1, kMax, z1, z2, eps)
        else:
            i1 =     iMin
            i2 =     iMax
            j1 = max(jMin - 1,      0)
            j2 = min(jMax,     Ny - 1)
            k1 = max(kMin - 1,      0)
            k2 = min(kMax,     Nz - 1)

        for i in range(i1, i2 + 1):
            for j in range(j1, j2 + 1):
                for k in range(k1, k2 + 1):
                    if sol.geometry.inside(Xn[i], Yc[j], Zc[k], s, p, eps):
                        n = (Ni * i) + (Nj * j) + (Nk * k) + N0
                        iHx[n] = m

        # Hy
        if s == 1:
            j1, j2 = _getSpan(Yn, Ny + 1, jMin,     jMax, y1, y2, eps)
            k1, k2 = _getSpan(Zc, Nz,     kMin - 1, kMax, z1, z2, eps)
            i1, i2 = _getSpan(Xc, Nx,     iMin - 1, iMax, x1, x2, eps)
        else:
            j1 =     jMin
            j2 =     jMax
            k1 = max(kMin - 1,      0)
            k2 = min(kMax,     Nz - 1)
            i1 = max(iMin - 1,      0)
            i2 = min(iMax,     Nx - 1)

        for i in range(i1, i2 + 1):
            for j in range(j1, j2 + 1):
                for k in range(k1, k2 + 1):
                    if sol.geometry.inside(Xc[i], Yn[j], Zc[k], s, p, eps):
                        n = (Ni * i) + (Nj * j) + (Nk * k) + N0
                        iHy[n] = m

        # Hz
        if s == 1:
            k1, k2 = _getSpan(Zn, Nz + 1, kMin,     kMax, z1, z2, eps)
            i1, i2 = _getSpan(Xc, Nx,     iMin - 1, iMax, x1, x2, eps)
            j1, j2 = _getSpan(Yc, Ny,     jMin - 1, jMax, y1, y2, eps)
        else:
            k1 =     kMin
            k2 =     kMax
            i1 = max(iMin - 1,      0)
            i2 = min(iMax,     Nx - 1)
            j1 = max(jMin - 1,      0)
            j2 = min(jMax,     Ny - 1)

        for i in range(i1, i2 + 1):
            for j in range(j1, j2 + 1):
                for k in range(k1, k2 + 1):
                    if sol.geometry.inside(Xc[i], Yc[j], Zn[k], s, p, eps):
                        n = (Ni * i) + (Nj * j) + (Nk * k) + N0
                        iHz[n] = m

    return iEx, iEy, iEz, iHx, iHy, iHz

# vector用物性値配列を作成する
def alloc_vector(Parm, NN):

    f_dtype = Parm['f_dtype']

    K1Ex = np.zeros(NN, f_dtype)
    K2Ex = np.zeros(NN, f_dtype)
    K1Ey = np.zeros(NN, f_dtype)
    K2Ey = np.zeros(NN, f_dtype)
    K1Ez = np.zeros(NN, f_dtype)
    K2Ez = np.zeros(NN, f_dtype)
    K1Hx = np.zeros(NN, f_dtype)
    K2Hx = np.zeros(NN, f_dtype)
    K1Hy = np.zeros(NN, f_dtype)
    K2Hy = np.zeros(NN, f_dtype)
    K1Hz = np.zeros(NN, f_dtype)
    K2Hz = np.zeros(NN, f_dtype)

    return \
        K1Ex, K2Ex, K1Ey, K2Ey, K1Ez, K2Ez, \
        K1Hx, K2Hx, K1Hy, K2Hy, K1Hz, K2Hz

# vector用物性値配列を作成する
@jit(cache=True, nopython=True)
def setup_vector(
    iEx, iEy, iEz, iHx, iHy, iHz,
    C1E, C2E, C1H, C2H,
    K1Ex, K2Ex, K1Ey, K2Ey, K1Ez, K2Ez,
    K1Hx, K2Hx, K1Hy, K2Hy, K1Hz, K2Hz,
    iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0):

    # Ex
    for i in range(iMin, iMax + 0):
        for j in range(jMin, jMax + 1):
            for k in range(kMin, kMax + 1):
                n = (Ni * i) + (Nj * j) + (Nk * k) + N0
                m = iEx[n]
                K1Ex[n] = C1E[m]
                K2Ex[n] = C2E[m]

    # Ey
    for i in range(iMin, iMax + 1):
        for j in range(jMin, jMax + 0):
            for k in range(kMin, kMax + 1):
                n = (Ni * i) + (Nj * j) + (Nk * k) + N0
                m = iEy[n]
                K1Ey[n] = C1E[m]
                K2Ey[n] = C2E[m]

    # Ez
    for i in range(iMin, iMax + 1):
        for j in range(jMin, jMax + 1):
            for k in range(kMin, kMax + 0):
                n = (Ni * i) + (Nj * j) + (Nk * k) + N0
                m = iEz[n]
                K1Ez[n] = C1E[m]
                K2Ez[n] = C2E[m]

    # Hx
    for i in range(iMin, iMax + 1):
        for j in range(jMin, jMax + 0):
            for k in range(kMin, kMax + 0):
                n = (Ni * i) + (Nj * j) + (Nk * k) + N0
                m = iHx[n]
                K1Hx[n] = C1H[m]
                K2Hx[n] = C2H[m]

    # Hy
    for i in range(iMin, iMax + 0):
        for j in range(jMin, jMax + 1):
            for k in range(kMin, kMax + 0):
                n = (Ni * i) + (Nj * j) + (Nk * k) + N0
                m = iHy[n]
                K1Hy[n] = C1H[m]
                K2Hy[n] = C2H[m]

    # Hz
    for i in range(iMin, iMax + 0):
        for j in range(jMin, jMax + 0):
            for k in range(kMin, kMax + 1):
                n = (Ni * i) + (Nj * j) + (Nk * k) + N0
                m = iHz[n]
                K1Hz[n] = C1H[m]
                K2Hz[n] = C2H[m]

    return \
        K1Ex, K2Ex, K1Ey, K2Ey, K1Ez, K2Ez, \
        K1Hx, K2Hx, K1Hy, K2Hy, K1Hz, K2Hz

# correct surface index
# correct curved surface (E <- H)
@jit(cache=True)
def correct_surface(
    Nx, Ny, Nz, iMaterial,
    iEx, iEy, iEz, iHx, iHy, iHz,
    iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0):

    # Ex
    for i in range(iMin, iMax + 0):
        for j in range(jMin, jMax + 1):
            for k in range(kMin, kMax + 1):
                if (j > 0) and(j < Ny) and (k > 0) and (k < Nz):
                    n   = Ni * (i    ) + Nj * (j    ) + Nk * (k    ) + N0
                    nj1 = Ni * (i    ) + Nj * (j - 1) + Nk * (k    ) + N0
                    nk1 = Ni * (i    ) + Nj * (j    ) + Nk * (k - 1) + N0
                    iEx[n] = _highest(iMaterial, \
                        iEx[n], iHy[n], iHy[nk1], iHz[n], iHz[nj1])

    # Ey
    for i in range(iMin, iMax + 1):
        for j in range(jMin, jMax + 0):
            for k in range(kMin, kMax + 1):
                if (k > 0) and (k < Nz) and (i > 0) and (i < Nx):
                    n   = Ni * (i    ) + Nj * (j    ) + Nk * (k    ) + N0
                    nk1 = Ni * (i    ) + Nj * (j    ) + Nk * (k - 1) + N0
                    ni1 = Ni * (i - 1) + Nj * (j    ) + Nk * (k    ) + N0
                    iEy[n] = _highest(iMaterial, \
                        iEy[n], iHz[n], iHz[ni1], iHx[n], iHx[nk1])

    # Ez
    for i in range(iMin, iMax + 1):
        for j in range(jMin, jMax + 1):
            for k in range(kMin, kMax + 0):
                if (i > 0) and (i < Nx) and (j > 0) and (j < Ny):
                    n   = Ni * (i    ) + Nj * (j    ) + Nk * (k    ) + N0
                    ni1 = Ni * (i - 1) + Nj * (j    ) + Nk * (k    ) + N0
                    nj1 = Ni * (i    ) + Nj * (j - 1) + Nk * (k    ) + N0
                    iEz[n] = _highest(iMaterial, \
                        iEz[n], iHx[n], iHx[nj1], iHy[n], iHy[ni1])

# 物性値係数を計算する
def factor(Parm, iMaterial, fMaterial):

    f_dtype = Parm['f_dtype']
    nmaterial = len(iMaterial)

    C1E = np.zeros(nmaterial, f_dtype)
    C2E = np.zeros(nmaterial, f_dtype)
    C1H = np.zeros(nmaterial, f_dtype)
    C2H = np.zeros(nmaterial, f_dtype)

    # air (1,1)
    C1E[0] = 1; C2E[0] = 1
    C1H[0] = 1; C2H[0] = 1

    # PEC (0,0)
    C1E[1] = 0; C2E[1] = 0
    C1H[1] = 0; C2H[1] = 0

    for m in range(2, nmaterial):
        # 誘電率係数
        if   iMaterial[m] == 1:
            # 通常媒質
            epsr = fMaterial[m, 0]
            esgm = fMaterial[m, 1]
            edenom = epsr + (esgm * Parm['ETA0'] * Parm['C'] * Parm['dt'])
            C1E[m] = epsr / edenom
            C2E[m] = 1 / edenom
        elif iMaterial[m] == 2:
            # 分散性媒質
            einf = fMaterial[m, 4]
            ae   = fMaterial[m, 5]
            be   = fMaterial[m, 6]
            ce   = fMaterial[m, 7]
            ke = math.exp(-ce * Parm['dt'])
            xi0 = (ae * Parm['dt']) + (be / ce) * (1 - ke)
            edenom = einf + xi0
            C1E[m] = einf / edenom
            C2E[m] = 1    / edenom

        # 透磁率係数
        amur = fMaterial[m, 2]
        msgm = fMaterial[m, 3]
        mdenom = amur + (msgm / Parm['ETA0'] * Parm['C'] * Parm['dt'])
        C1H[m] = amur / mdenom
        C2H[m] = 1 / mdenom
    #print(C1E, C2E)
    #print(C1H, C2H)

    return C1E, C2E, C1H, C2H

# (private) highest material
# type = 1 or 2
@jit(cache=True)
def _highest(iMaterial, id0, id1, id2, id3, id4):
    PEC = 1#Parm['PEC']

    ret = id0

    if   (id1 == id0) and (id2 == id0) and (id3 == id0) and (id4 == id0):
        pass
    elif (id0 == PEC) or (id1 == PEC) or (id2 == PEC) or (id3 == PEC) or (id4 == PEC):
        ret = PEC
    elif iMaterial[id1] == 2:
        ret = id1
    elif iMaterial[id2] == 2:
        ret = id2
    elif iMaterial[id3] == 2:
        ret = id3
    elif iMaterial[id4] == 2:
        ret = id4

    return ret

# p[i1] <= p[n1] <= p1 <= p2 <= p[n2] <= p[i2]
# output : n1, n2
@jit(cache=True)
def _getSpan(p, n, i1, i2, p1, p2, eps):
    i1 = max(0, min(i1, n - 1))
    i2 = max(0, min(i2, n - 1))

    # k1 <= k2
    k1 = min(i1, i2)
    k2 = max(i1, i2)

    # q1 <= q2
    q1 = min(p1, p2)
    q2 = max(p1, p2)

    n1 = k1
    n2 = k2

    if n1 == n2:
        return n1, n2

    if   (q1 < p[k1] - eps) and (q2 < p[k1] - eps):
        # p1, p2 < p[k1] -> n1 > n2
        n1 = k1 - 1
        n2 = k1 - 2
    elif (q1 > p[k2] + eps) and (q2 > p[k2] + eps):
        # p1, p2 > p[k2] -> n1 > n2
        n1 = k2 + 2
        n2 = k2 + 1
    else:
        # p[n1] <= p1
        for k in range(k1, k2):
            if (q1 > p[k] - eps) and (q1 < p[k + 1] - eps):
                n1 = k
                #printf("A %d" % n1)
                break
        if abs(q1 - p[k2]) < eps:
            n1 = k2
            #print("A2 %d" % n1)
        # p2 <= p[n2]
        for k in range(k2, k1, -1):
            if (q2 < p[k] + eps) and (q2 > p[k - 1] + eps):
                n2 = k
                #print("B %d" % n2)
                break
        if abs(q2 - p[k1]) < eps:
            n2 = k1
            #printf("B2 %d" % n2)

    #assert((i1 <= n1) and (n1 <= n2) and (n2 <= i2))

    return n1, n2
"""
# debug
if __name__ == "__main__":
    p = np.arange(8).astype('f8')
    n1, n2 = _getSpan(p, len(p), 0, 6, 1.0, 6.1, 1e-6)
    print(n1, n2)
"""