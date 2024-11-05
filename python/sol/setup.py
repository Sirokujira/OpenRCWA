# -*- coding: utf-8 -*-
"""
setup.py
各種準備作業
"""

import math, sys
import numpy as np
from numba import jit, prange
import sol.material, sol.dft
import sol.setupMurHx, sol.setupMurHy, sol.setupMurHz
import sol.setupPmlEx, sol.setupPmlEy, sol.setupPmlEz, sol.setupPmlHx, sol.setupPmlHy, sol.setupPmlHz
import sol.setupDispEx, sol.setupDispEy, sol.setupDispEz

def setup(
    Parm, fPlanewave, Nx, Ny, Nz, Xn, Yn, Zn, Xc, Yc, Zc,
    iMaterial, fMaterial, iGeometry, fGeometry, Freq2,
    iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0, NN):

    # 物性値係数を計算する
    C1E, C2E, C1H, C2H \
    = sol.material.factor(Parm, iMaterial, fMaterial)

    # 物性値番号の配列を作成する
    i_dtype = Parm['i_dtype']
    iEx = np.zeros(NN, i_dtype)
    iEy = np.zeros(NN, i_dtype)
    iEz = np.zeros(NN, i_dtype)
    iHx = np.zeros(NN, i_dtype)
    iHy = np.zeros(NN, i_dtype)
    iHz = np.zeros(NN, i_dtype)

    # 物性値番号を計算する
    sol.material.setup(
        Nx, Ny, Nz, Xn, Yn, Zn, Xc, Yc, Zc, iGeometry, fGeometry,
        iEx, iEy, iEz, iHx, iHy, iHz,
        iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)

    #print(np.sum(iEx), np.sum(iEy), np.sum(iEz), np.sum(iHx), np.sum(iHy), np.sum(iHz))
    # 表面の物性値番号を調整する(精度向上のため)
    sol.material.correct_surface(
        Nx, Ny, Nz, iMaterial,
        iEx, iEy, iEz, iHx, iHy, iHz,
        iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)
    #print(np.sum(iEx), np.sum(iEy), np.sum(iEz), np.sum(iHx), np.sum(iHy), np.sum(iHz))

    # vectorモードの物性値係数
    K1Ex = K2Ex = K1Ey = K2Ey = K1Ez = K2Ez = None
    K1Hx = K2Hx = K1Hy = K2Hy = K1Hz = K2Hz = None
    if Parm['vector']:
        K1Ex, K2Ex, K1Ey, K2Ey, K1Ez, K2Ez, \
        K1Hx, K2Hx, K1Hy, K2Hy, K1Hz, K2Hz \
        = sol.material.alloc_vector(Parm, NN)
        sol.material.setup_vector(
            iEx, iEy, iEz, iHx, iHy, iHz,
            C1E, C2E, C1H, C2H,
            K1Ex, K2Ex, K1Ey, K2Ey, K1Ez, K2Ez,
            K1Hx, K2Hx, K1Hy, K2Hy, K1Hz, K2Hz,
            iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)

    # 分散性媒質の準備
    _dispersion(iMaterial, fMaterial, Parm['dt'])

    # 分散性媒質の係数を計算する
    iDispEx, fDispEx \
    = sol.setupDispEx.setEx(
        iMaterial, fMaterial, iEx,
        iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)
    iDispEy, fDispEy \
    = sol.setupDispEy.setEy(
        iMaterial, fMaterial, iEy,
        iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)
    iDispEz, fDispEz \
    = sol.setupDispEz.setEz(
        iMaterial, fMaterial, iEz,
        iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)

    # メッシュ因子を計算する
    RXn, RYn, RZn, RXc, RYc, RZc \
    = _meshfactor(Parm, Nx, Ny, Nz, Xn, Yn, Zn)
    #print(RXc)

    # ABC
    fMurHx = fMurHy = fMurHz = iMurHx = iMurHy = iMurHz = None
    iPmlEx = iPmlEy = iPmlEz = iPmlHx = iPmlHy = iPmlHz = None
    gPmlXn = gPmlYn = gPmlZn = gPmlXc = gPmlYc = gPmlZc = None
    rPmlE = rPmlH = rPml = None
    if   Parm['abc'][0] == 0:
        # Mur
        fMurHx, fMurHy, fMurHz, iMurHx, iMurHy, iMurHz, \
        = _mur(Parm, Nx, Ny, Nz, Xn, Yn, Zn, fMaterial, iEx, iEy, iEz,
            iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)
    elif Parm['abc'][0] == 1:
        # PML
        iPmlEx, iPmlEy, iPmlEz, iPmlHx, iPmlHy, iPmlHz \
        = _pml(
            Parm, Nx, Ny, Nz, iEx, iEy, iEz, iHx, iHy, iHz,
            iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)
        gPmlXn, gPmlYn, gPmlZn, gPmlXc, gPmlYc, gPmlZc, rPmlE, rPmlH, rPml \
        = _pml_factor(
            Parm, Nx, Ny, Nz, Xn, Yn, Zn, Xc, Yc, Zc, fMaterial)

    # DFT用の係数を準備する
    cEdft, cHdft, cFdft \
    = sol.dft.setup(Parm, fPlanewave, Xn, Yn, Zn, Freq2)

    return \
    iEx, iEy, iEz, iHx, iHy, iHz, \
    C1E, C2E, C1H, C2H, \
    K1Ex, K2Ex, K1Ey, K2Ey, K1Ez, K2Ez, \
    K1Hx, K2Hx, K1Hy, K2Hy, K1Hz, K2Hz, \
    RXn, RYn, RZn, RXc, RYc, RZc, \
    fMurHx, fMurHy, fMurHz, iMurHx, iMurHy, iMurHz, \
    iPmlEx, iPmlEy, iPmlEz, iPmlHx, iPmlHy, iPmlHz, \
    gPmlXn, gPmlYn, gPmlZn, gPmlXc, gPmlYc, gPmlZc, rPmlE, rPmlH, rPml, \
    iDispEx, iDispEy, iDispEz, fDispEx, fDispEy, fDispEz, \
    cEdft, cHdft, cFdft

# 電磁界配列のindex
# 領域分割しないときは Npx=Npy=Npz=npx=npy=npz=1, comm_rank=0
def getIndex(Parm, Nx, Ny, Nz, Npx, Npy, Npz, npx, npy, npz, comm_rank):

    # too many process (MPI)
    if (npx > Nx) or (npy > Ny) or (npz > Nz):
        if comm_rank == 0:
            print("*** too many process = %dx%dx%d (limit = %dx%dx%d)" % (npx, npy, npz, Nx, Ny, Nz))
        sys.exit()

    # ABC関係のパラメーター
    abc = Parm['abc'][0]  # 0(Mur) / 1(PML)
    pml = Parm['abc'][1]  # PMLの層数
    lx = 1 if (abc == 0) else pml
    ly = 1 if (abc == 0) else pml
    lz = 1 if (abc == 0) else pml

    # 計算領域のX/Y/Z方向のindexの下限と上限
    # 冗長分を含む、正味は電磁界成分で異なる
    # MPIを想定して別途変数(iMin,...)を用意する
    # Ex[iMin - 0: iMax + 0, jMin - 0: jMax + 1, kMin - 0: kMax + 1]
    # Ey[iMin - 0: iMax + 1, jMin - 0: jMax + 0, kMin - 0: kMax + 1]
    # Ez[iMin - 0: iMax + 1, jMin - 0: jMax + 1, kMin - 0: kMax + 0]
    # Hx[iMin - 1, iMax + 2, jMin - 1, jMax + 1, kMin - 1, kMax + 1]
    # Hy[iMin - 1, iMax + 1, jMin - 1, jMax + 2, kMin - 1, kMax + 1]
    # Hz[iMin - 1, iMax + 1, jMin - 1, jMax + 1, kMin - 1, kMax + 2]

    # MPI : 領域番号(Ipx, Ipy, Ipz)を取得する
    Ipx = Ipy = Ipz = 0
    ip = 0
    for i in range(npx):
        for j in range(npy):
            for k in range(npz):
                if comm_rank == ip:
                    Ipx = i
                    Ipy = j
                    Ipz = k
                ip += 1
    

    # min, max
    iMin, iMax = _idminmax(Nx, npx, Npx, Ipx)
    jMin, jMax = _idminmax(Ny, npy, Npy, Ipy)
    kMin, kMax = _idminmax(Nz, npz, Npz, Ipz)
    #print(comm_rank, lx, ly, lz)
    #print(comm_rank, iMin, iMax, jMin, jMax, kMin, kMax)

    # 電磁界の1次元配列計算用の係数
    # 配列の番号 = Ni*i + Nj*k + Nk*k + N0
    Nk = 1
    Nj = (kMax - kMin + (2 * lz) + 1)
    Ni = (jMax - jMin + (2 * ly) + 1) * Nj
    N0 = -((iMin - lx) * Ni + (jMin - ly) * Nj + (kMin - lz) * Nk)

    # 電磁界配列の最大数(成分で共通)
    NN = (iMax + lx) * Ni + (jMax + ly) * Nj + (kMax + lz) * Nk + N0 + 1

    #print(comm_rank, Ni, Nj, Nk, N0, NN)
    #print(comm_rank, (Ni * (iMin - 1)) + (Nj * (jMin - 1)) + (Nk * (kMin - 1)) + N0)
    #print(comm_rank, (Ni * (iMax + 1)) + (Nj * (jMax + 1)) + (Nk * (kMax + 1)) + N0 + 1)
    assert((Ni * (iMin - lx)) + (Nj * (jMin - ly)) + (Nk * (kMin - lz)) + N0 == 0)
    assert((Ni * (iMax + lx)) + (Nj * (jMax + ly)) + (Nk * (kMax + lz)) + N0 == NN - 1)

    return iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0, NN, Ipx, Ipy, Ipz

# === 以下はprivate関数 ===

# 領域indexの上下限を求める
def _idminmax(n, _np, np, ip):

    idmin = 0
    idmax = n

    if _np > 1:
        # MPI
        nc = max(n // np, 1)
        idmin = (ip + 0) * nc
        idmax = (ip + 1) * nc
        if ip == np - 1:
            idmax = n

    return idmin, idmax

# Mur ABC
def _mur(Parm, Nx, Ny, Nz, Xn, Yn, Zn, fMaterial, iEx, iEy, iEz,
    iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0):

    f_dtype = Parm['f_dtype']
    cdt = Parm['C'] * Parm['dt']

    # 配列の型を宣言する(jitに必要)
    fMurHx = np.zeros((0, 2), f_dtype)
    fMurHy = np.zeros((0, 2), f_dtype)
    fMurHz = np.zeros((0, 2), f_dtype)
    iMurHx = np.zeros((0, 6), 'i4')
    iMurHy = np.zeros((0, 6), 'i4')
    iMurHz = np.zeros((0, 6), 'i4')

    # 配列の大きさを取得する
    num = [0] * 3
    num[0] = sol.setupMurHx.setHx(0,
        Ny, Nz, Yn, Zn, fMaterial, iEy, iEz, fMurHx, iMurHx, cdt,
        iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)
    num[1] = sol.setupMurHy.setHy(0,
        Nx, Nz, Xn, Zn, fMaterial, iEz, iEx, fMurHy, iMurHy, cdt,
        iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)
    num[2] = sol.setupMurHz.setHz(0,
        Nx, Ny, Xn, Yn, fMaterial, iEx, iEy, fMurHz, iMurHz, cdt,
        iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)
    #print(num)

    # 配列の大きさを宣言する
    fMurHx = np.zeros((num[0], 2), f_dtype)
    fMurHy = np.zeros((num[1], 2), f_dtype)
    fMurHz = np.zeros((num[2], 2), f_dtype)
    iMurHx = np.zeros((num[0], 6), 'i4')
    iMurHy = np.zeros((num[1], 6), 'i4')
    iMurHz = np.zeros((num[2], 6), 'i4')

    # 配列にデータを代入する
    sol.setupMurHx.setHx(1,
        Ny, Nz, Yn, Zn, fMaterial, iEy, iEz, fMurHx, iMurHx, cdt,
        iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)
    sol.setupMurHy.setHy(1,
        Nx, Nz, Xn, Zn, fMaterial, iEz, iEx, fMurHy, iMurHy, cdt,
        iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)
    sol.setupMurHz.setHz(1,
        Nx, Ny, Xn, Yn, fMaterial, iEx, iEy, fMurHz, iMurHz, cdt,
        iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)

    #print(num, np.sum(fMurHx), np.sum(fMurHy), np.sum(fMurHz))
    return fMurHx, fMurHy, fMurHz, iMurHx, iMurHy, iMurHz

# PML ABC
def _pml(Parm, Nx, Ny, Nz, iEx, iEy, iEz, iHx, iHy, iHz,
    iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0):

    # PML層数L
    pml_l = Parm['abc'][1]
    
    # jitを使う場合はダミーで初期化する必要がある(None引数は不可)
    iPmlEx = np.zeros((0, 4), 'i4')
    iPmlEy = np.zeros((0, 4), 'i4')
    iPmlEz = np.zeros((0, 4), 'i4')
    iPmlHx = np.zeros((0, 4), 'i4')
    iPmlHy = np.zeros((0, 4), 'i4')
    iPmlHz = np.zeros((0, 4), 'i4')

    # 配列の大きさを取得する
    num = [0] * 6
    num[0] = sol.setupPmlEx.setEx(0, pml_l, Nx, Ny, Nz, iEx, iPmlEx,
        iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)
    num[1] = sol.setupPmlEy.setEy(0, pml_l, Nx, Ny, Nz, iEy, iPmlEy,
        iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)
    num[2] = sol.setupPmlEz.setEz(0, pml_l, Nx, Ny, Nz, iEz, iPmlEz,
        iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)
    num[3] = sol.setupPmlHx.setHx(0, pml_l, Nx, Ny, Nz, iHx, iPmlHx,
        iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)
    num[4] = sol.setupPmlHy.setHy(0, pml_l, Nx, Ny, Nz, iHy, iPmlHy,
        iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)
    num[5] = sol.setupPmlHz.setHz(0, pml_l, Nx, Ny, Nz, iHz, iPmlHz,
        iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)

    # 配列を宣言する(i, j, k, m)
    iPmlEx = np.zeros((num[0], 4), 'i4')
    iPmlEy = np.zeros((num[1], 4), 'i4')
    iPmlEz = np.zeros((num[2], 4), 'i4')
    iPmlHx = np.zeros((num[3], 4), 'i4')
    iPmlHy = np.zeros((num[4], 4), 'i4')
    iPmlHz = np.zeros((num[5], 4), 'i4')

    # 配列にデータを代入する
    sol.setupPmlEx.setEx(1, pml_l, Nx, Ny, Nz, iEx, iPmlEx,
        iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)
    sol.setupPmlEy.setEy(1, pml_l, Nx, Ny, Nz, iEy, iPmlEy,
        iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)
    sol.setupPmlEz.setEz(1, pml_l, Nx, Ny, Nz, iEz, iPmlEz,
        iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)
    sol.setupPmlHx.setHx(1, pml_l, Nx, Ny, Nz, iHx, iPmlHx,
        iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)
    sol.setupPmlHy.setHy(1, pml_l, Nx, Ny, Nz, iHy, iPmlHy,
        iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)
    sol.setupPmlHz.setHz(1, pml_l, Nx, Ny, Nz, iHz, iPmlHz,
        iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)

    return iPmlEx, iPmlEy, iPmlEz, iPmlHx, iPmlHy, iPmlHz

# PML因子
def _pml_factor(Parm, Nx, Ny, Nz, Xn, Yn, Zn, Xc, Yc, Zc, fMaterial):
    #fc = d = 0
    f_dtype = Parm['f_dtype']

    l  = Parm['abc'][1]
    m  = Parm['abc'][2]
    r0 = Parm['abc'][3]
    kpml = (m + 1) / (2.0 * l) * math.log(1.0 / r0)

    cdt = Parm['C'] * Parm['dt']

    # 作業配列
    f = np.zeros(2 * l + 1, 'f8')
    f[0] = f[1] = 0
    for n in range(2, 2 * l + 1):
        f[n] = kpml * math.pow((n - 1) / (2.0 * l), m)

    # E
    gPmlXn = np.zeros(Nx + 2 * l, f_dtype)
    gPmlYn = np.zeros(Ny + 2 * l, f_dtype)
    gPmlZn = np.zeros(Nz + 2 * l, f_dtype)
    for i in range(-l + 1, Nx + l):
        if   i <= 0:
            fc = f[- 2 * (i     )]; d = Xn[1] - Xn[0]
        elif i >= Nx:
            fc = f[+ 2 * (i - Nx)]; d = Xn[-1] - Xn[-2]
        else:
            fc = 0;                 d = Xc[i] - Xc[i - 1]
        gPmlXn[i + l] = cdt / d * fc
    for j in range(-l + 1, Ny + l):
        if   j <= 0:
            fc = f[- 2 * (j     )]; d = Yn[1] - Yn[0]
        elif j >= Ny:
            fc = f[+ 2 * (j - Ny)]; d = Yn[-1] - Yn[-2]
        else:
            fc = 0;                 d = Yc[j] - Yc[j - 1]
        gPmlYn[j + l] = cdt / d * fc
    for k in range(-l + 1, Nz + l):
        if   k <= 0:
            fc = f[- 2 * (k     )]; d = Zn[1] - Zn[0]
        elif k >= Nz:
            fc = f[+ 2 * (k - Nz)]; d = Zn[-1] - Zn[-2]
        else:
            fc = 0;                 d = Zc[k] - Zc[k - 1]
        gPmlZn[k + l] = cdt / d * fc

    # H
    gPmlXc = np.zeros(Nx + 2 * l, f_dtype)
    gPmlYc = np.zeros(Ny + 2 * l, f_dtype)
    gPmlZc = np.zeros(Nz + 2 * l, f_dtype)
    for i in range(-l, Nx + l):
        if   i <  0:
            fc = f[- 2 * (i     ) - 1]; d = Xn[1] - Xn[0]
        elif i >= Nx:
            fc = f[+ 2 * (i - Nx) + 1]; d = Xn[-1] - Xn[-2]
        else:
            fc = 0;                     d = Xn[i + 1] - Xn[i]
        gPmlXc[i + l] = cdt / d * fc
    for j in range(-l, Ny + l):
        if   j <  0:
            fc = f[- 2 * (j     ) - 1]; d = Yn[1] - Yn[0]
        elif j >= Ny:
            fc = f[+ 2 * (j - Ny) + 1]; d = Yn[-1] - Yn[-2]
        else:
            fc = 0;                     d = Yn[j + 1] - Yn[j]
        gPmlYc[j + l] = cdt / d * fc
    for k in range(-l, Nz + l):
        if   k <  0:
            fc = f[- 2 * (k     ) - 1]; d = Zn[1] - Zn[0]
        elif k >= Nz:
            fc = f[+ 2 * (k - Nz) + 1]; d = Zn[-1] - Zn[-2]
        else:
            fc = 0;                     d = Zn[k + 1] - Zn[k]
        gPmlZc[k + l] = cdt / d * fc

    # free
    f = None

    # 境界媒質の係数(の逆数)
    nmaterial = fMaterial.shape[0]
    rPmlE = np.zeros(nmaterial, f_dtype)
    rPmlH = np.zeros(nmaterial, f_dtype)
    rPml  = np.zeros(nmaterial, f_dtype)
    for m in range(nmaterial):
        epsr = fMaterial[m, 0]
        amur = fMaterial[m, 2]
        rPmlE[m] = 1 / epsr
        rPmlH[m] = 1 / amur
        rPml[m] = 1 / math.sqrt(epsr * amur)

    # 以下は必要
    PEC = Parm['PEC']
    rPmlE[PEC] = 0
    rPmlH[PEC] = 0
    rPml[PEC] = 0

    return gPmlXn, gPmlYn, gPmlZn, gPmlXc, gPmlYc, gPmlZc, rPmlE, rPmlH, rPml

# メッシュ因子
def _meshfactor(Parm, Nx, Ny, Nz, Xn, Yn, Zn):
    # mesh factor : c * dt / d
    f_dtype = Parm['f_dtype']
    cdt = Parm['C'] * Parm['dt']

    # alloc
    RXn = np.zeros(Nx + 1, f_dtype)
    RYn = np.zeros(Ny + 1, f_dtype)
    RZn = np.zeros(Nz + 1, f_dtype)
    RXc = np.empty(Nx + 0, f_dtype)
    RYc = np.empty(Ny + 0, f_dtype)
    RZc = np.empty(Nz + 0, f_dtype)

    # 節点
    RXn[1:-1] = cdt / ((Xn[2:] - Xn[:-2]) / 2)
    RYn[1:-1] = cdt / ((Yn[2:] - Yn[:-2]) / 2)
    RZn[1:-1] = cdt / ((Zn[2:] - Zn[:-2]) / 2)
    RXn[0]  = cdt / (Xn[1]  - Xn[0])
    RYn[0]  = cdt / (Yn[1]  - Yn[0])
    RZn[0]  = cdt / (Zn[1]  - Zn[0])
    RXn[-1] = cdt / (Xn[-1] - Xn[-2])
    RYn[-1] = cdt / (Yn[-1] - Yn[-2])
    RZn[-1] = cdt / (Zn[-1] - Zn[-2])

    # セル中心
    RXc = cdt / (Xn[1:] - Xn[:-1])
    RYc = cdt / (Yn[1:] - Yn[:-1])
    RZc = cdt / (Zn[1:] - Zn[:-1])

    return RXn, RYn, RZn, RXc, RYc, RZc

# 分散性媒質の係数を用意する
@jit(cache=True, nogil=True, parallel=True, nopython=True)
def _dispersion(iMaterial, fMaterial, dt):

    for m in prange(fMaterial.shape[0]):
        if iMaterial[m] == 2:
            einf = fMaterial[m, 4]
            ae   = fMaterial[m, 5]
            be   = fMaterial[m, 6]
            ce   = fMaterial[m, 7]
            ke = math.exp(-ce * dt)
            xi0 = (ae * dt) + (be / ce) * (1 - ke)
            dxi0 = (be / ce) * (1 - ke) * (1 - ke)
            #print(einf, ae, be, ce)
            #print(ce * dt, ke, xi0, dxi0)
            fMaterial[m,  8] = 1 / (einf + xi0)
            fMaterial[m,  9] = dxi0
            fMaterial[m, 10] = ke
    #print(fMaterial)
