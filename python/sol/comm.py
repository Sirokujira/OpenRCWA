# -*- coding: utf-8 -*-
"""
comm.py (MPI)
"""

import numpy as np
from numba import jit
from mpi4py import MPI
import sol.setup

# 入力データを全プロセスで共有する
def broadcast(Parm,
    Nx, Ny, Nz, Xn, Yn, Zn,
    iMaterial, fMaterial, iGeometry, fGeometry,
    iFeed, fFeed, iPoint, fPoint, iInductor, fInductor, Freq1, Freq2):

    i_buf = None
    d_buf = None
    i_num = np.zeros(1, 'i4')
    d_num = np.zeros(1, 'i4')

    # rootに変数を格納する
    if Parm['comm_rank'] == 0:
        #print(iMaterial.size, iGeometry.size)
        #print(fMaterial.size, fGeometry.size)
        #print(Xn.size, Yn.size, Zn.size)
        #print(iFeed.size, fFeed.size)
        #print(Freq1.size, Freq2.size)

        # データ数
        i_num[0] = 20 \
            + iMaterial.size + iGeometry.size \
            + iFeed.size + iPoint.size + iInductor.size
        d_num[0] = 8 \
            + Xn.size + Yn.size + Zn.size + fMaterial.size + fGeometry.size \
            + fFeed.size + fPoint.size + fInductor.size \
            + Freq1.size + Freq2.size

        # alloc
        i_buf = np.zeros(i_num, 'i4')
        d_buf = np.zeros(d_num, 'f8')

        i_id = 0
        d_id = 0

        # 整数データ(20個)
        i_buf[i_id] = Nx;                   i_id += 1
        i_buf[i_id] = Ny;                   i_id += 1
        i_buf[i_id] = Nz;                   i_id += 1
        i_buf[i_id] = fMaterial.shape[0];   i_id += 1
        i_buf[i_id] = fGeometry.shape[0];   i_id += 1
        i_buf[i_id] = fFeed.shape[0];       i_id += 1
        i_buf[i_id] = fPoint.shape[0];      i_id += 1
        i_buf[i_id] = fInductor.shape[0];   i_id += 1
        i_buf[i_id] = len(Freq1);           i_id += 1
        i_buf[i_id] = len(Freq2);           i_id += 1
        i_buf[i_id] = Parm['source'];       i_id += 1
        i_buf[i_id] = Parm['planewave'][2]; i_id += 1
        i_buf[i_id] = Parm['abc'][0];       i_id += 1
        i_buf[i_id] = Parm['abc'][1];       i_id += 1
        i_buf[i_id] = Parm['pbc'][0];       i_id += 1
        i_buf[i_id] = Parm['pbc'][1];       i_id += 1
        i_buf[i_id] = Parm['pbc'][2];       i_id += 1
        i_buf[i_id] = Parm['prop'];         i_id += 1
        i_buf[i_id] = Parm['solver'][0];    i_id += 1
        i_buf[i_id] = Parm['solver'][1];    i_id += 1

        # 実数データ(8個)
        d_buf[d_id] = Parm['planewave'][0]; d_id += 1
        d_buf[d_id] = Parm['planewave'][1]; d_id += 1
        d_buf[d_id] = Parm['rfeed'];        d_id += 1
        d_buf[d_id] = Parm['abc'][2];       d_id += 1
        d_buf[d_id] = Parm['abc'][3];       d_id += 1
        d_buf[d_id] = Parm['solver'][2];    d_id += 1
        d_buf[d_id] = Parm['dt'];           d_id += 1
        d_buf[d_id] = Parm['tw'];           d_id += 1

        for i in range(len(Xn)):
            d_buf[d_id] = Xn[i]; d_id += 1
        for j in range(len(Yn)):
            d_buf[d_id] = Yn[j]; d_id += 1
        for k in range(len(Zn)):
            d_buf[d_id] = Zn[k]; d_id += 1

        for n in range(iMaterial.shape[0]):
            i_buf[i_id] = iMaterial[n]; i_id += 1
        for n in range(fMaterial.shape[0]):
            for m in range(fMaterial.shape[1]):
                d_buf[d_id] = fMaterial[n, m]; d_id += 1

        for n in range(iGeometry.shape[0]):
            for m in range(iGeometry.shape[1]):
                i_buf[i_id] = iGeometry[n, m]; i_id += 1
        for n in range(fGeometry.shape[0]):
            for m in range(fGeometry.shape[1]):
                d_buf[d_id] = fGeometry[n, m]; d_id += 1

        for n in range(iFeed.shape[0]):
            for m in range(iFeed.shape[1]):
                i_buf[i_id] = iFeed[n, m]; i_id += 1
        for n in range(fFeed.shape[0]):
            for m in range(fFeed.shape[1]):
                d_buf[d_id] = fFeed[n, m]; d_id += 1

        for n in range(iPoint.shape[0]):
            for m in range(iPoint.shape[1]):
                i_buf[i_id] = iPoint[n, m]; i_id += 1
        for n in range(fPoint.shape[0]):
            for m in range(fPoint.shape[1]):
                d_buf[d_id] = fPoint[n, m]; d_id += 1

        for n in range(iInductor.shape[0]):
            for m in range(iInductor.shape[1]):
                i_buf[i_id] = iInductor[n, m]; i_id += 1
        for n in range(fInductor.shape[0]):
            for m in range(fInductor.shape[1]):
                d_buf[d_id] = fInductor[n, m]; d_id += 1

        for _, freq1 in enumerate(Freq1):
            d_buf[d_id] = freq1; d_id += 1

        for _, freq2 in enumerate(Freq2):
            d_buf[d_id] = freq2; d_id += 1

        # check
        assert(i_id == i_num[0])
        assert(d_id == d_num[0])

    # broadcast (root to non-root)
    MPI.COMM_WORLD.Bcast(i_num)
    MPI.COMM_WORLD.Bcast(d_num)
    #print(Parm['comm_size'], Parm['comm_rank'], i_num[0], d_num[0])

    # alloc
    if Parm['comm_rank'] > 0:
        i_buf = np.zeros(i_num[0], 'i4')
        d_buf = np.zeros(d_num[0], 'f8')

    MPI.COMM_WORLD.Bcast(i_buf)
    MPI.COMM_WORLD.Bcast(d_buf)

    # 受信したデータを変数に代入する (非root)
    if Parm['comm_rank'] > 0:
        i_id = 0
        d_id = 0

        # パラメーター初期化
        Parm['source']     = 0
        Parm['planewave']  = [0] * 3
        Parm['rfeed']      = 0
        Parm['abc']        = [0] * 4
        Parm['pbc']        = [0] * 3
        Parm['prop']       = 0
        Parm['solver']     = [1000, 50, 1e-3]
        Parm['dt']         = 0
        Parm['tw']         = 0

        # 整数データ(20個)
        Nx = i_buf[i_id];                             i_id += 1
        Ny = i_buf[i_id];                             i_id += 1
        Nz = i_buf[i_id];                             i_id += 1
        iMaterial = np.zeros( i_buf[i_id],     'i4'); i_id += 1
        iGeometry = np.zeros((i_buf[i_id], 2), 'i4'); i_id += 1
        iFeed     = np.zeros((i_buf[i_id], 4), 'i4'); i_id += 1
        iPoint    = np.zeros((i_buf[i_id], 4), 'i4'); i_id += 1
        iInductor = np.zeros((i_buf[i_id], 4), 'i4'); i_id += 1
        Freq1     = np.zeros( i_buf[i_id],     'f8'); i_id += 1
        Freq2     = np.zeros( i_buf[i_id],     'f8'); i_id += 1
        Parm['source']       = i_buf[i_id];           i_id += 1
        Parm['planewave'][2] = i_buf[i_id];           i_id += 1
        Parm['abc'][0]       = i_buf[i_id];           i_id += 1
        Parm['abc'][1]       = i_buf[i_id];           i_id += 1
        Parm['pbc'][0]       = i_buf[i_id];           i_id += 1
        Parm['pbc'][1]       = i_buf[i_id];           i_id += 1
        Parm['pbc'][2]       = i_buf[i_id];           i_id += 1
        Parm['prop']         = i_buf[i_id];           i_id += 1
        Parm['solver'][0]    = i_buf[i_id];           i_id += 1
        Parm['solver'][1]    = i_buf[i_id];           i_id += 1

        fMaterial = np.zeros((iMaterial.shape[0], 11), 'f8')
        fGeometry = np.zeros((iGeometry.shape[0],  8), 'f8')
        fFeed     = np.zeros((iFeed.shape[0],      9), 'f8')
        fPoint    = np.zeros((iPoint.shape[0],     6), 'f8')
        fInductor = np.zeros((iInductor.shape[0], 10), 'f8')

        # 実数データ(8個)
        Parm['planewave'][0] = d_buf[d_id]; d_id += 1
        Parm['planewave'][1] = d_buf[d_id]; d_id += 1
        Parm['rfeed']        = d_buf[d_id]; d_id += 1
        Parm['abc'][2]       = d_buf[d_id]; d_id += 1
        Parm['abc'][3]       = d_buf[d_id]; d_id += 1
        Parm['solver'][2]    = d_buf[d_id]; d_id += 1
        Parm['dt']           = d_buf[d_id]; d_id += 1
        Parm['tw']           = d_buf[d_id]; d_id += 1
        #print(Parm['comm_rank'], Parm)

        Xn = np.zeros(Nx + 1, 'f8')
        Yn = np.zeros(Ny + 1, 'f8')
        Zn = np.zeros(Nz + 1, 'f8')
        for i in range(Nx + 1):
            Xn[i] = d_buf[d_id]; d_id += 1
        for j in range(Ny + 1):
            Yn[j] = d_buf[d_id]; d_id += 1
        for k in range(Nz + 1):
            Zn[k] = d_buf[d_id]; d_id += 1
        #print(Parm['comm_rank'], Xn, Yn, Zn)

        for n in range(iMaterial.shape[0]):
            iMaterial[n] = i_buf[i_id]; i_id += 1
        for n in range(fMaterial.shape[0]):
            for m in range(fMaterial.shape[1]):
                fMaterial[n, m] = d_buf[d_id]; d_id += 1
        #print(Parm['comm_rank'], iMaterial, fMaterial)

        for n in range(iGeometry.shape[0]):
            for m in range(iGeometry.shape[1]):
                iGeometry[n, m] = i_buf[i_id]; i_id += 1
        for n in range(fGeometry.shape[0]):
            for m in range(fGeometry.shape[1]):
                fGeometry[n, m] = d_buf[d_id]; d_id += 1
        #print(Parm['comm_rank'], iGeometry, fGeometry)

        for n in range(iFeed.shape[0]):
            for m in range(iFeed.shape[1]):
                iFeed[n, m] = i_buf[i_id]; i_id += 1
        for n in range(fFeed.shape[0]):
            for m in range(fFeed.shape[1]):
                fFeed[n, m] = d_buf[d_id]; d_id += 1

        for n in range(iPoint.shape[0]):
            for m in range(iPoint.shape[1]):
                iPoint[n, m] = i_buf[i_id]; i_id += 1
        for n in range(fPoint.shape[0]):
            for m in range(fPoint.shape[1]):
                fPoint[n, m] = d_buf[d_id]; d_id += 1

        for n in range(iInductor.shape[0]):
            for m in range(iInductor.shape[1]):
                iInductor[n, m] = i_buf[i_id]; i_id += 1
        for n in range(fInductor.shape[0]):
            for m in range(fInductor.shape[1]):
                fInductor[n, m] = d_buf[d_id]; d_id += 1

        for n in range(len(Freq1)):
            Freq1[n] = d_buf[d_id]; d_id += 1

        for n in range(len(Freq2)):
            Freq2[n] = d_buf[d_id]; d_id += 1
        #print(Parm['comm_rank'], Freq1, Freq2)

        # check
        assert(i_id == i_num[0])
        assert(d_id == d_num[0])

    # free
    i_buf = None
    d_buf = None

    return \
    Nx, Ny, Nz, Xn, Yn, Zn, \
    iMaterial, fMaterial, iGeometry, fGeometry, \
    iFeed, fFeed, iPoint, fPoint, iInductor, fInductor, Freq1, Freq2

# X境界
def boundary_X(Parm, Npx, Npy, Npz, Ipx, Ipy, Ipz, iMin, iMax, jMin, jMax, kMin, kMax):

    la = 1 if Parm['abc'][0] == 0 else Parm['abc'][1]
    f_dtype = Parm['f_dtype']

    Bx_jhy = np.zeros(2, 'i4')
    Bx_jhz = np.zeros(2, 'i4')
    Bx_khy = np.zeros(2, 'i4')
    Bx_khz = np.zeros(2, 'i4')
    Bx_jhy[0] = jMin + ((- la + 1) if (Ipy == 0      ) else 0)
    Bx_jhy[1] = jMax + ((+ la - 1) if (Ipy == Npy - 1) else 0)
    Bx_khy[0] = kMin + ((- la    ) if (Ipz == 0      ) else 0)
    Bx_khy[1] = kMax + ((+ la - 1) if (Ipz == Npz - 1) else -1)
    Bx_jhz[0] = jMin + ((- la    ) if (Ipy == 0      ) else 0)
    Bx_jhz[1] = jMax + ((+ la - 1) if (Ipy == Npy - 1) else -1)
    Bx_khz[0] = kMin + ((- la + 1) if (Ipz == 0      ) else 0)
    Bx_khz[1] = kMax + ((+ la - 1) if (Ipz == Npz - 1) else 0)

    numhy = (Bx_jhy[1] - Bx_jhy[0] + 1) \
          * (Bx_khy[1] - Bx_khy[0] + 1)
    numhz = (Bx_jhz[1] - Bx_jhz[0] + 1) \
          * (Bx_khz[1] - Bx_khz[0] + 1)
    #size_t size_hy_x = Bx_numhy_x * sizeof(real_t);
    #size_t size_hz_x = Bx_numhz_x * sizeof(real_t);
    SendBuf_Bx_hy = np.zeros(numhy, f_dtype)
    RecvBuf_Bx_hy = np.zeros(numhy, f_dtype)
    SendBuf_Bx_hz = np.zeros(numhz, f_dtype)
    RecvBuf_Bx_hz = np.zeros(numhz, f_dtype)
    #Bid.ip[0] = iMin - 1;  // - boundary - 1 (recv)
    #Bid.ip[1] = iMin;      // - boundary     (send)
    #Bid.ip[2] = iMax - 1;  // + boundary - 1 (send)
    #Bid.ip[3] = iMax;      // + boundary     (recv)
    #print(Bx_jhy, Bx_jhz, Bx_khy, Bx_khz)

    return \
    Bx_jhy, Bx_jhz, Bx_khy, Bx_khz, \
    SendBuf_Bx_hy, SendBuf_Bx_hz, RecvBuf_Bx_hy, RecvBuf_Bx_hz

# Y境界
def boundary_Y(Parm, Npx, Npy, Npz, Ipx, Ipy, Ipz, iMin, iMax, jMin, jMax, kMin, kMax):

    la = 1 if Parm['abc'][0] == 0 else Parm['abc'][1]
    f_dtype = Parm['f_dtype']

    By_khz = np.zeros(2, 'i4')
    By_khx = np.zeros(2, 'i4')
    By_ihz = np.zeros(2, 'i4')
    By_ihx = np.zeros(2, 'i4')
    By_khz[0] = kMin + ((- la + 1) if (Ipz == 0      ) else 0)
    By_khz[1] = kMax + ((+ la - 1) if (Ipz == Npz - 1) else 0)
    By_ihz[0] = iMin + ((- la    ) if (Ipx == 0      ) else 0)
    By_ihz[1] = iMax + ((+ la - 1) if (Ipx == Npx - 1) else -1)
    By_khx[0] = kMin + ((- la    ) if (Ipz == 0      ) else 0)
    By_khx[1] = kMax + ((+ la - 1) if (Ipz == Npz - 1) else -1)
    By_ihx[0] = iMin + ((- la + 1) if (Ipx == 0      ) else 0)
    By_ihx[1] = iMax + ((+ la - 1) if (Ipx == Npx - 1) else 0)

    numhz = (By_khz[1] - By_khz[0] + 1) \
          * (By_ihz[1] - By_ihz[0] + 1)
    numhx = (By_khx[1] - By_khx[0] + 1) \
          * (By_ihx[1] - By_ihx[0] + 1)
    SendBuf_By_hz = np.zeros(numhz, f_dtype)
    RecvBuf_By_hz = np.zeros(numhz, f_dtype)
    SendBuf_By_hx = np.zeros(numhx, f_dtype)
    RecvBuf_By_hx = np.zeros(numhx, f_dtype)

    return \
    By_khz, By_khx, By_ihz, By_ihx, \
    SendBuf_By_hz, SendBuf_By_hx, RecvBuf_By_hz, RecvBuf_By_hx

# Z境界
def boundary_Z(Parm, Npx, Npy, Npz, Ipx, Ipy, Ipz, iMin, iMax, jMin, jMax, kMin, kMax):

    la = 1 if Parm['abc'][0] == 0 else Parm['abc'][1]
    f_dtype = Parm['f_dtype']

    Bz_ihx = np.zeros(2, 'i4')
    Bz_ihy = np.zeros(2, 'i4')
    Bz_jhx = np.zeros(2, 'i4')
    Bz_jhy = np.zeros(2, 'i4')
    Bz_ihx[0] = iMin + ((- la + 1) if (Ipx == 0      ) else 0)
    Bz_ihx[1] = iMax + ((+ la - 1) if (Ipx == Npx - 1) else 0)
    Bz_jhx[0] = jMin + ((- la    ) if (Ipy == 0      ) else 0)
    Bz_jhx[1] = jMax + ((+ la - 1) if (Ipy == Npy - 1) else -1)
    Bz_ihy[0] = iMin + ((- la    ) if (Ipx == 0      ) else 0)
    Bz_ihy[1] = iMax + ((+ la - 1) if (Ipx == Npx - 1) else -1)
    Bz_jhy[0] = jMin + ((- la + 1) if (Ipy == 0      ) else 0)
    Bz_jhy[1] = jMax + ((+ la - 1) if (Ipy == Npy - 1) else 0)

    numhx = (Bz_ihx[1] - Bz_ihx[0] + 1) \
          * (Bz_jhx[1] - Bz_jhx[0] + 1)
    numhy = (Bz_ihy[1] - Bz_ihy[0] + 1) \
          * (Bz_jhy[1] - Bz_jhy[0] + 1)
    SendBuf_Bz_hx = np.zeros(numhx, f_dtype)
    RecvBuf_Bz_hx = np.zeros(numhx, f_dtype)
    SendBuf_Bz_hy = np.zeros(numhy, f_dtype)
    RecvBuf_Bz_hy = np.zeros(numhy, f_dtype)

    return \
    Bz_ihx, Bz_ihy, Bz_jhx, Bz_jhy, \
    SendBuf_Bz_hx, SendBuf_Bz_hy, RecvBuf_Bz_hx, RecvBuf_Bz_hy

# 和(スカラー, 全プロセスで共有)
def sum_scalar(var):
    sendbuf = np.array([var], 'f8')
    recvbuf = np.zeros(1, 'f8')
    MPI.COMM_WORLD.Allreduce(sendbuf, recvbuf)

    return recvbuf[0]

# 全領域の電磁界をrootに集める
def gather_near3d(Parm, Freq2,
    cEx, cEy, cEz, cHx, cHy, cHz,
    Nx, Ny, Nz, iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0, NN):

    comm_size = Parm['comm_size']
    comm_rank = Parm['comm_rank']
    #print(comm_size, comm_rank)

    if (NN <= 0) or (len(Freq2) <= 0):
        return

    g_cEx = g_cEy = g_cEz = g_cHx = g_cHy = g_cHz = None
    
    # root
    if comm_rank == 0:
        # rootの部分領域のindexを保存しておく
        imin = iMin
        imax = iMax
        jmin = jMin
        jmax = jMax
        kmin = kMin
        kmax = kMax
        ni = Ni
        nj = Nj
        nk = Nk
        n0 = N0
        nn = NN
        #print(iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0, NN)

        # root : 全領域用のindex
        iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0, NN, _, _, _ \
        = sol.setup.getIndex(Parm, Nx, Ny, Nz, 1, 1, 1, 1, 1, 1, 0)
        #print(iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0, NN)

        fsize = NN * len(Freq2)
        g_cEx = np.zeros(fsize, 'c8')
        g_cEy = np.zeros(fsize, 'c8')
        g_cEz = np.zeros(fsize, 'c8')
        g_cHx = np.zeros(fsize, 'c8')
        g_cHy = np.zeros(fsize, 'c8')
        g_cHz = np.zeros(fsize, 'c8')

        # 部分配列を全体配列にコピーする
        for ifreq in range(len(Freq2)):
            g_start = ifreq * NN  # 全体配列の開始アドレス
            l_start = ifreq * nn  # 部分配列の全体配列内での開始アドレス

            """for i in range(imin, imax + 0):
                for j in range(jmin, jmax + 1):
                    for k in range(kmin, kmax + 1):
                        g_n = (Ni * i) + (Nj * j) + (Nk * k) + N0 + (ifreq * NN)
                        n   = (ni * i) + (nj * j) + (nk * k) + n0 + (ifreq * nn)
                        g_cEx[g_n] = cEx[n]"""
            # Ex
            _copy3d(g_cEx, cEx,
                imin + 0, imax + 0, jmin + 0, jmax + 1, kmin + 0, kmax + 1,
                Ni, Nj, Nk, N0, ni, nj, nk, n0, g_start, l_start)
            # Ey
            _copy3d(g_cEy, cEy,
                imin + 0, imax + 1, jmin + 0, jmax + 0, kmin + 0, kmax + 1,
                Ni, Nj, Nk, N0, ni, nj, nk, n0, g_start, l_start)
            # Ez
            _copy3d(g_cEz, cEz,
                imin + 0, imax + 1, jmin + 0, jmax + 1, kmin + 0, kmax + 0,
                Ni, Nj, Nk, N0, ni, nj, nk, n0, g_start, l_start)
            # Hx
            _copy3d(g_cHx, cHx,
                imin - 0, imax + 1, jmin - 1, jmax + 1, kmin - 1, kmax + 1,
                Ni, Nj, Nk, N0, ni, nj, nk, n0, g_start, l_start)
            # Hy
            _copy3d(g_cHy, cHy,
                imin - 1, imax + 1, jmin - 0, jmax + 1, kmin - 1, kmax + 1,
                Ni, Nj, Nk, N0, ni, nj, nk, n0, g_start, l_start)
            # Hz
            _copy3d(g_cHz, cHz,
                imin - 1, imax + 1, jmin - 1, jmax + 1, kmin - 0, kmax + 1,
                Ni, Nj, Nk, N0, ni, nj, nk, n0, g_start, l_start)

            """
            g_n = ifreq * NN
            n   = ifreq * NN_root
            g_cEx[g_n: g_n + NN_root] = cEx[n: n + NN_root].copy()
            g_cEy[g_n: g_n + NN_root] = cEy[n: n + NN_root].copy()
            g_cEz[g_n: g_n + NN_root] = cEz[n: n + NN_root].copy()
            g_cHx[g_n: g_n + NN_root] = cHx[n: n + NN_root].copy()
            g_cHy[g_n: g_n + NN_root] = cHy[n: n + NN_root].copy()
            g_cHz[g_n: g_n + NN_root] = cHz[n: n + NN_root].copy()
            """

    for ifreq in range(len(Freq2)):

        if comm_rank == 0:
            # root : 受信

            for rank in range(1, comm_size):

                # 受信
                irecv = np.zeros(11, 'i4')
                MPI.COMM_WORLD.Recv(irecv, source=rank)

                imin, imax, jmin, jmax, kmin, kmax, ni, nj, nk, n0, nn = irecv
                #print(imin, imax, jmin, jmax, kmin, kmax, ni, nj, nk, n0, nn)
                
                # 受信用バッファー作成
                recv = np.zeros(nn, 'c8')

                # 開始アドレス
                g_start = ifreq * NN
                l_start = 0

                # 受信後全体配列に代入する

                # Ex
                MPI.COMM_WORLD.Recv(recv, source=rank)
                _copy3d(g_cEx, recv,
                    imin + 0, imax + 0, jmin + 0, jmax + 1, kmin + 0, kmax + 1,
                    Ni, Nj, Nk, N0, ni, nj, nk, n0, g_start, l_start)
                """
                for i in range(imin, imax + 0):
                    for j in range(jmin, jmax + 1):
                        for k in range(kmin, kmax + 1):
                            g_n = (Ni * i) + (Nj * j) + (Nk * k) + N0 + (ifreq * NN)
                            n   = (ni * i) + (nj * j) + (nk * k) + n0
                            g_cEx[g_n] = recv[n]"""

                # Ey
                MPI.COMM_WORLD.Recv(recv, source=rank)
                _copy3d(g_cEy, recv,
                    imin + 0, imax + 1, jmin + 0, jmax + 0, kmin + 0, kmax + 1,
                    Ni, Nj, Nk, N0, ni, nj, nk, n0, g_start, l_start)
                """
                for i in range(imin, imax + 1):
                    for j in range(jmin, jmax + 0):
                        for k in range(kmin, kmax + 1):
                            g_n = (Ni * i) + (Nj * j) + (Nk * k) + N0 + (ifreq * NN)
                            n   = (ni * i) + (nj * j) + (nk * k) + n0
                            g_cEy[g_n] = recv[n]"""

                # Ey
                MPI.COMM_WORLD.Recv(recv, source=rank)
                _copy3d(g_cEz, recv,
                    imin + 0, imax + 1, jmin + 0, jmax + 1, kmin + 0, kmax + 0,
                    Ni, Nj, Nk, N0, ni, nj, nk, n0, g_start, l_start)
                """
                for i in range(imin, imax + 1):
                    for j in range(jmin, jmax + 1):
                        for k in range(kmin, kmax + 0):
                            g_n = (Ni * i) + (Nj * j) + (Nk * k) + N0 + (ifreq * NN)
                            n   = (ni * i) + (nj * j) + (nk * k) + n0
                            g_cEz[g_n] = recv[n]"""

                # Hx
                MPI.COMM_WORLD.Recv(recv, source=rank)
                _copy3d(g_cHx, recv,
                    imin - 0, imax + 1, jmin - 1, jmax + 1, kmin - 1, kmax + 1,
                    Ni, Nj, Nk, N0, ni, nj, nk, n0, g_start, l_start)
                """
                for i in range(imin - 0, imax + 1):
                    for j in range(jmin - 1, jmax + 1):
                        for k in range(kmin - 1, kmax + 1):
                            g_n = (Ni * i) + (Nj * j) + (Nk * k) + N0 + (ifreq * NN)
                            n   = (ni * i) + (nj * j) + (nk * k) + n0
                            g_cHx[g_n] = recv[n]"""

                # Hy
                MPI.COMM_WORLD.Recv(recv, source=rank)
                _copy3d(g_cHy, recv,
                    imin - 1, imax + 1, jmin - 0, jmax + 1, kmin - 1, kmax + 1,
                    Ni, Nj, Nk, N0, ni, nj, nk, n0, g_start, l_start)
                """
                for i in range(imin - 1, imax + 1):
                    for j in range(jmin - 0, jmax + 1):
                        for k in range(kmin - 1, kmax + 1):
                            g_n = (Ni * i) + (Nj * j) + (Nk * k) + N0 + (ifreq * NN)
                            n   = (ni * i) + (nj * j) + (nk * k) + n0
                            g_cHy[g_n] = recv[n]"""

                # Hz
                MPI.COMM_WORLD.Recv(recv, source=rank)
                _copy3d(g_cHz, recv,
                    imin - 1, imax + 1, jmin - 1, jmax + 1, kmin - 0, kmax + 1,
                    Ni, Nj, Nk, N0, ni, nj, nk, n0, g_start, l_start)
                """
                for i in range(imin - 1, imax + 1):
                    for j in range(jmin - 1, jmax + 1):
                        for k in range(kmin - 0, kmax + 1):
                            g_n = (Ni * i) + (Nj * j) + (Nk * k) + N0 + (ifreq * NN)
                            n   = (ni * i) + (nj * j) + (nk * k) + n0
                            g_cHz[g_n] = recv[n]"""

                # free
                recv = None

        else:
            # 非root : rootに送信する

            # indexを送信する(11個)
            isend = np.array([iMin, iMax, jMin, jMax, kMin, kMax, \
                              Ni, Nj, Nk, N0, NN]).astype('i4')
            MPI.COMM_WORLD.Send(isend, 0)

            # データ本体を送信する
            #print(NN)
            n0 = (ifreq + 0) * NN
            n1 = (ifreq + 1) * NN
            MPI.COMM_WORLD.Send(cEx[n0: n1], 0)
            MPI.COMM_WORLD.Send(cEy[n0: n1], 0)
            MPI.COMM_WORLD.Send(cEz[n0: n1], 0)
            MPI.COMM_WORLD.Send(cHx[n0: n1], 0)
            MPI.COMM_WORLD.Send(cHy[n0: n1], 0)
            MPI.COMM_WORLD.Send(cHz[n0: n1], 0)
    
    
    # root : ポインターをコピーする
    if comm_rank == 0:
        #cEx = cEy = cEz = cHx = cHy = cHz = None
        cEx = g_cEx
        cEy = g_cEy
        cEz = g_cEz
        cHx = g_cHx
        cHy = g_cHy
        cHz = g_cHz

    return \
        cEx, cEy, cEz, cHx, cHy, cHz, \
        iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0, NN

# 3D配列のコピー
@jit(cache=True, nopython=True)
def _copy3d(g_a, l_a,
    i_min, i_max, j_min, j_max, k_min, k_max,
    Ni, Nj, Nk, N0, ni, nj, nk, n0, g_start, l_start):

    for i in range(i_min, i_max):
        for j in range(j_min, j_max):
            for k in range(k_min, k_max):
                g_n = (Ni * i) + (Nj * j) + (Nk * k) + N0 + g_start
                l_n = (ni * i) + (nj * j) + (nk * k) + n0 + l_start
                g_a[g_n] = l_a[l_n]

# send feed waveform to root process
def gather_feed(Parm, iFeed, VFeed, IFeed,
    Nx, Ny, Nz, iMin, iMax, jMin, jMax, kMin, kMax):

    comm_size = Parm['comm_size']
    comm_rank = Parm['comm_rank']

    count = Parm['solver'][0] + 1

    for ifeed in range(iFeed.shape[0]):
        i, j, k = iFeed[ifeed, 1:4]
        b = _mpi_proc(comm_size, comm_rank, i, j, k,
            Nx, Ny, Nz, iMin, iMax, jMin, jMax, kMin, kMax)
        # 非rootにあるときのみ通信する
        if   (comm_rank == 0) and (not b):
            MPI.COMM_WORLD.Recv(VFeed[ifeed, :count], source=MPI.ANY_SOURCE)
            MPI.COMM_WORLD.Recv(IFeed[ifeed, :count], source=MPI.ANY_SOURCE)
        elif (comm_rank  > 0) and b:
            MPI.COMM_WORLD.Send(VFeed[ifeed, :count], 0)
            MPI.COMM_WORLD.Send(IFeed[ifeed, :count], 0)

        MPI.COMM_WORLD.Barrier()

# send point waveform to root process
def gather_point(Parm, iPoint, VPoint,
    Nx, Ny, Nz, iMin, iMax, jMin, jMax, kMin, kMax):

    comm_size = Parm['comm_size']
    comm_rank = Parm['comm_rank']

    count = Parm['solver'][0] + 1

    for ipoint in range(iPoint.shape[0]):
        i, j, k = iPoint[ipoint, 1:4]
        b = _mpi_proc(comm_size, comm_rank, i, j, k,
            Nx, Ny, Nz, iMin, iMax, jMin, jMax, kMin, kMax)
        # 非rootにあるときのみ通信する
        if   (comm_rank == 0) and (not b):
            MPI.COMM_WORLD.Recv(VPoint[ipoint, :count], source=MPI.ANY_SOURCE)
        elif (comm_rank  > 0) and b:
            MPI.COMM_WORLD.Send(VPoint[ipoint, :count], 0)

        MPI.COMM_WORLD.Barrier()

# MPI buffer size
# TODO 削除
def buffer_size(Parm, Ny, Nz, iMin, iMax, Ni, Nj, Nk, N0):

    # ABC層数
    #lx = 1 if (Parm['abc'][0] == 0) else Parm['abc'][1]
    ly = 1 if (Parm['abc'][0] == 0) else Parm['abc'][1]
    lz = 1 if (Parm['abc'][0] == 0) else Parm['abc'][1]

    jmin_hy = -ly + 1
    jmin_hz = -ly + 0
    kmin_hy = -lz + 0
    kmin_hz = -lz + 1
    #print(Parm['comm_rank'], jmin_hy, jmin_hz, kmin_hy, kmin_hz)

    jmax_hy = Ny + ly
    jmax_hz = Ny + ly
    kmax_hy = Nz + lz
    kmax_hz = Nz + lz
    #print(Parm['comm_rank'], jmax_hy, jmax_hz, kmax_hy, kmax_hz)

    # [0] : -X boundary - 1 (recv)
    # [1] : -X boundary (send)
    # [2] : +X boundary (send)
    # [3] : +X boundary + 1 (recv)

    offset_hy = np.zeros(4, 'i4')
    offset_hz = np.zeros(4, 'i4')

    #Offset_Hy[0] = NA(iMin - 1, jmin_hy, kmin_hy);
    #Offset_Hy[1] = NA(iMin,     jmin_hy, kmin_hy);
    #Offset_Hy[2] = NA(iMax - 1, jmin_hy, kmin_hy);
    #Offset_Hy[3] = NA(iMax,     jmin_hy, kmin_hy);
    offset_hy[0] = Ni * (iMin - 1) + Nj * jmin_hy + Nk * kmin_hy + N0
    offset_hy[1] = Ni * (iMin    ) + Nj * jmin_hy + Nk * kmin_hy + N0
    offset_hy[2] = Ni * (iMax - 1) + Nj * jmin_hy + Nk * kmin_hy + N0
    offset_hy[3] = Ni * (iMax    ) + Nj * jmin_hy + Nk * kmin_hy + N0

    #Offset_Hz[0] = NA(iMin - 1, jmin_hz, kmin_hz);
    #Offset_Hz[1] = NA(iMin,     jmin_hz, kmin_hz);
    #Offset_Hz[2] = NA(iMax - 1, jmin_hz, kmin_hz);
    #Offset_Hz[3] = NA(iMax,     jmin_hz, kmin_hz);
    offset_hz[0] = Ni * (iMin - 1) + Nj * jmin_hz + Nk * kmin_hz + N0
    offset_hz[1] = Ni * (iMin    ) + Nj * jmin_hz + Nk * kmin_hz + N0
    offset_hz[2] = Ni * (iMax - 1) + Nj * jmin_hz + Nk * kmin_hz + N0
    offset_hz[3] = Ni * (iMax    ) + Nj * jmin_hz + Nk * kmin_hz + N0

    length_hy = Nj * (jmax_hy - jmin_hy) + Nk * (kmax_hy - kmin_hy)
    length_hz = Nj * (jmax_hz - jmin_hz) + Nk * (kmax_hz - kmin_hz)

    #lbuf = length_hy + length_hz
    #print(Parm['comm_rank'], length_hy, length_hz, lbuf)
    #buf = np.zeros(lbuf, Parm['f_dtype'])

    return offset_hy, offset_hz, length_hy, length_hz

# 自分の領域は(i,j,k)を含むか
def _mpi_proc(comm_size, comm_rank, i, j, k,
    Nx, Ny, Nz, iMin, iMax, jMin, jMax, kMin, kMax):

    b1 = (iMin <= i) and (i < iMax) and \
         (jMin <= j) and (j < jMax) and \
         (kMin <= k) and (k < kMax)
    b2 = (comm_rank == comm_size - 1) and \
         ((i == Nx) or (j == Ny) or (k == Nz))
    return b1 or b2

    #return ((i >= iMin) and (i < iMax)) \
    #    or ((comm_rank == comm_size - 1) and (i == Nx))
"""
def mpi_proc(Parm, Nx, iMin, iMax):

    iProc = np.zeros(Nx + 1, 'i4')

    for i in range(Nx + 1):
        iProc[i] = ((i >= iMin) and (i < iMax)) \
            or ((Parm['comm_rank'] == Parm['comm_size'] - 1) and (i == Nx))
        #print(Parm['comm_size'], Parm['comm_rank'], i, iProc[i])

    return iProc
"""