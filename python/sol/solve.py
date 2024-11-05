# -*- coding: utf-8 -*-
"""
solve.py
RCWA法の主要部:反復計算
"""

import numpy as np
import sol.updateEx, sol.updateEy, sol.updateEz, sol.updateHx, sol.updateHy, sol.updateHz
import sol.MurH
import sol.PmlEx, sol.PmlEy, sol.PmlEz, sol.PmlHx, sol.PmlHy, sol.PmlHz
import sol.DispEx, sol.DispEy, sol.DispEz
import sol.PbcX, sol.PbcY, sol.PbcZ
import sol.feed, sol.inductor, sol.pointV, sol.Near3d, sol.average
import sol.monitor, sol.cputime
import sol.comm, sol.comm_X, sol.comm_Y, sol.comm_Z

def iteration(
    io, fp,
    Parm, fPlanewave, iFeed, fFeed, iPoint, fPoint, iInductor, fInductor, Freq2,
    Nx, Ny, Nz, Npx, Npy, Npz, Ipx, Ipy, Ipz,
    Xn, Yn, Zn, Xc, Yc, Zc,
    iEx, iEy, iEz, iHx, iHy, iHz,
    C1E, C2E, C1H, C2H,
    K1Ex, K2Ex, K1Ey, K2Ey, K1Ez, K2Ez,
    K1Hx, K2Hx, K1Hy, K2Hy, K1Hz, K2Hz,
    RXn, RYn, RZn, RXc, RYc, RZc,
    fMurHx, fMurHy, fMurHz, iMurHx, iMurHy, iMurHz,
    iPmlEx, iPmlEy, iPmlEz, iPmlHx, iPmlHy, iPmlHz,
    gPmlXn, gPmlYn, gPmlZn, gPmlXc, gPmlYc, gPmlZc, rPmlE, rPmlH, rPml,
    iDispEx, fDispEx, iDispEy, fDispEy, iDispEz, fDispEz,
    cEdft, cHdft,
    iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0, NN):

    # 計算制御パラメーター
    maxiter = Parm['solver'][0]
    nout    = Parm['solver'][1]
    converg = Parm['solver'][2]
    dt      = Parm['dt']
    abc     = Parm['abc'][0]
    pml_l   = Parm['abc'][1]
    pbc     = Parm['pbc']
    source  = Parm['source']
    #print(Parm['comm_rank'], maxiter, nout, converg, dt, abc, pml_l, pbc, source)

    # MPI用buffer作成
    #print(Npx, Npy, Npz)
    if Npx > 1:
        Bx_jhy, Bx_jhz, Bx_khy, Bx_khz, \
        SendBuf_Bx_hy, SendBuf_Bx_hz, RecvBuf_Bx_hy, RecvBuf_Bx_hz, \
        = sol.comm.boundary_X(Parm, Npx, Npy, Npz, Ipx, Ipy, Ipz, iMin, iMax, jMin, jMax, kMin, kMax)
    if Npy > 1:
        By_khz, By_khx, By_ihz, By_ihx, \
        SendBuf_By_hz, SendBuf_By_hx, RecvBuf_By_hz, RecvBuf_By_hx, \
        = sol.comm.boundary_Y(Parm, Npx, Npy, Npz, Ipx, Ipy, Ipz, iMin, iMax, jMin, jMax, kMin, kMax)
    if Npz > 1:
        Bz_ihx, Bz_ihy, Bz_jhx, Bz_jhy, \
        SendBuf_Bz_hx, SendBuf_Bz_hy, RecvBuf_Bz_hx, RecvBuf_Bz_hy  \
        = sol.comm.boundary_Z(Parm, Npx, Npy, Npz, Ipx, Ipy, Ipz, iMin, iMax, jMin, jMax, kMin, kMax)

    # 計算に必要な配列を確保して初期化する
    Ex, Ey, Ez, Hx, Hy, Hz, \
    Exy, Exz, Eyz, Eyx, Ezx, Ezy, \
    Hxy, Hxz, Hyz, Hyx, Hzx, Hzy, \
    cEx, cEy, cEz, cHx, cHy, cHz, \
    VFeed, IFeed, VPoint, \
    Eiter, Hiter, Iiter \
    = _initfield( \
    Parm, NN, iFeed, iPoint, Freq2, iPmlEx, iPmlEy, iPmlEz, iPmlHx, iPmlHy, iPmlHz)

    # 反復計算の変数を初期化する
    emax = hmax = 0
    converged = False
    Niter = 0
    tdft = 0

    # タイムステップに関する反復計算
    t = 0
    for itime in range(maxiter + 1):
        # 磁界を更新する
        t += 0.5 * dt
        #print(Parm['comm_rank'], itime, t)
        sol.updateHx.calHx(
            Parm, t, fPlanewave, Xn, Yc, Zc, Hx, Ey, Ez, iHx,
            C1H, C2H, K1Hx, K2Hx, RYc, RZc,
            iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)
        sol.updateHy.calHy(
            Parm, t, fPlanewave, Yn, Zc, Xc, Hy, Ez, Ex, iHy,
            C1H, C2H, K1Hy, K2Hy, RZc, RXc,
            iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)
        sol.updateHz.calHz(
            Parm, t, fPlanewave, Zn, Xc, Yc, Hz, Ex, Ey, iHz,
            C1H, C2H, K1Hz, K2Hz, RXc, RYc,
            iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)

        # 磁界に吸収境界条件を適用する
        if   abc == 0:
            # Mur
            sol.MurH.calcH(
                Hx, fMurHx, iMurHx, Ni, Nj, Nk, N0)
            sol.MurH.calcH(
                Hy, fMurHy, iMurHy, Ni, Nj, Nk, N0)
            sol.MurH.calcH(
                Hz, fMurHz, iMurHz, Ni, Nj, Nk, N0)
        elif abc == 1:
            # PML
            sol.PmlHx.calHx(
                Ny, Nz, Hx, Ey, Ez, Hxy, Hxz, RYc, RZc,
                iPmlHx, gPmlYc, gPmlZc, rPmlH, rPml, pml_l,
                Ni, Nj, Nk, N0)
            sol.PmlHy.calHy(
                Nz, Nx, Hy, Ez, Ex, Hyz, Hyx, RZc, RXc,
                iPmlHy, gPmlZc, gPmlXc, rPmlH, rPml, pml_l,
                Ni, Nj, Nk, N0)
            sol.PmlHz.calHz(
                Nx, Ny, Hz, Ex, Ey, Hzx, Hzy, RXc, RYc,
                iPmlHz, gPmlXc, gPmlYc, rPmlH, rPml, pml_l,
                Ni, Nj, Nk, N0)

        # 磁界に周期境界条件を適用する
        if pbc[0] == 1:  # X
            if Npx > 1:  # MPI
                sol.comm_X.share(1, Hy, Hz,
                    SendBuf_Bx_hy, SendBuf_Bx_hz, RecvBuf_Bx_hy, RecvBuf_Bx_hz,
                    Bx_jhy, Bx_jhz, Bx_khy, Bx_khz,
                    Npx, Npy, Npz, Ipx, Ipy, Ipz,
                    iMin, iMax, Ni, Nj, Nk, N0)
            else:
                sol.PbcX.x(Nx, Hy, Hz, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)
        if pbc[1] == 1:  # Y
            if Npy > 1:  # MPI
                sol.comm_Y.share(1, Hz, Hx,
                    SendBuf_By_hz, SendBuf_By_hx, RecvBuf_By_hz, RecvBuf_By_hx,
                    By_khz, By_khx, By_ihz, By_ihx,
                    Npx, Npy, Npz, Ipx, Ipy, Ipz,
                    jMin, jMax, Ni, Nj, Nk, N0)
            else:
                sol.PbcY.y(Ny, Hz, Hx, kMin, kMax, iMin, iMax, Ni, Nj, Nk, N0)
        if pbc[2] == 1:  # Z
            if Npz > 1:  # MPI
                sol.comm_Z.share(1, Hx, Hy,
                    SendBuf_Bz_hx, SendBuf_Bz_hy, RecvBuf_Bz_hx, RecvBuf_Bz_hy,
                    Bz_ihx, Bz_ihy, Bz_jhx, Bz_jhy,
                    Npx, Npy, Npz, Ipx, Ipy, Ipz,
                    kMin, kMax, Ni, Nj, Nk, N0)
            else:
                sol.PbcZ.z(Nz, Hx, Hy, iMin, iMax, jMin, jMax, Ni, Nj, Nk, N0)

        # MPI: 境界のHを共有する
        if Npx > 1:
            sol.comm_X.share(0, Hy, Hz,
                SendBuf_Bx_hy, SendBuf_Bx_hz, RecvBuf_Bx_hy, RecvBuf_Bx_hz,
                Bx_jhy, Bx_jhz, Bx_khy, Bx_khz,
                Npx, Npy, Npz, Ipx, Ipy, Ipz,
                iMin, iMax, Ni, Nj, Nk, N0)
        if Npy > 1:
            sol.comm_Y.share(0, Hz, Hx,
                SendBuf_By_hz, SendBuf_By_hx, RecvBuf_By_hz, RecvBuf_By_hx,
                By_khz, By_khx, By_ihz, By_ihx,
                Npx, Npy, Npz, Ipx, Ipy, Ipz,
                jMin, jMax, Ni, Nj, Nk, N0)
        if Npz > 1:
            sol.comm_Z.share(0, Hx, Hy,
                SendBuf_Bz_hx, SendBuf_Bz_hy, RecvBuf_Bz_hx, RecvBuf_Bz_hy,
                Bz_ihx, Bz_ihy, Bz_jhx, Bz_jhy,
                Npx, Npy, Npz, Ipx, Ipy, Ipz,
                kMin, kMax, Ni, Nj, Nk, N0)

        # 電界を更新する
        t += 0.5 * dt
        sol.updateEx.calEx(
            Parm, t, fPlanewave, Xc, Yn, Zn, Ex, Hy, Hz, iEx,
            C1E, C2E, K1Ex, K2Ex, RYn, RZn,
            iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)
        sol.updateEy.calEy(
            Parm, t, fPlanewave, Yc, Zn, Xn, Ey, Hz, Hx, iEy,
            C1E, C2E, K1Ey, K2Ey, RZn, RXn,
            iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)
        sol.updateEz.calEz(
            Parm, t, fPlanewave, Zc, Xn, Yn, Ez, Hx, Hy, iEz,
            C1E, C2E, K1Ez, K2Ez, RXn, RYn,
            iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)

        # 分散性媒質の電界を更新する
        if iDispEx.shape[0] > 0:
            sol.DispEx.calEx(
                source, t, fPlanewave, Xc, Yn, Zn, Ex, iDispEx, fDispEx,
                Ni, Nj, Nk, N0)
        if iDispEy.shape[0] > 0:
            sol.DispEy.calEy(
                source, t, fPlanewave, Yc, Zn, Xn, Ey, iDispEy, fDispEy,
                Ni, Nj, Nk, N0)
        if iDispEz.shape[0] > 0:
            sol.DispEz.calEz(
                source, t, fPlanewave, Zc, Xn, Yn, Ez, iDispEz, fDispEz, Ni, Nj, Nk, N0)

        # 電界に吸収境界条件を適用する
        if abc == 1:
            # PML
            sol.PmlEx.calEx(
                Ny, Nz, Ex, Hy, Hz, Exy, Exz, RYn, RZn,
                iPmlEx, gPmlYn, gPmlZn, rPmlE, rPml, pml_l,
                Ni, Nj, Nk, N0)
            sol.PmlEy.calEy(
                Nz, Nx, Ey, Hz, Hx, Eyz, Eyx, RZn, RXn,
                iPmlEy, gPmlZn, gPmlXn, rPmlE, rPml, pml_l,
                Ni, Nj, Nk, N0)
            sol.PmlEz.calEz(
                Nx, Ny, Ez, Hx, Hy, Ezx, Ezy, RXn, RYn,
                iPmlEz, gPmlXn, gPmlYn, rPmlE, rPml, pml_l,
                Ni, Nj, Nk, N0)

        # 給電点に電圧を印加し、V/I時間波形を保存する
        if iFeed.shape[0] > 0:
            sol.feed.evi(
                itime, Parm, iFeed, fFeed, VFeed, IFeed,
                Ex, Ey, Ez, Hx, Hy, Hz, iEx, iEy, iEz,
                iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)

        # inductorを計算する
        if iInductor.shape[0] > 0:
            sol.inductor.calcL(
                Parm, iInductor, fInductor,
                Ex, Ey, Ez, Hx, Hy, Hz,
                iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)

        # 観測点の電圧を保存する
        if iPoint.shape[0] > 0:
            sol.pointV.v(
                itime, Parm, fPlanewave, Xn, Yn, Zn, Xc, Yc, Zc, iPoint, fPoint, VPoint,
                Ex, Ey, Ez,
                iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)

        # 全領域の電磁界のDFTを更新する(計算時間がかかる)
        t0 = sol.cputime.t(Parm['comm_size'], 0)
        sol.Near3d.dft(
            itime, Freq2, cEdft, cHdft,
            Ex, Ey, Ez, Hx, Hy, Hz,
            cEx, cEy, cEz, cHx, cHy, cHz,
            iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0, NN)
        tdft += sol.cputime.t(Parm['comm_size'], 0) - t0


        # 平均電磁界を計算して収束判定する
        if (itime % nout == 0) or (itime == maxiter):
            # 電磁界の和を計算する
            esum, hsum = sol.average.calcA(
                Ex, Ey, Ez, Hx, Hy, Hz,
                iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)
            #print(Parm['comm_rank'], esum, hsum)

            # MPI:　平均電磁界の全プロセスの和を全プロセスで共有する(収束判定のため)
            if Parm['comm_size'] > 1:
                esum = sol.comm.sum_scalar(esum)
                hsum = sol.comm.sum_scalar(hsum)

            # 平均電磁界
            eave = esum / (4.0 * Nx * Ny * Nz)
            have = hsum / (2.0 * Nx * Ny * Nz)

            # 平均電磁界を保存する (ポスト処理用)
            Eiter[Niter] = eave
            Hiter[Niter] = have
            Iiter[Niter] = itime
            Niter += 1

            # 収束状況をlogに出力する
            if io:
                ostr = "%7d %.6f %.6f" % (itime, eave, have)
                sol.monitor.log1(fp, ostr)

            # 収束判定する
            emax = max(emax, eave)
            hmax = max(hmax, have)
            if (eave < emax * converg) and \
               (have < hmax * converg):
                converged = True
                break

    # 収束結果をlogに出力する
    if io:
        ostr = "    --- %s ---" % ("converged" if converged else "max steps")
        sol.monitor.log1(fp, ostr)

    # タイムステップ数
    Ntime = itime + 1

    # MPI: 後処理に必要な電磁界をrootに集める
    if Parm['comm_size'] > 1:
        # 給電点時間波形
        if iFeed.shape[0] > 0:
            sol.comm.gather_feed(Parm, iFeed, VFeed, IFeed,
                Nx, Ny, Nz, iMin, iMax, jMin, jMax, kMin, kMax)

        # 観測点時間波形
        if iPoint.shape[0] > 0:
            sol.comm.gather_point(Parm, iPoint, VPoint,
                Nx, Ny, Nz, iMin, iMax, jMin, jMax, kMin, kMax)

        # 全領域の電磁界をrootに集める
        cEx, cEy, cEz, cHx, cHy, cHz, \
        iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0, NN \
        = sol.comm.gather_near3d(Parm, Freq2,
            cEx, cEy, cEz, cHx, cHy, cHz,
            Nx, Ny, Nz, iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0, NN)

    # free
    Ex = Ey = Ez = Hx = Hy = Hz = \
    Exy = Exz = Eyz = Eyx = Ezx = Ezy = \
    Hxy = Hxz = Hyz = Hyx = Hzx = Hzy = \
    iEx = iEy = iEz = iHx = iHy = iHz = None

    # MPIではiMin等は変更されているので返す必要がある

    return \
    cEx, cEy, cEz, cHx, cHy, cHz, \
    VFeed, IFeed, VPoint, Ntime, \
    Eiter, Hiter, Iiter, Niter, \
    iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0, NN, \
    tdft

# (private) 計算に必要な配列を確保して初期化する
def _initfield(
    Parm, NN, iFeed, iPoint, Freq2, iPmlEx, iPmlEy, iPmlEz, iPmlHx, iPmlHy, iPmlHz):

    f_dtype = Parm['f_dtype']
    maxiter = Parm['solver'][0]
    nout    = Parm['solver'][1]

    nfeed  = iFeed.shape[0]
    npoint = iPoint.shape[0]

    Ex = np.zeros(NN, f_dtype)
    Ey = np.zeros(NN, f_dtype)
    Ez = np.zeros(NN, f_dtype)
    Hx = np.zeros(NN, f_dtype)
    Hy = np.zeros(NN, f_dtype)
    Hz = np.zeros(NN, f_dtype)

    Exy = Exz = Eyz = Eyx = Ezx = Ezy = None
    Hxy = Hxz = Hyz = Hyx = Hzx = Hzy = None

    if Parm['abc'][0] == 1:
        Exy = np.zeros(iPmlEx.shape[0], f_dtype)
        Exz = np.zeros(iPmlEx.shape[0], f_dtype)
        Eyz = np.zeros(iPmlEy.shape[0], f_dtype)
        Eyx = np.zeros(iPmlEy.shape[0], f_dtype)
        Ezx = np.zeros(iPmlEz.shape[0], f_dtype)
        Ezy = np.zeros(iPmlEz.shape[0], f_dtype)

        Hxy = np.zeros(iPmlHx.shape[0], f_dtype)
        Hxz = np.zeros(iPmlHx.shape[0], f_dtype)
        Hyz = np.zeros(iPmlHy.shape[0], f_dtype)
        Hyx = np.zeros(iPmlHy.shape[0], f_dtype)
        Hzx = np.zeros(iPmlHz.shape[0], f_dtype)
        Hzy = np.zeros(iPmlHz.shape[0], f_dtype)

    dsize = NN * len(Freq2)

    # 単精度複素数
    cEx = np.zeros(dsize, 'c8')
    cEy = np.zeros(dsize, 'c8')
    cEz = np.zeros(dsize, 'c8')
    cHx = np.zeros(dsize, 'c8')
    cHy = np.zeros(dsize, 'c8')
    cHz = np.zeros(dsize, 'c8')

    VFeed  = np.zeros((nfeed,  maxiter + 1), 'f8')
    IFeed  = np.zeros((nfeed,  maxiter + 1), 'f8')
    VPoint = np.zeros((npoint, maxiter + 1), 'f8')

    Eiter = np.zeros((maxiter + 1) // nout + 1, 'f8')
    Hiter = np.zeros((maxiter + 1) // nout + 1, 'f8')
    Iiter = np.zeros((maxiter + 1) // nout + 1, 'i4')

    return \
    Ex, Ey, Ez, Hx, Hy, Hz, \
    Exy, Exz, Eyz, Eyx, Ezx, Ezy, \
    Hxy, Hxz, Hyz, Hyx, Hzx, Hzy, \
    cEx, cEy, cEz, cHx, cHy, cHz, \
    VFeed, IFeed, VPoint, \
    Eiter, Hiter, Iiter
