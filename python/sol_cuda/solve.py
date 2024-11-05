# -*- coding: utf-8 -*-
"""
solve.py (CUDA)
RCWA法の主要部:反復計算
"""

import numpy as np
from numba import cuda
import sol_cuda.updateEx, sol_cuda.updateEy, sol_cuda.updateEz
import sol_cuda.updateHx, sol_cuda.updateHy, sol_cuda.updateHz
import sol_cuda.MurH
import sol_cuda.PmlEx, sol_cuda.PmlEy, sol_cuda.PmlEz
import sol_cuda.PmlHx, sol_cuda.PmlHy, sol_cuda.PmlHz
import sol_cuda.PbcX, sol_cuda.PbcY, sol_cuda.PbcZ
import sol_cuda.DispEx, sol_cuda.DispEy, sol_cuda.DispEz
import sol_cuda.feed, sol_cuda.inductor, sol_cuda.pointV
import sol_cuda.average
import sol_cuda.Near3d
import sol_cuda.comm_X, sol_cuda.comm_Y, sol_cuda.comm_Z
import sol.comm, sol.monitor, sol.cputime

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

    # MPI用buffer作成
    SendBuf_Bx_hy = SendBuf_Bx_hz = RecvBuf_Bx_hy = RecvBuf_Bx_hz = \
    SendBuf_By_hz = SendBuf_By_hx = RecvBuf_By_hz = RecvBuf_By_hx = \
    SendBuf_Bz_hx = SendBuf_Bz_hy = RecvBuf_Bz_hx = RecvBuf_Bz_hy = None
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
    #print(Bx_jhy, Bx_jhz, Bx_khy, Bx_khz)

    # 計算に必要な配列を確保して初期化する
    Ex, Ey, Ez, Hx, Hy, Hz, \
    Exy, Exz, Eyz, Eyx, Ezx, Ezy, \
    Hxy, Hxz, Hyz, Hyx, Hzx, Hzy, \
    cEx, cEy, cEz, cHx, cHy, cHz, \
    VFeed, IFeed, VPoint, \
    Eiter, Hiter, Iiter \
    = _initfield(
        Parm, NN, iFeed, iPoint, Freq2,
        iPmlEx, iPmlEy, iPmlEz, iPmlHx, iPmlHy, iPmlHz)
    
    # host memory を　device memory　にコピーする(変数の頭に"d_"がつく)
    d_Xn, d_Yn, d_Zn, d_Xc, d_Yc, d_Zc, \
    d_RXn, d_RYn, d_RZn, d_RXc, d_RYc, d_RZc, \
    d_iEx, d_iEy, d_iEz, d_iHx, d_iHy, d_iHz, \
    d_C1E, d_C2E, d_C1H, d_C2H, \
    d_K1Ex, d_K2Ex, d_K1Ey, d_K2Ey, d_K1Ez, d_K2Ez, \
    d_K1Hx, d_K2Hx, d_K1Hy, d_K2Hy, d_K1Hz, d_K2Hz, \
    d_Ex, d_Ey, d_Ez, d_Hx, d_Hy, d_Hz, \
    d_Exy, d_Exz, d_Eyz, d_Eyx, d_Ezx, d_Ezy, \
    d_Hxy, d_Hxz, d_Hyz, d_Hyx, d_Hzx, d_Hzy, \
    d_cEx, d_cEy, d_cEz, d_cHx, d_cHy, d_cHz, \
    d_fMurHx, d_fMurHy, d_fMurHz, d_iMurHx, d_iMurHy, d_iMurHz, \
    d_iPmlEx, d_iPmlEy, d_iPmlEz, d_iPmlHx, d_iPmlHy, d_iPmlHz, \
    d_gPmlXn, d_gPmlYn, d_gPmlZn, d_gPmlXc, d_gPmlYc, d_gPmlZc, \
    d_rPmlE, d_rPmlH, d_rPml, \
    d_fPlanewave, d_iFeed, d_fFeed, d_iPoint, d_fPoint, d_iInductor, d_fInductor, \
    d_iDispEx, d_fDispEx, d_iDispEy, d_fDispEy, d_iDispEz, d_fDispEz, \
    d_VFeed, d_IFeed, d_VPoint, \
    d_SendBuf_Bx_hy, d_SendBuf_Bx_hz, d_RecvBuf_Bx_hy, d_RecvBuf_Bx_hz, \
    d_SendBuf_By_hz, d_SendBuf_By_hx, d_RecvBuf_By_hz, d_RecvBuf_By_hx, \
    d_SendBuf_Bz_hx, d_SendBuf_Bz_hy, d_RecvBuf_Bz_hx, d_RecvBuf_Bz_hy  \
    = _alloc_device_memory(
        Xn, Yn, Zn, Xc, Yc, Zc,
        RXn, RYn, RZn, RXc, RYc, RZc,
        iEx, iEy, iEz, iHx, iHy, iHz,
        C1E, C2E, C1H, C2H,
        K1Ex, K2Ex, K1Ey, K2Ey, K1Ez, K2Ez,
        K1Hx, K2Hx, K1Hy, K2Hy, K1Hz, K2Hz,
        Ex, Ey, Ez, Hx, Hy, Hz,
        Exy, Exz, Eyz, Eyx, Ezx, Ezy,
        Hxy, Hxz, Hyz, Hyx, Hzx, Hzy,
        cEx, cEy, cEz, cHx, cHy, cHz,
        fMurHx, fMurHy, fMurHz, iMurHx, iMurHy, iMurHz,
        iPmlEx, iPmlEy, iPmlEz, iPmlHx, iPmlHy, iPmlHz,
        gPmlXn, gPmlYn, gPmlZn, gPmlXc, gPmlYc, gPmlZc,
        rPmlE, rPmlH, rPml,
        fPlanewave, iFeed, fFeed, iPoint, fPoint, iInductor, fInductor,
        iDispEx, fDispEx, iDispEy, fDispEy, iDispEz, fDispEz,
        VFeed, IFeed, VPoint,
        SendBuf_Bx_hy, SendBuf_Bx_hz, RecvBuf_Bx_hy, RecvBuf_Bx_hz,
        SendBuf_By_hz, SendBuf_By_hx, RecvBuf_By_hz, RecvBuf_By_hx,
        SendBuf_Bz_hx, SendBuf_Bz_hy, RecvBuf_Bz_hx, RecvBuf_Bz_hy)

    # 平均電磁界用の device memory
    average_array = np.zeros((iMax - iMin) * (jMax - jMin) * (kMax - kMin), Parm['f_dtype'])
    d_average_array = cuda.to_device(average_array)

    # block size
    block3d = (32, 4, 1)
    block2d = (16, 16)
    block1d = 128

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
        sol_cuda.updateHx.calHx(block3d,
            Parm, t, d_fPlanewave, d_Xn, d_Yc, d_Zc, d_Hx, d_Ey, d_Ez, d_iHx,
            d_C1H, d_C2H, d_K1Hx, d_K2Hx, d_RYc, d_RZc,
            iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)
        sol_cuda.updateHy.calHy(block3d,
            Parm, t, d_fPlanewave, d_Yn, d_Zc, d_Xc, d_Hy, d_Ez, d_Ex, d_iHy,
            d_C1H, d_C2H, d_K1Hy, d_K2Hy, d_RZc, d_RXc,
            iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)
        sol_cuda.updateHz.calHz(block3d,
            Parm, t, d_fPlanewave, d_Zn, d_Xc, d_Yc, d_Hz, d_Ex, d_Ey, d_iHz,
            d_C1H, d_C2H, d_K1Hz, d_K2Hz, d_RXc, d_RYc,
            iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)

        # 磁界に吸収境界条件を適用する
        if   abc == 0:
            # Mur
            sol_cuda.MurH.calcH(block1d,
                d_Hx, d_fMurHx, d_iMurHx, Ni, Nj, Nk, N0)
            sol_cuda.MurH.calcH(block1d,
                d_Hy, d_fMurHy, d_iMurHy, Ni, Nj, Nk, N0)
            sol_cuda.MurH.calcH(block1d,
                d_Hz, d_fMurHz, d_iMurHz, Ni, Nj, Nk, N0)
        elif abc == 1:
            # PML
            sol_cuda.PmlHx.calHx(block1d,
                Ny, Nz, d_Hx, d_Ey, d_Ez, d_Hxy, d_Hxz, d_RYc, d_RZc,
                d_iPmlHx, d_gPmlYc, d_gPmlZc, d_rPmlH, d_rPml, pml_l,
                Ni, Nj, Nk, N0)
            sol_cuda.PmlHy.calHy(block1d,
                Nz, Nx, d_Hy, d_Ez, d_Ex, d_Hyz, d_Hyx, d_RZc, d_RXc,
                d_iPmlHy, d_gPmlZc, d_gPmlXc, d_rPmlH, d_rPml, pml_l,
                Ni, Nj, Nk, N0)
            sol_cuda.PmlHz.calHz(block1d,
                Nx, Ny, d_Hz, d_Ex, d_Ey, d_Hzx, d_Hzy, d_RXc, d_RYc,
                d_iPmlHz, d_gPmlXc, d_gPmlYc, d_rPmlH, d_rPml, pml_l,
                Ni, Nj, Nk, N0)

        # 磁界に周期境界条件を適用する(mode=1)
        if pbc[0] == 1:  # X
            if Npx > 1:  # MPI
                sol_cuda.comm_X.share(1, d_Hy, d_Hz,
                    SendBuf_Bx_hy, SendBuf_Bx_hz, RecvBuf_Bx_hy, RecvBuf_Bx_hz,
                    d_SendBuf_Bx_hy, d_SendBuf_Bx_hz, d_RecvBuf_Bx_hy, d_RecvBuf_Bx_hz,
                    Bx_jhy, Bx_jhz, Bx_khy, Bx_khz,
                    Npx, Npy, Npz, Ipx, Ipy, Ipz,
                    iMin, iMax, Ni, Nj, Nk, N0)
            else:
                sol_cuda.PbcX.x(block2d, d_Hy, d_Hz, Nx, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)
        if pbc[1] == 1:  # Y
            if Npy > 1:  # MPI
                sol_cuda.comm_Y.share(1, d_Hz, d_Hx,
                    SendBuf_By_hz, SendBuf_By_hx, RecvBuf_By_hz, RecvBuf_By_hx,
                    d_SendBuf_By_hz, d_SendBuf_By_hx, d_RecvBuf_By_hz, d_RecvBuf_By_hx,
                    By_khz, By_khx, By_ihz, By_ihx,
                    Npx, Npy, Npz, Ipx, Ipy, Ipz,
                    jMin, jMax, Ni, Nj, Nk, N0)
            else:
                sol_cuda.PbcY.y(block2d, d_Hz, d_Hx, Ny, kMin, kMax, iMin, iMax, Ni, Nj, Nk, N0)
        if pbc[2] == 1:  # Z
            if Npz > 1:  # MPI
                sol_cuda.comm_Z.share(1, d_Hx, d_Hy,
                    SendBuf_Bz_hx, SendBuf_Bz_hy, RecvBuf_Bz_hx, RecvBuf_Bz_hy,
                    d_SendBuf_Bz_hx, d_SendBuf_Bz_hy, d_RecvBuf_Bz_hx, d_RecvBuf_Bz_hy,
                    Bz_ihx, Bz_ihy, Bz_jhx, Bz_jhy,
                    Npx, Npy, Npz, Ipx, Ipy, Ipz,
                    kMin, kMax, Ni, Nj, Nk, N0)
            else:
                sol_cuda.PbcZ.z(block2d, d_Hx, d_Hy, Nz, iMin, iMax, jMin, jMax, Ni, Nj, Nk, N0)

        # MPI: 境界の磁界を共有する
        if Npx > 1:
            sol_cuda.comm_X.share(0, d_Hy, d_Hz,
                SendBuf_Bx_hy, SendBuf_Bx_hz, RecvBuf_Bx_hy, RecvBuf_Bx_hz,
                d_SendBuf_Bx_hy, d_SendBuf_Bx_hz, d_RecvBuf_Bx_hy, d_RecvBuf_Bx_hz,
                Bx_jhy, Bx_jhz, Bx_khy, Bx_khz,
                Npx, Npy, Npz, Ipx, Ipy, Ipz,
                iMin, iMax, Ni, Nj, Nk, N0)
        if Npy > 1:
            sol_cuda.comm_Y.share(0, d_Hz, d_Hx,
                SendBuf_By_hz, SendBuf_By_hx, RecvBuf_By_hz, RecvBuf_By_hx,
                d_SendBuf_By_hz, d_SendBuf_By_hx, d_RecvBuf_By_hz, d_RecvBuf_By_hx,
                By_khz, By_khx, By_ihz, By_ihx,
                Npx, Npy, Npz, Ipx, Ipy, Ipz,
                jMin, jMax, Ni, Nj, Nk, N0)
        if Npz > 1:
            sol_cuda.comm_Z.share(0, d_Hx, d_Hy,
                SendBuf_Bz_hx, SendBuf_Bz_hy, RecvBuf_Bz_hx, RecvBuf_Bz_hy,
                d_SendBuf_Bz_hx, d_SendBuf_Bz_hy, d_RecvBuf_Bz_hx, d_RecvBuf_Bz_hy,
                Bz_ihx, Bz_ihy, Bz_jhx, Bz_jhy,
                Npx, Npy, Npz, Ipx, Ipy, Ipz,
                kMin, kMax, Ni, Nj, Nk, N0)

        # 電界を更新する
        t += 0.5 * dt
        sol_cuda.updateEx.calEx(block3d,
            Parm, t, d_fPlanewave, d_Xc, d_Yn, d_Zn, d_Ex, d_Hy, d_Hz, d_iEx,
            d_C1E, d_C2E, d_K1Ex, d_K2Ex, d_RYn, d_RZn,
            iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)
        sol_cuda.updateEy.calEy(block3d,
            Parm, t, d_fPlanewave, d_Yc, d_Zn, d_Xn, d_Ey, d_Hz, d_Hx, d_iEy,
            d_C1E, d_C2E, d_K1Ey, d_K2Ey, d_RZn, d_RXn,
            iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)
        sol_cuda.updateEz.calEz(block3d,
            Parm, t, d_fPlanewave, d_Zc, d_Xn, d_Yn, d_Ez, d_Hx, d_Hy, d_iEz,
            d_C1E, d_C2E, d_K1Ez, d_K2Ez, d_RXn, d_RYn,
            iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)

        # 分散性媒質の電界を更新する
        if iDispEx.shape[0] > 0:
            sol_cuda.DispEx.calEx(block1d,
                source, t, d_fPlanewave, d_Xc, d_Yn, d_Zn, d_Ex, d_iDispEx, d_fDispEx,
                Ni, Nj, Nk, N0)
        if iDispEy.shape[0] > 0:
            sol_cuda.DispEy.calEy(block1d,
                source, t, d_fPlanewave, d_Yc, d_Zn, d_Xn, d_Ey, d_iDispEy, d_fDispEy,
                Ni, Nj, Nk, N0)
        if iDispEz.shape[0] > 0:
            sol_cuda.DispEz.calEz(block1d,
                source, t, d_fPlanewave, d_Zc, d_Xn, d_Yn, d_Ez, d_iDispEz, d_fDispEz,
                Ni, Nj, Nk, N0)

        # 電界に吸収境界条件を適用する
        if abc == 1:
            # PML
            sol_cuda.PmlEx.calEx(block1d,
                Ny, Nz, d_Ex, d_Hy, d_Hz, d_Exy, d_Exz, d_RYn, d_RZn,
                d_iPmlEx, d_gPmlYn, d_gPmlZn, d_rPmlE, d_rPml, pml_l,
                Ni, Nj, Nk, N0)
            sol_cuda.PmlEy.calEy(block1d,
                Nz, Nx, d_Ey, d_Hz, d_Hx, d_Eyz, d_Eyx, d_RZn, d_RXn,
                d_iPmlEy, d_gPmlZn, d_gPmlXn, d_rPmlE, d_rPml, pml_l,
                Ni, Nj, Nk, N0)
            sol_cuda.PmlEz.calEz(block1d,
                Nx, Ny, d_Ez, d_Hx, d_Hy, d_Ezx, d_Ezy, d_RXn, d_RYn,
                d_iPmlEz, d_gPmlXn, d_gPmlYn, d_rPmlE, d_rPml, pml_l,
                Ni, Nj, Nk, N0)

        # 給電点に電圧を印加し、V/I時間波形を保存する
        if iFeed.shape[0] > 0:
            sol_cuda.feed.evi(block1d,
                itime, Parm, d_iFeed, d_fFeed, d_VFeed, d_IFeed,
                d_Ex, d_Ey, d_Ez, d_Hx, d_Hy, d_Hz, d_iEx, d_iEy, d_iEz,
                iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)

        # inductorを計算する
        if iInductor.shape[0] > 0:
            sol_cuda.inductor.calcL(block1d,
                Parm, d_iInductor, d_fInductor,
                d_Ex, d_Ey, d_Ez, d_Hx, d_Hy, d_Hz,
                iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)

        # 観測点の電圧波形を保存する
        if iPoint.shape[0] > 0:
            sol_cuda.pointV.v(block1d,
                itime, Parm, d_fPlanewave, d_Xn, d_Yn, d_Zn, d_Xc, d_Yc, d_Zc,
                d_iPoint, d_fPoint, d_VPoint,
                d_Ex, d_Ey, d_Ez,
                iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)

        # 全領域の電磁界のDFTを更新する(計算時間がかかる)
        t0 = sol.cputime.t(Parm['comm_size'], 1)
        sol_cuda.Near3d.dft(block3d,
            itime, Freq2, cEdft, cHdft,
            d_Ex, d_Ey, d_Ez, d_Hx, d_Hy, d_Hz,
            d_cEx, d_cEy, d_cEz, d_cHx, d_cHy, d_cHz,
            iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0, NN)
        tdft += sol.cputime.t(Parm['comm_size'], 1) - t0

        # 平均電磁界を計算して収束判定する
        if (itime % nout == 0) or (itime == maxiter):
            # 電磁界の和を計算する
            esum, hsum = sol_cuda.average.calcA(block3d, d_average_array,
                d_Ex, d_Ey, d_Ez, d_Hx, d_Hy, d_Hz,
                iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)

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

    # device から host にコピーする
    _copy_to_host(
        d_cEx, d_cEy, d_cEz, d_cHx, d_cHy, d_cHz,
        d_VFeed, d_IFeed, d_VPoint,
        cEx, cEy, cEz, cHx, cHy, cHz,
        VFeed, IFeed, VPoint)

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

        # 全領域の電磁界をrootに集める(iMin等は変更されている)
        cEx, cEy, cEz, cHx, cHy, cHz, \
        iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0, NN \
        = sol.comm.gather_near3d(Parm, Freq2,
            cEx, cEy, cEz, cHx, cHy, cHz,
            Nx, Ny, Nz, iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0, NN)

    # free
    Ex = Ey = Ez = Hx = Hy = Hz = \
    Exy = Exz = Eyz = Eyx = Ezx = Ezy = \
    Hxy = Hxz = Hyz = Hyx = Hzx = Hzy = \
    iEx = iEy = iEz = iHx = iHy = iHz = \
    K1Ex = K2Ex = K1Ey = K2Ey = K1Ez = K2Ez = \
    K1Hx = K2Hx = K1Hy = K2Hy = K1Hz = K2Hz = None

    # MPIではiMin等は変更されているので返す必要がある

    return \
    cEx, cEy, cEz, cHx, cHy, cHz, \
    VFeed, IFeed, VPoint, Ntime, \
    Eiter, Hiter, Iiter, Niter, \
    iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0, NN, \
    tdft

# (private) 計算に必要な配列を確保して初期化する
def _initfield(
    Parm, NN, iFeed, iPoint, Freq2,
    iPmlEx, iPmlEy, iPmlEz, iPmlHx, iPmlHy, iPmlHz):

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

# (private) host memory を device memory にコピーする
def _alloc_device_memory(
    Xn, Yn, Zn, Xc, Yc, Zc,
    RXn, RYn, RZn, RXc, RYc, RZc,
    iEx, iEy, iEz, iHx, iHy, iHz,
    C1E, C2E, C1H, C2H,
    K1Ex, K2Ex, K1Ey, K2Ey, K1Ez, K2Ez,
    K1Hx, K2Hx, K1Hy, K2Hy, K1Hz, K2Hz,
    Ex, Ey, Ez, Hx, Hy, Hz,
    Exy, Exz, Eyz, Eyx, Ezx, Ezy,
    Hxy, Hxz, Hyz, Hyx, Hzx, Hzy,
    cEx, cEy, cEz, cHx, cHy, cHz,
    fMurHx, fMurHy, fMurHz, iMurHx, iMurHy, iMurHz,
    iPmlEx, iPmlEy, iPmlEz, iPmlHx, iPmlHy, iPmlHz,
    gPmlXn, gPmlYn, gPmlZn, gPmlXc, gPmlYc, gPmlZc,
    rPmlE, rPmlH, rPml,
    fPlanewave, iFeed, fFeed, iPoint, fPoint, iInductor, fInductor,
    iDispEx, fDispEx, iDispEy, fDispEy, iDispEz, fDispEz,
    VFeed, IFeed, VPoint,
    SendBuf_Bx_hy, SendBuf_Bx_hz, RecvBuf_Bx_hy, RecvBuf_Bx_hz,
    SendBuf_By_hz, SendBuf_By_hx, RecvBuf_By_hz, RecvBuf_By_hx,
    SendBuf_Bz_hx, SendBuf_Bz_hy, RecvBuf_Bz_hx, RecvBuf_Bz_hy):

    d_Xn = cuda.to_device(Xn)
    d_Yn = cuda.to_device(Yn)
    d_Zn = cuda.to_device(Zn)
    d_Xc = cuda.to_device(Xc)
    d_Yc = cuda.to_device(Yc)
    d_Zc = cuda.to_device(Zc)

    d_RXn = cuda.to_device(RXn)
    d_RYn = cuda.to_device(RYn)
    d_RZn = cuda.to_device(RZn)
    d_RXc = cuda.to_device(RXc)
    d_RYc = cuda.to_device(RYc)
    d_RZc = cuda.to_device(RZc)

    d_iEx = cuda.to_device(iEx)
    d_iEy = cuda.to_device(iEy)
    d_iEz = cuda.to_device(iEz)
    d_iHx = cuda.to_device(iHx)
    d_iHy = cuda.to_device(iHy)
    d_iHz = cuda.to_device(iHz)

    d_C1E = cuda.to_device(C1E)
    d_C2E = cuda.to_device(C2E)
    d_C1H = cuda.to_device(C1H)
    d_C2H = cuda.to_device(C2H)

    d_K1Ex = cuda.to_device(K1Ex)
    d_K2Ex = cuda.to_device(K2Ex)
    d_K1Ey = cuda.to_device(K1Ey)
    d_K2Ey = cuda.to_device(K2Ey)
    d_K1Ez = cuda.to_device(K1Ez)
    d_K2Ez = cuda.to_device(K2Ez)
    d_K1Hx = cuda.to_device(K1Hx)
    d_K2Hx = cuda.to_device(K2Hx)
    d_K1Hy = cuda.to_device(K1Hy)
    d_K2Hy = cuda.to_device(K2Hy)
    d_K1Hz = cuda.to_device(K1Hz)
    d_K2Hz = cuda.to_device(K2Hz)

    d_Ex = cuda.to_device(Ex)
    d_Ey = cuda.to_device(Ey)
    d_Ez = cuda.to_device(Ez)
    d_Hx = cuda.to_device(Hx)
    d_Hy = cuda.to_device(Hy)
    d_Hz = cuda.to_device(Hz)

    d_Exy = cuda.to_device(Exy)
    d_Exz = cuda.to_device(Exz)
    d_Eyz = cuda.to_device(Eyz)
    d_Eyx = cuda.to_device(Eyx)
    d_Ezx = cuda.to_device(Ezx)
    d_Ezy = cuda.to_device(Ezy)

    d_Hxy = cuda.to_device(Hxy)
    d_Hxz = cuda.to_device(Hxz)
    d_Hyz = cuda.to_device(Hyz)
    d_Hyx = cuda.to_device(Hyx)
    d_Hzx = cuda.to_device(Hzx)
    d_Hzy = cuda.to_device(Hzy)

    d_cEx = cuda.to_device(cEx)
    d_cEy = cuda.to_device(cEy)
    d_cEz = cuda.to_device(cEz)
    d_cHx = cuda.to_device(cHx)
    d_cHy = cuda.to_device(cHy)
    d_cHz = cuda.to_device(cHz)

    d_fMurHx = cuda.to_device(fMurHx)
    d_fMurHy = cuda.to_device(fMurHy)
    d_fMurHz = cuda.to_device(fMurHz)
    d_iMurHx = cuda.to_device(iMurHx)
    d_iMurHy = cuda.to_device(iMurHy)
    d_iMurHz = cuda.to_device(iMurHz)

    d_iPmlEx = cuda.to_device(iPmlEx)
    d_iPmlEy = cuda.to_device(iPmlEy)
    d_iPmlEz = cuda.to_device(iPmlEz)
    d_iPmlHx = cuda.to_device(iPmlHx)
    d_iPmlHy = cuda.to_device(iPmlHy)
    d_iPmlHz = cuda.to_device(iPmlHz)
    d_gPmlXn = cuda.to_device(gPmlXn)
    d_gPmlYn = cuda.to_device(gPmlYn)
    d_gPmlZn = cuda.to_device(gPmlZn)
    d_gPmlXc = cuda.to_device(gPmlXc)
    d_gPmlYc = cuda.to_device(gPmlYc)
    d_gPmlZc = cuda.to_device(gPmlZc)
    d_rPmlE  = cuda.to_device(rPmlE)
    d_rPmlH  = cuda.to_device(rPmlH)
    d_rPml   = cuda.to_device(rPml)

    d_fPlanewave = cuda.to_device(fPlanewave)

    d_iFeed     = cuda.to_device(iFeed)
    d_fFeed     = cuda.to_device(fFeed)
    d_iPoint    = cuda.to_device(iPoint)
    d_fPoint    = cuda.to_device(fPoint)
    d_iInductor = cuda.to_device(iInductor)
    d_fInductor = cuda.to_device(fInductor)

    d_iDispEx = cuda.to_device(iDispEx)
    d_fDispEx = cuda.to_device(fDispEx)
    d_iDispEy = cuda.to_device(iDispEy)
    d_fDispEy = cuda.to_device(fDispEy)
    d_iDispEz = cuda.to_device(iDispEz)
    d_fDispEz = cuda.to_device(fDispEz)

    d_VFeed  = cuda.to_device(VFeed)
    d_IFeed  = cuda.to_device(IFeed)
    d_VPoint = cuda.to_device(VPoint)

    d_SendBuf_Bx_hy = cuda.to_device(SendBuf_Bx_hy)
    d_SendBuf_Bx_hz = cuda.to_device(SendBuf_Bx_hz)
    d_RecvBuf_Bx_hy = cuda.to_device(RecvBuf_Bx_hy)
    d_RecvBuf_Bx_hz = cuda.to_device(RecvBuf_Bx_hz)
    d_SendBuf_By_hz = cuda.to_device(SendBuf_By_hz)
    d_SendBuf_By_hx = cuda.to_device(SendBuf_By_hx)
    d_RecvBuf_By_hz = cuda.to_device(RecvBuf_By_hz)
    d_RecvBuf_By_hx = cuda.to_device(RecvBuf_By_hx)
    d_SendBuf_Bz_hx = cuda.to_device(SendBuf_Bz_hx)
    d_SendBuf_Bz_hy = cuda.to_device(SendBuf_Bz_hy)
    d_RecvBuf_Bz_hx = cuda.to_device(RecvBuf_Bz_hx)
    d_RecvBuf_Bz_hy = cuda.to_device(RecvBuf_Bz_hy)

    return \
    d_Xn, d_Yn, d_Zn, d_Xc, d_Yc, d_Zc, \
    d_RXn, d_RYn, d_RZn, d_RXc, d_RYc, d_RZc, \
    d_iEx, d_iEy, d_iEz, d_iHx, d_iHy, d_iHz, \
    d_C1E, d_C2E, d_C1H, d_C2H, \
    d_K1Ex, d_K2Ex, d_K1Ey, d_K2Ey, d_K1Ez, d_K2Ez, \
    d_K1Hx, d_K2Hx, d_K1Hy, d_K2Hy, d_K1Hz, d_K2Hz, \
    d_Ex, d_Ey, d_Ez, d_Hx, d_Hy, d_Hz, \
    d_Exy, d_Exz, d_Eyz, d_Eyx, d_Ezx, d_Ezy, \
    d_Hxy, d_Hxz, d_Hyz, d_Hyx, d_Hzx, d_Hzy, \
    d_cEx, d_cEy, d_cEz, d_cHx, d_cHy, d_cHz, \
    d_fMurHx, d_fMurHy, d_fMurHz, d_iMurHx, d_iMurHy, d_iMurHz, \
    d_iPmlEx, d_iPmlEy, d_iPmlEz, d_iPmlHx, d_iPmlHy, d_iPmlHz, \
    d_gPmlXn, d_gPmlYn, d_gPmlZn, d_gPmlXc, d_gPmlYc, d_gPmlZc, \
    d_rPmlE, d_rPmlH, d_rPml, \
    d_fPlanewave, d_iFeed, d_fFeed, d_iPoint, d_fPoint, d_iInductor, d_fInductor, \
    d_iDispEx, d_fDispEx, d_iDispEy, d_fDispEy, d_iDispEz, d_fDispEz, \
    d_VFeed, d_IFeed, d_VPoint, \
    d_SendBuf_Bx_hy, d_SendBuf_Bx_hz, d_RecvBuf_Bx_hy, d_RecvBuf_Bx_hz, \
    d_SendBuf_By_hz, d_SendBuf_By_hx, d_RecvBuf_By_hz, d_RecvBuf_By_hx, \
    d_SendBuf_Bz_hx, d_SendBuf_Bz_hy, d_RecvBuf_Bz_hx, d_RecvBuf_Bz_hy

# (private) device memory を host memory にコピーする
def _copy_to_host(
    d_cEx, d_cEy, d_cEz, d_cHx, d_cHy, d_cHz,
    d_VFeed, d_IFeed, d_VPoint,
    cEx, cEy, cEz, cHx, cHy, cHz,
    VFeed, IFeed, VPoint):

    d_cEx.copy_to_host(cEx)
    d_cEy.copy_to_host(cEy)
    d_cEz.copy_to_host(cEz)
    d_cHx.copy_to_host(cHx)
    d_cHy.copy_to_host(cHy)
    d_cHz.copy_to_host(cHz)

    d_VFeed.copy_to_host(VFeed)
    d_IFeed.copy_to_host(IFeed)
    d_VPoint.copy_to_host(VPoint)
