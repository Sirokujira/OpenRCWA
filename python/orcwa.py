# -*- coding: utf-8 -*-
"""
OpenRCWA (Python)
Version 4.2.0
ost.py : solver
"""

import sys
import numpy as np
import numba
from mpi4py import MPI
import warnings
import sol.input_data, sol.input_setup, sol.setup, sol.monitor, sol.cputime
import sol.geometry, sol.plot3d_geom, sol.outputChars, sol.save_bin
import sol.solve
import sol.comm
import sol_cuda.solve

def main(argv):
    # ロゴ
    version = 'OpenRCWA (Python) Version 4.2.0'

    # 計算モード
    GPU = 0    # 0=CPU/1=GPU
    VECTOR = 0  # VECTOR: 0=OFF/1=ON

    # 入力データファイル名(orcwa.pyからの相対パス)
    orcwa_in = "python.orcwa"
    #orcwa_in = "../data/sample/1st_sample.orcwa"
    #orcwa_in = "../data/benchmark/benchmark100.orcwa"

    # Numbaスレッド数
    thread = 4

    # MPI領域分割数
    Npx = 1
    Npy = 1
    Npz = 1

    # 型宣言(適当に変更する)
    f_dtype = 'f4'  # 'f4' or 'f8'  (単精度/倍精度)
    i_dtype = 'u1'  # 'u1' or 'i4'  (属性数256以下/以上)

    # 終了時のプロンプト
    prompt = 0

    # 引数処理(引数があるときは引数優先)
    if len(argv) > 1:
        GPU, VECTOR, thread, Npx, Npy, Npz, f_dtype, prompt, orcwa_in \
        = _args(argv, version)
    #print(GPU, VECTOR, thread, Npx, Npy, Npz, f_dtype, prompt, orcwa_in)

    # GPU時のwarningを抑制
    if GPU:
        warnings.filterwarnings('ignore')

    # MPI
    comm = MPI.COMM_WORLD
    comm_size = comm.Get_size()  # 非MPI時は1
    comm_rank = comm.Get_rank()  # 非MPI時は0
    if comm_size == 1:
        Npx = Npy = Npz = 1  # 非MPI時は領域分割なし
    elif Npx * Npy * Npz != comm_size:
        Npx = comm_size  # 分割数が正しくないときはすべてX方向に分割する
        Npy = 1
        Npz = 1

    # io : ON/OFF
    io = (comm_rank == 0)

    # 出力ファイル名
    fn_log = 'orcwa.log'
    fn_out = 'orcwa.npz'

    # Numbaスレッド数設定
    numba.set_num_threads(thread)

    # cpu time
    cpu = [0] * 5
    cpu[0] = sol.cputime.t(comm_size, GPU)

    # logファイルを開く
    fp_log = None
    if io:
        fp_log = open(fn_log, 'wt', encoding='utf-8')

    # 経過表示 (1)
    if io:
        # ロゴ
        logo = '<<< %s >>>\n%s, process=%dx%dx%d=%d, thread=%d, vector=%s %s' % \
            (version, ('GPU' if GPU else 'CPU'), Npx, Npy, Npz, comm_size, \
            thread, ('on' if VECTOR else 'off'), ('(single)' if (f_dtype == 'f4') else '(double)'))
        sol.monitor.log1(fp_log, logo)

    # パラメーター (辞書型 : スカラーまたは固定サイズの配列)
    Parm = sol.input_data.const()

    # パラメーター追加
    Parm['vector'] = VECTOR
    Parm['f_dtype'] = f_dtype
    Parm['i_dtype'] = i_dtype
    Parm['comm_size'] = comm_size
    Parm['comm_rank'] = comm_rank
    Parm['Npx'] = Npx
    Parm['Npy'] = Npy
    Parm['Npz'] = Npz
    #print(comm_size, comm_rank, Parm)

    # [1] データ入力
    Nx = Ny = Nz = 0
    Xn = Yn = Zn = \
    iMaterial = fMaterial = iGeometry = fGeometry = \
    iFeed = fFeed = iPoint = fPoint = iInductor = fInductor = \
    Freq1 = Freq2 = None
    if io:
        (Nx, Ny, Nz, Xn, Yn, Zn,
        iMaterial, fMaterial, iGeometry, fGeometry, iFeed, fFeed,
        iPoint, fPoint, iInductor, fInductor, Freq1, Freq2) \
        = sol.input_data.read(orcwa_in, Parm)

    # 物体形状の線分を作成する (図形出力用)
    if io:
        gline, mline = sol.geometry.lines(iGeometry, fGeometry)

    # 物体形状を図形出力して終了する
    if io and (Parm['plot3dgeom'] == 1):
        sol.plot3d_geom.shape(Parm, Nx, Ny, Nz, Xn, Yn, Zn, iFeed, gline, mline)

    # (MPI) broadcast
    if Parm['comm_size'] > 1:
        (Nx, Ny, Nz, Xn, Yn, Zn,
        iMaterial, fMaterial, iGeometry, fGeometry,
        iFeed, fFeed, iPoint, fPoint, iInductor, fInductor, Freq1, Freq2) \
        = sol.comm.broadcast(Parm,
            Nx, Ny, Nz, Xn, Yn, Zn,
            iMaterial, fMaterial, iGeometry, fGeometry,
            iFeed, fFeed, iPoint, fPoint, iInductor, fInductor, Freq1, Freq2)

    # [2] 計算の準備作業

    # セル中心
    Xc, Yc, Zc = sol.input_setup.cell_center(Xn, Yn, Zn)

    # 各種準備
    sol.input_setup.setup(
        Parm, Xn, Yn, Zn, Xc, Yc, Zc,
        iMaterial, fMaterial, iGeometry, fGeometry, iFeed, fFeed,
        iPoint, fPoint, iInductor, fInductor, Freq1, Freq2)

    # 平面波入射データ
    fPlanewave = np.zeros(15, 'f8')
    if Parm['source'] == 1:
        sol.input_setup.planewave(Parm, Xn, Yn, Zn, Freq2, fPlanewave)

    # 電磁界配列計算用の係数
    iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0, NN, Ipx, Ipy, Ipz \
    = sol.setup.getIndex(Parm, Nx, Ny, Nz, Npx, Npy, Npz, Npx, Npy, Npz, Parm['comm_rank'])
    #print(comm_size, comm_rank, iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0, NN, Ipx, Ipy, Ipz)
    Parm['Ipx'] = Ipx
    Parm['Ipy'] = Ipy
    Parm['Ipz'] = Ipz

    # 各種準備
    (iEx, iEy, iEz, iHx, iHy, iHz,
    C1E, C2E, C1H, C2H,
    K1Ex, K2Ex, K1Ey, K2Ey, K1Ez, K2Ez,
    K1Hx, K2Hx, K1Hy, K2Hy, K1Hz, K2Hz,
    RXn, RYn, RZn, RXc, RYc, RZc,
    fMurHx, fMurHy, fMurHz, iMurHx, iMurHy, iMurHz,
    iPmlEx, iPmlEy, iPmlEz, iPmlHx, iPmlHy, iPmlHz,
    gPmlXn, gPmlYn, gPmlZn, gPmlXc, gPmlYc, gPmlZc, rPmlE, rPmlH, rPml,
    iDispEx, iDispEy, iDispEz, fDispEx, fDispEy, fDispEz,
    cEdft, cHdft, cFdft) \
    = sol.setup.setup(
        Parm, fPlanewave, Nx, Ny, Nz, Xn, Yn, Zn, Xc, Yc, Zc,
        iMaterial, fMaterial, iGeometry, fGeometry, Freq2,
        iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0, NN)

    # Yee電界点を図形出力して終了する（時間がかかる）
    if io and (Parm['plot3dgeom'] == 2):
        sol.plot3d_geom.cell(Parm, Nx, Ny, Nz, Xn, Yn, Zn,
            iEx, iEy, iEz, Ni, Nj, Nk, N0)

    # 経過表示 (2)
    if io:
        sol.monitor.log2(fp_log,
            GPU, Parm, Nx, Ny, Nz, NN,
            iMaterial, iGeometry, iFeed, iPoint, Freq1, Freq2,
            iPmlEx, iPmlEy, iPmlEz, iPmlHx, iPmlHy, iPmlHz)

    # [3] 計算の主要部

    cpu[1] = sol.cputime.t(comm_size, GPU)
    tdft = 0

    # RCWA反復計算ファイル名
    SOLVER = sol.solve.iteration if GPU == 0 else sol_cuda.solve.iteration

    # RCWA反復計算(RCWA計算の主要部)
    (cEx, cEy, cEz, cHx, cHy, cHz,
    VFeed, IFeed, VPoint, Ntime,
    Eiter, Hiter, Iiter, Niter,
    iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0, NN,
    tdft) \
    = SOLVER(
        io, fp_log,
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
        iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0, NN)

    cpu[3] = sol.cputime.t(comm_size, GPU)
    cpu[2] = cpu[3] - tdft

    # [4] 出力

    if io:
        # 計算結果の一部を計算しlogに出力する
        Zin, Ref, Pin, Spara, Coupling, fSurface, cSurface \
        = sol.outputChars.out(
            fp_log, Parm, Nx, Ny, Nz, Xn, Yn, Zn, Xc, Yc, Zc,
            iFeed, fFeed, iPoint, Freq1, Freq2,
            cEx, cEy, cEz, cHx, cHy, cHz,
            VFeed, IFeed, VPoint, Ntime, cFdft,
            Ni, Nj, Nk, N0, NN)

        # 経過表示 (3)
        sol.monitor.log3(fp_log, fn_log, fn_out)

        # 計算結果をファイルに出力する
        sol.save_bin.save(fn_out,
            Parm['title'], Parm['dt'], Parm['source'], fPlanewave, fFeed[:, 5],
            Ni, Nj, Nk, N0, NN,
            Xn, Yn, Zn, Freq1, Freq2,
            cEx, cEy, cEz, cHx, cHy, cHz,
            Eiter, Hiter, Iiter, Niter,
            VFeed, IFeed, VPoint, Ntime,
            Zin, Ref, Pin, Spara, Coupling,
            fSurface, cSurface, gline)

    # free
    cEx = cEy = cEz = cHx = cHy = cHz = None

    # cpu time
    cpu[4] = sol.cputime.t(comm_size, GPU)

    if io:
        # 経過表示 (4)
        sol.monitor.log4(fp_log, cpu)

        # logファイルを閉じる
        fp_log.close()

    # prompt
    if io and prompt:
        input()

# (private) 引数処理
def _args(argv, version):

    usage = 'Usage : python orcwa.py [-cpu|-gpu] [-n <thread>] [-p <x> <y> <z>] [-no-vector|-vector] [-single|-double] <datafile>'
    GPU = 0
    VECTOR = 0
    thread = 1
    Npx = Npy = Npz = 1
    f_dtype = 'f4'
    prompt = 0
    orcwa_in = ''

    i = 1
    while i < len(argv):
        arg = argv[i].lower()
        if   arg == '-gpu':
            GPU = 1
            i += 1
        elif arg == '-cpu':
            GPU = 0
            i += 1
        elif arg == '-n':
            thread = int(argv[i + 1])
            i += 2
        elif arg == '-p':
            Npx = int(argv[i + 1])
            Npy = int(argv[i + 2])
            Npx = int(argv[i + 1])
            Npy = int(argv[i + 2])
            Npz = int(argv[i + 3])
            i += 4
        elif arg == '-vector':
            VECTOR = 1
            i += 1
        elif arg == '-no-vector':
            VECTOR = 0
            i += 1
        elif arg == '-single':
            f_dtype = 'f4'
            i += 1
        elif arg == '-double':
            f_dtype = 'f8'
            i += 1
        elif arg == '-prompt':
            prompt = 1
            i += 1
        elif arg == '--help':
            print(usage)
            sys.exit()
        elif arg == '--version':
            print(version)
            sys.exit()
        else:
            orcwa_in = argv[i]
            i += 1

    return GPU, VECTOR, thread, Npx, Npy, Npz, f_dtype, prompt, orcwa_in

# enyry point
if __name__ == "__main__":
    main(sys.argv)
