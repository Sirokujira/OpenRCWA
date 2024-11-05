# -*- coding: utf-8 -*-
"""
monitor.py
"""

import sys
import datetime

def log1(fp, msg):
    log1_(fp,         msg)
    log1_(sys.stdout, msg)

def log2(fp,
        GPU, Parm, Nx, Ny, Nz, NN,
        iMaterial, iGeometry, iFeed, iPoint, Freq1, Freq2,
        iPmlEx, iPmlEy, iPmlEz, iPmlHx, iPmlHy, iPmlHz):
    log2_(fp,
        GPU, Parm, Nx, Ny, Nz, NN,
        iMaterial, iGeometry, iFeed, iPoint, Freq1, Freq2,
        iPmlEx, iPmlEy, iPmlEz, iPmlHx, iPmlHy, iPmlHz)
    log2_(sys.stdout,
        GPU, Parm, Nx, Ny, Nz, NN,
        iMaterial, iGeometry, iFeed, iPoint, Freq1, Freq2,
        iPmlEx, iPmlEy, iPmlEz, iPmlHx, iPmlHy, iPmlHz)

def log3(fp, fn_log, fn_out):
    log3_(fp,         fn_log, fn_out)
    log3_(sys.stdout, fn_log, fn_out)

def log4(fp, cpu):
    log4_(fp,         cpu) 
    log4_(sys.stdout, cpu)

# (private) message
def log1_(fp, msg):
    fp.write('%s\n' % msg)
    fp.flush()

# (private) condition
def log2_(fp,
    GPU, Parm, Nx, Ny, Nz, NN,
    iMaterial, iGeometry, iFeed, iPoint, Freq1, Freq2,
    iPmlEx, iPmlEy, iPmlEz, iPmlHx, iPmlHy, iPmlHz):

    cpu_mem, gpu_mem = _memory_size(Parm, Nx, Ny, Nz, Freq2, NN,
        iPmlEx, iPmlEy, iPmlEz, iPmlHx, iPmlHy, iPmlHz)
    out = _output_size(NN, Freq2)

    fp.write("%s\n" % datetime.datetime.now().ctime())
    fp.write("Title = %s\n" % Parm['title'])
    fp.write("Source = %s\n" % ("feed" if Parm['source'] == 0 else "plane wave"))
    fp.write("Cells = %d x %d x %d = %d\n" % (Nx, Ny, Nz, Nx * Ny * Nz))
    fp.write("No. of Materials  = %d\n" % len(iMaterial))
    fp.write("No. of Geometries = %d\n" % iGeometry.shape[0])
    if Parm['source'] == 0:
        fp.write("No. of Feeds      = %d\n" % iFeed.shape[0])
    fp.write("No. of Points     = %d\n" % iPoint.shape[0])
    fp.write("No. of Freq.s (1) = %d\n" % len(Freq1))
    fp.write("No. of Freq.s (2) = %d\n" % len(Freq2))
    fp.write("CPU Memory size   = %d [MB]\n" % cpu_mem)
    if GPU:
        fp.write("GPU Memory size   = %d [MB]\n" % gpu_mem)
    fp.write("Output filesize   = %d [MB]\n" % out)
    if   Parm['abc'][0] == 0:
        fp.write("ABC = Mur-1st\n")
    elif Parm['abc'][0] == 1:
        fp.write("ABC = PML (L=%d, M=%.2f, R0=%.2e)\n" % (Parm['abc'][1], Parm['abc'][2], Parm['abc'][3]) )
    if (sum(Parm['pbc']) > 0):
        fp.write("PBC : %s%s%s\n" % (("X" if Parm['pbc'][0] == 1 else ""),
                                     ("Y" if Parm['pbc'][1] == 1 else ""),
                                     ("Z" if Parm['pbc'][2] == 1 else "")))
    fp.write("Dt[sec] = %.4e, Tw[sec] = %.4e, Tw/Dt = %.3f\n" % (Parm['dt'], Parm['tw'], Parm['tw'] / Parm['dt']))
    fp.write("Iterations = %d, Convergence = %.3e\n" % (Parm['solver'][0], Parm['solver'][2]))
    fp.write("=== iteration start ===\n")
    fp.write("   step   <E>      <H>\n")
    fp.flush()

# (private) output files
def log3_(fp, fn_log, fn_out):
    fp.write("%s\n" % "=== output files ===")
    fp.write("%s, %s\n" % (fn_log, fn_out))
    fp.flush()

# (private) cpu time
def log4_(fp, cpu):
    fp.write("=== cpu time [sec] ===\n")
    fp.write("  part-1 : %11.3f\n" % (cpu[1] - cpu[0]))
    fp.write("  part-2 : %11.3f\n" % (cpu[2] - cpu[1]))
    fp.write("  part-3 : %11.3f\n" % (cpu[3] - cpu[2]))
    fp.write("  part-4 : %11.3f\n" % (cpu[4] - cpu[3]))
    fp.write("  --------------------\n")
    fp.write("  total  : %11.3f\n" % (cpu[4] - cpu[0]))
    fp.write("=== normal end ===\n")
    fp.write("%s\n" % datetime.datetime.now().ctime())
    fp.flush()

# (private) メモリーサイズ [MB]
def _memory_size(Parm, Nx, Ny, Nz, Freq2, NN,
    iPmlEx, iPmlEy, iPmlEz, iPmlHx, iPmlHy, iPmlHz):

    f_dsize = 4 if Parm['f_dtype'] == 'f4' else 8
    i_dsize = 1 if Parm['i_dtype'] == 'u1' else 4

    mem  = NN * 6 * f_dsize         # Ex, Ey, Ez, Hx, Hy, Hz
    mem += NN * 6 * i_dsize         # iEx, iEy, iEz, iHx, iHy, iHz

    if Parm['vector']:
        mem += NN * 12 * f_dsize    # K1Ex, K2Ex, ... K1Hx, K2Hx, ...

    # PML
    if Parm['abc'][0] == 1:
        num = iPmlEx.shape[0] + iPmlEy.shape[0] + iPmlEz.shape[0] \
            + iPmlHx.shape[0] + iPmlHy.shape[0] + iPmlHz.shape[0]
        mem += num * ((4 * 4) + (2 * 6 * f_dsize))  # i, j, k, m, Exy, Exz, ...

    cpu_mem = mem + len(Freq2) * Nx * Ny * Nz * 6 * 8  # 'c8': cEx, cEy, cEz, cHx, cHy, cHz

    gpu_mem =     mem / (1024 * 1024) + 1
    cpu_mem = cpu_mem / (1024 * 1024) + 1

    return cpu_mem, gpu_mem

# (private) 出力ファイルサイズ [MB]
def _output_size(NN, Freq2):

    mem = 0
    mem += len(Freq2) * NN * 6 * 8  # 'c8': cEx, cEy, cEz, cHx, cHy, cHz

    return mem // (1024 * 1024) + 1
