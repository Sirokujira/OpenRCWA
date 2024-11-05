# -*- coding: utf-8 -*-
"""
cputime.py
"""

import time
from mpi4py import MPI
from numba import cuda

def t(comm_size, gpu):

    if (comm_size > 1) and gpu:
        cuda.synchronize()
        MPI.COMM_WORLD.Barrier()
        t = MPI.Wtime()
    elif (comm_size > 1):
        MPI.COMM_WORLD.Barrier()
        t = MPI.Wtime()
    elif gpu:
        cuda.synchronize()
        t = time.time()
    else:
        t = time.time()

    return t
