# -*- coding: utf-8 -*-
"""
planewave.py (CUDA)
"""

import math
from numba import cuda

# 入射波電界とその時間微分(device関数)
@cuda.jit(device=True, cache=True, inline=True)
def f(x, y, z, t, f0, p):

    t -= ((x - p[ 9]) * p[6] \
        + (y - p[10]) * p[7] \
        + (z - p[11]) * p[8]) / p[14]

    at = p[12] * t
    ex = math.exp(-at**2) if at**2 < 16 else 0

    fi = at * ex * f0
    dfi = p[13] * p[12] * (1 - 2 * at**2) * ex * f0

    return fi, dfi
