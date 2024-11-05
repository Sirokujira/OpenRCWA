# -*- coding: utf-8 -*-
"""
save_bin.py
"""

import numpy as np

def save(fn,
    title, dt, source, fPlanewave, z0,
    Ni, Nj, Nk, N0, NN,
    Xn, Yn, Zn, Freq1, Freq2,
    cEx, cEy, cEz, cHx, cHy, cHz,
    Eiter, Hiter, Iiter, Niter,
    VFeed, IFeed, VPoint, Ntime,
    Zin, Ref, Pin, Spara, Coupling,
    fSurface, cSurface, gline):

    np.savez(fn,
    title, dt, source, fPlanewave, z0,
    Ni, Nj, Nk, N0, NN,
    Xn, Yn, Zn, Freq1, Freq2,
    cEx, cEy, cEz, cHx, cHy, cHz,
    Eiter, Hiter, Iiter, Niter,
    VFeed, IFeed, VPoint, Ntime,
    Zin, Ref, Pin, Spara, Coupling,
    fSurface, cSurface, gline)
