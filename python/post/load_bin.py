# -*- coding: utf-8 -*-
"""
load_bin.py
"""

import numpy as np

def load(fn):
    d = np.load(fn)
    title      = d[d.files[ 0]]
    dt         = d[d.files[ 1]]
    source     = d[d.files[ 2]]
    fPlanewave = d[d.files[ 3]]
    z0         = d[d.files[ 4]]
    Ni         = d[d.files[ 5]]
    Nj         = d[d.files[ 6]]
    Nk         = d[d.files[ 7]]
    N0         = d[d.files[ 8]]
    NN         = d[d.files[ 9]]
    Xn         = d[d.files[10]]
    Yn         = d[d.files[11]]
    Zn         = d[d.files[12]]
    Freq1      = d[d.files[13]]
    Freq2      = d[d.files[14]]
    cEx        = d[d.files[15]]
    cEy        = d[d.files[16]]
    cEz        = d[d.files[17]]
    cHx        = d[d.files[18]]
    cHy        = d[d.files[19]]
    cHz        = d[d.files[20]]
    Eiter      = d[d.files[21]]
    Hiter      = d[d.files[22]]
    Iiter      = d[d.files[23]]
    Niter      = d[d.files[24]]
    VFeed      = d[d.files[25]]
    IFeed      = d[d.files[26]]
    VPoint     = d[d.files[27]]
    Ntime      = d[d.files[28]]
    Zin        = d[d.files[29]]
    Ref        = d[d.files[30]]
    Pin        = d[d.files[31]]
    Spara      = d[d.files[32]]
    Coupling   = d[d.files[33]]
    fSurface   = d[d.files[34]]
    cSurface   = d[d.files[35]]
    gline      = d[d.files[36]]

    return \
    title, dt, source, fPlanewave, z0, \
    Ni, Nj, Nk, N0, NN, \
    Xn, Yn, Zn, Freq1, Freq2, \
    cEx, cEy, cEz, cHx, cHy, cHz, \
    Eiter, Hiter, Iiter, Niter, \
    VFeed, IFeed, VPoint, Ntime, \
    Zin, Ref, Pin, Spara, Coupling, \
    fSurface, cSurface, gline
