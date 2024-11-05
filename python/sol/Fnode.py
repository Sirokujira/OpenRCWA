# -*- coding: utf-8 -*-
"""
Fnode.py
"""

from numba import jit

# 指定した周波数と節点(i,j,k)の電界(複素数)
@jit(cache=True, nopython=True)
def e(ifreq, i, j, k,
    Nx, Ny, Nz, cEx, cEy, cEz, Ni, Nj, Nk, N0, NN):

    n0 = ifreq * NN

    # Ex
    if   i <= 0:
        i = 0
        n1 = n0 + (Ni * (i    )) + (Nj * (j    )) + (Nk * (k    )) + N0
        n2 = n0 + (Ni * (i + 1)) + (Nj * (j    )) + (Nk * (k    )) + N0
        #n1 = n0 + NEX(i,     j,     k,     Nx, Ny, Nz)
        #n2 = n0 + NEX(i + 1, j,     k,     Nx, Ny, Nz)
        cex = (cEx[n1] * 3 - cEx[n2] * 1) / 2
    elif i >= Nx:
        i = Nx
        n1 = n0 + (Ni * (i - 1)) + (Nj * (j    )) + (Nk * (k    )) + N0
        n2 = n0 + (Ni * (i - 2)) + (Nj * (j    )) + (Nk * (k    )) + N0
        #n1 = n0 + NEX(i - 1, j,     k,     Nx, Ny, Nz)
        #n2 = n0 + NEX(i - 2, j,     k,     Nx, Ny, Nz)
        cex = (cEx[n1] * 3 - cEx[n2] * 1) / 2
    else:
        n1 = n0 + (Ni * (i    )) + (Nj * (j    )) + (Nk * (k    )) + N0
        n2 = n0 + (Ni * (i - 1)) + (Nj * (j    )) + (Nk * (k    )) + N0
        #n1 = n0 + NEX(i,     j,     k,     Nx, Ny, Nz)
        #n2 = n0 + NEX(i - 1, j,     k,     Nx, Ny, Nz)
        cex = (cEx[n1] + cEx[n2]) / 2

    # Ey
    if   j <= 0:
        j = 0
        n1 = n0 + (Ni * (i    )) + (Nj * (j    )) + (Nk * (k    )) + N0
        n2 = n0 + (Ni * (i    )) + (Nj * (j + 1)) + (Nk * (k    )) + N0
        #n1 = n0 + NEY(i,     j,     k,     Nx, Ny, Nz)
        #n2 = n0 + NEY(i,     j + 1, k,     Nx, Ny, Nz )
        cey = (cEy[n1] * 3 - cEy[n2] * 1) / 2
    elif j >= Ny:
        j = Ny
        n1 = n0 + (Ni * (i    )) + (Nj * (j - 1)) + (Nk * (k    )) + N0
        n2 = n0 + (Ni * (i    )) + (Nj * (j - 2)) + (Nk * (k    )) + N0
        #n1 = n0 + NEY(i,     j - 1, k,     Nx, Ny, Nz)
        #n2 = n0 + NEY(i,     j - 2, k,     Nx, Ny, Nz)
        cey = (cEy[n1] * 3 - cEy[n2] * 1) / 2
    else:
        n1 = n0 + (Ni * (i    )) + (Nj * (j    )) + (Nk * (k    )) + N0
        n2 = n0 + (Ni * (i    )) + (Nj * (j - 1)) + (Nk * (k    )) + N0
        #n1 = n0 + NEY(i,     j,     k,     Nx, Ny, Nz)
        #n2 = n0 + NEY(i,     j - 1, k,     Nx, Ny, Nz)
        cey = (cEy[n1] + cEy[n2]) / 2

    # Ez
    if   k <= 0:
        k = 0
        n1 = n0 + (Ni * (i    )) + (Nj * (j    )) + (Nk * (k    )) + N0
        n2 = n0 + (Ni * (i    )) + (Nj * (j    )) + (Nk * (k + 1)) + N0
        #n1 = n0 + NEZ(i,     j,     k,     Nx, Ny, Nz)
        #n2 = n0 + NEZ(i,     j,     k + 1, Nx, Ny, Nz)
        cez = (cEz[n1] * 3 - cEz[n2] * 1) / 2
    elif k >= Nz:
        k = Nz
        n1 = n0 + (Ni * (i    )) + (Nj * (j    )) + (Nk * (k - 1)) + N0
        n2 = n0 + (Ni * (i    )) + (Nj * (j    )) + (Nk * (k - 2)) + N0
        #n1 = n0 + NEZ(i,     j,     k - 1, Nx, Ny, Nz)
        #n2 = n0 + NEZ(i,     j,     k - 2, Nx, Ny, Nz)
        cez = (cEz[n1] * 3 - cEz[n2] * 1) / 2
    else:
        n1 = n0 + (Ni * (i    )) + (Nj * (j    )) + (Nk * (k    )) + N0
        n2 = n0 + (Ni * (i    )) + (Nj * (j    )) + (Nk * (k - 1)) + N0
        #n1 = n0 + NEZ(i,     j,     k,     Nx, Ny, Nz)
        #n2 = n0 + NEZ(i,     j,     k - 1, Nx, Ny, Nz)
        cez = (cEz[n1] + cEz[n2]) / 2
    
    return cex, cey, cez

# 指定した周波数と節点(i,j,k)の磁界(複素数)
@jit(cache=True)
def h(ifreq, i, j, k,
    cHx, cHy, cHz, Ni, Nj, Nk, N0, NN):

    n0 = ifreq * NN

    # Hx
    n00 = n0 + (Ni * (i    )) + (Nj * (j    )) + (Nk * (k    )) + N0
    n01 = n0 + (Ni * (i    )) + (Nj * (j    )) + (Nk * (k - 1)) + N0
    n10 = n0 + (Ni * (i    )) + (Nj * (j - 1)) + (Nk * (k    )) + N0
    n11 = n0 + (Ni * (i    )) + (Nj * (j - 1)) + (Nk * (k - 1)) + N0
    #n00 = n0 + NHX(i,     j,     k,     Nx, Ny, Nz)
    #n01 = n0 + NHX(i,     j,     k - 1, Nx, Ny, Nz)
    #n10 = n0 + NHX(i,     j - 1, k,     Nx, Ny, Nz)
    #n11 = n0 + NHX(i,     j - 1, k - 1, Nx, Ny, Nz)
    chx = (cHx[n00] + cHx[n01] + cHx[n10] + cHx[n11]) / 4

    # Hy
    n00 = n0 + (Ni * (i    )) + (Nj * (j    )) + (Nk * (k    )) + N0
    n01 = n0 + (Ni * (i - 1)) + (Nj * (j    )) + (Nk * (k    )) + N0
    n10 = n0 + (Ni * (i    )) + (Nj * (j    )) + (Nk * (k - 1)) + N0
    n11 = n0 + (Ni * (i - 1)) + (Nj * (j    )) + (Nk * (k - 1)) + N0
    #n00 = n0 + NHY(i,     j,     k,     Nx, Ny, Nz)
    #n01 = n0 + NHY(i - 1, j,     k,     Nx, Ny, Nz)
    #n10 = n0 + NHY(i,     j,     k - 1, Nx, Ny, Nz)
    #n11 = n0 + NHY(i - 1, j - 1, k - 1, Nx, Ny, Nz)
    chy = (cHy[n00] + cHy[n01] + cHy[n10] + cHy[n11]) / 4

    # Hz
    n00 = n0 + (Ni * (i    )) + (Nj * (j    )) + (Nk * (k    )) + N0
    n01 = n0 + (Ni * (i    )) + (Nj * (j - 1)) + (Nk * (k    )) + N0
    n10 = n0 + (Ni * (i - 1)) + (Nj * (j    )) + (Nk * (k    )) + N0
    n11 = n0 + (Ni * (i - 1)) + (Nj * (j - 1)) + (Nk * (k    )) + N0
    #n00 = n0 + NHZ(i,     j,     k,     Nx, Ny, Nz)
    #n01 = n0 + NHZ(i,     j - 1, k,     Nx, Ny, Nz)
    #n10 = n0 + NHZ(i - 1, j,     k,     Nx, Ny, Nz)
    #n11 = n0 + NHZ(i - 1, j - 1, k,     Nx, Ny, Nz)
    chz = (cHz[n00] + cHz[n01] + cHz[n10] + cHz[n11]) / 4

    return chx, chy, chz
