# -*- coding: utf-8 -*-
"""
orcwa_datalib.py
OpenRCWA データ作成ライブラリ (Python)
"""

version = 'OpenRCWA 4 2'

# === 計算 ===

_title = None
_xsection = None
_ysection = None
_zsection = None
_m = []
_g = []
_feed = []
_planewave = None
_load = []
_point = []
_rfeed = None
_pulsewidth = None
_timestep = None
_abc = None
_pbc = None
_freq1 = None
_freq2 = None
_solver = None

def title(title):
    global _title
    _title = title

def xsection(*args):
   if len(args) % 2 == 0:
       print('*** xsectionの引数が偶数個です')
       return
   global _xsection
   _xsection = args

def ysection(*args):
   if len(args) % 2 == 0:
       print('*** ysectionの引数が偶数個です')
       return
   global _ysection
   _ysection = args

def zsection(*args):
   if len(args) % 2 == 0:
       print('*** zsectionの引数が偶数個です')
       return
   global _zsection
   _zsection = args

def material(epsr, esgm, amur, msgm):
    global _m
    _m.append([1, epsr, esgm, amur, msgm])

def material_dispersion(einf, ae, be, ce):
    global _m
    _m.append([2, einf, ae, be, ce])

def geometry(*args):
    s = int(args[1])
    ng = 8 if ((s == 31) or (s == 32) or (s == 33)
            or (s == 41) or (s == 42) or (s == 43)
            or (s == 51) or (s == 52) or (s == 53)) else 6
    if len(args) < 2 + ng:
       print('*** geometryの引数が少ないです')
       return
    global _g
    _g.append(list(args))

def geometry_array(material, shape, g):
    s = shape
    ng = 8 if ((s == 31) or (s == 32) or (s == 33)
            or (s == 41) or (s == 42) or (s == 43)
            or (s == 51) or (s == 52) or (s == 53)) else 6
    if len(g) < ng:
       print('*** geometry_arrayの引数が少ないです')
       return
    global _g
    _g.append([material, shape] + g)

def feed(cdir, x, y, z, amp, delay, z0):
    global _feed
    _feed.append([cdir, x, y, z, amp, delay, z0])

def planewave(theta, phi, pol):
    global _planewave
    _planewave = [theta, phi, pol]

def load(cdir, x, y, z, ltype, rcl):
    global _load
    _load.append([cdir, x, y, z, ltype, rcl])

def point(cdir, x, y, z, prop=''):
    global _point
    _point.append([cdir, x, y, z, prop])

def rfeed(rfeed):
    global _rfeed
    _rfeed = rfeed

def pulsewidth(pulsewidth):
    global _pulsewidth
    _pulsewidth = pulsewidth

def timestep(timestep):
    global _timestep
    _timestep = timestep

def pml(l, m, r0):
    global _abc
    _abc = [1,  l, m, r0]

def pbc(pbcx, pbcy, pbcz):
    global _pbc
    _pbc = [pbcx, pbcy, pbcz]


def frequency1(fstart, fend, div):
    global _freq1
    _freq1 = [fstart, fend, div]

def frequency2(fstart, fend, div):
    global _freq2
    _freq2 = [fstart, fend, div]

def solver(maxiter, nout, converg):
    global _solver
    _solver = [maxiter, nout, converg]

# === ポスト処理 ===

_matchingloss = None
_plotiter = None
_plotfeed = None
_plotpoint = None
_plotsmith = None
_plotzin = None
_plotyin = None
_plotref = None
_plotspara = None
_plotcoupling = None
_plotfar0d = None
_freqdiv = None
_plotfar1d = []
_far1dstyle = None
_far1dcomponent = None
_far1ddb = None
_far1dnorm = None
_far1dscale = None
_plotfar2d = None
_far2dcomponent = None
_far2ddb = None
_far2dscale = None
_far2dobj = None
_plotnear1d = []
_near1ddb = None
_near1dnoinc = None
_near1dscale = None
_plotnear2d = []
_near2ddim = None
_near2dframe = None
_near2ddb = None
_near2dscale = None
_near2dcontour = None
_near2dobj = None
_near2dnoinc = None
_near2dzoom = None

# 周波数特性

def matchingloss(i0):
    global _matchingloss
    _matchingloss = i0

def plotiter(i0):
    global _plotiter
    _plotiter = i0

def plotfeed(i0):
    global _plotfeed
    _plotfeed = i0

def plotpoint(i0):
    global _plotpoint
    _plotpoint = i0

def plotsmith(i0):
    global _plotsmith
    _plotsmith = i0

def plotzin(scale, smin = 0, smax = 0, sdiv = 0):
    global _plotzin
    _plotzin = [scale, smin, smax, sdiv]

def plotyin(scale, smin = 0, smax = 0, sdiv = 0):
    global _plotyin
    _plotyin = [scale, smin, smax, sdiv]

def plotref(scale, smin = 0, smax = 0, sdiv = 0):
    global _plotref
    _plotref = [scale, smin, smax, sdiv]

def plotspara(scale, smin = 0, smax = 0, sdiv = 0):
    global _plotspara
    _plotspara = [scale, smin, smax, sdiv]

def plotcoupling(scale, smin = 0, smax = 0, sdiv = 0):
    global _plotcoupling
    _plotcoupling = [scale, smin, smax, sdiv]

def plotfar0d(theta, phi, scale, smin = 0, smax = 0, sdiv = 0):
    global _plotfar0d
    _plotfar0d = [theta, phi, scale, smin, smax, sdiv]

def freqdiv(freqdiv):
    global _freqdiv
    _freqdiv = freqdiv

# 遠方界1D

def plotfar1d(cdir, div, angle = 0):
    global _plotfar1d
    _plotfar1d.append([cdir, div, angle])

def far1dstyle(i0):
    global _far1dstyle
    _far1dstyle = i0

def far1dcomponent(i0, i1, i2):
    global _far1dcomponent
    _far1dcomponent = [i0, i1, i2]

def far1ddb(i0):
    global _far1ddb
    _far1ddb = i0

def far1dnorm(i0):
    global _far1dnorm
    _far1dnorm = i0

def far1dscale(smin, smax, sdiv):
    global _far1dscale
    _far1dscale = [smin, smax, sdiv]

# 遠方界2D

def plotfar2d(divtheta, divphi):
    global _plotfar2d
    _plotfar2d = [divtheta, divphi]

def far2dcomponent(i0, i1, i2, i3, i4, i5, i6):
    global _far2dcomponent
    _far2dcomponent = [i0, i1, i2, i3, i4, i5, i6]

def far2ddb(i0):
    global _far2ddb
    _far2ddb = i0

def far2dscale(smin, smax, sdiv):
    global _far2dscale
    _far2dscale = [smin, smax, sdiv]

def far2dobj(obj):
    global _far2dobj
    _far2dobj = obj

# 近傍界1D

def plotnear1d(component, cdir, p1, p2):
    global _plotnear1d
    _plotnear1d.append([component, cdir, p1, p2])

def near1ddb(i0):
    global _near1ddb
    _near1ddb = i0

def near1dnoinc(i0):
    global _near1dnoinc
    _near1dnoinc = i0

def near1dscale(smin, smax, sdiv):
    global _near1dscale
    _near1dscale = [smin, smax, sdiv]

# 近傍界2D

def plotnear2d(component, cdir, p):
    global _plotnear2d
    _plotnear2d.append([component, cdir, p])

def near2ddim(i0, i1):
    global _near2ddim
    _near2ddim = [i0, i1]

def near2dframe(i0):
    global _near2dframe
    _near2dframe = i0

def near2ddb(i0):
    global _near2ddb
    _near2ddb = i0

def near2dscale(smin, smax):
    global _near2dscale
    _near2dscale = [smin, smax]

def near2dcontour(i0):
    global _near2dcontour
    _near2dcontour = i0

def near2dobj(i0):
    global _near2dobj
    _near2dobj = i0

def near2dnoinc(i0):
    global _near2dnoinc
    _near2dnoinc = i0

def near2dzoom(p0, p1, p2, p3):
    global _near2dzoom
    _near2dzoom = [p0, p1, p2, p3]

# ファイル出力
def output(filename):
    fp = open(filename, 'wt', encoding='utf-8')

    # === 計算 ===

    # ヘッダー
    fp.write('%s\n' % version)

    # タイトル
    if _title != None:
        fp.write('title = %s\n' % _title)

    # メッシュ
    if _xsection != None:
        fp.write('xmesh =')
        for x in _xsection:
            fp.write(' %s' % x)
        fp.write('\n')
    if _ysection != None:
        fp.write('ymesh =')
        for y in _ysection:
            fp.write(' %s' % y)
        fp.write('\n')
    if _zsection != None:
        fp.write('zmesh =')
        for z in _zsection:
            fp.write(' %s' % z)
        fp.write('\n')
    
    # 物性値
    for m in _m:
        fp.write('material = %d %g %g %g %g\n' % (m[0], m[1], m[2], m[3], m[4]))

    # 物体形状
    for g in _g:
        fp.write('geometry = %d %d' % (int(g[0]), int(g[1])))
        for s in g[2:]:
            fp.write(' %g' % s)
        fp.write('\n')

    # 給電点
    for feed in _feed:
        fp.write('feed = %s %g %g %g %g %g %g\n' %
            (feed[0], feed[1], feed[2], feed[3], feed[4], feed[5], feed[6]))

    # 平面波入射
    if _planewave != None:
        fp.write('planewave = %g %g %d\n' % (_planewave[0], _planewave[1], int(_planewave[2])))

    # 負荷
    for load in _load:
        fp.write('load = %s %g %g %g %s %g\n' %
            (load[0], load[1], load[2], load[3], load[4], load[5]))

    # 観測点
    for point in _point:
        fp.write('point = %s %g %g %g %s\n' %
            (point[0], point[1], point[2], point[3], point[4]))

    # その他
    if _rfeed != None:
        fp.write('rfeed = %g\n' % _rfeed)
    if _pulsewidth != None:
        fp.write('pulsewidth = %g\n' % _pulsewidth)
    if _timestep != None:
        fp.write('timestep = %g\n' % _timestep)
    
    # PML/PBC
    if _abc != None:
        fp.write('abc = 1 %d %g %g\n' % (_abc[1], _abc[2], _abc[3]))
    if _pbc != None:
        fp.write('pbc = %d %d %d\n' % (_pbc[0], _pbc[1], _pbc[2]))
    
    # 周波数
    if _freq1 != None:
        fp.write('frequency1 = %g %g %d\n' % (_freq1[0], _freq1[1], int(_freq1[2])))
    if _freq2 != None:
        fp.write('frequency2 = %g %g %d\n' % (_freq2[0], _freq2[1], int(_freq2[2])))

    # 計算条件
    if _solver != None:
        fp.write('solver = %d %d %g\n' % (int(_solver[0]), int(_solver[1]), _solver[2]))

    # === ポスト処理 ===

    # 周波数特性

    if _matchingloss != None:
        fp.write('matchingloss = %d\n' % _matchingloss)

    if _plotiter != None:
        fp.write('plotiter = %d\n' %  _plotiter)

    if _plotfeed != None:
        fp.write('plotfeed = %d\n' % _plotfeed)

    if _plotpoint != None:
        fp.write('plotpoint = %d\n' % _plotpoint)

    if _plotsmith != None:
        fp.write('plotsmith = %d\n' % _plotsmith)

    if _plotzin != None:
        fp.write('plotzin = %d' % _plotzin[0])
        if _plotzin[0] == 2:
            fp.write(' %g %g %d' % (_plotzin[1], _plotzin[2], int(_plotzin[3])))
        fp.write('\n')

    if _plotyin != None:
        fp.write('plotyin = %d' % _plotyin[0])
        if _plotyin[0] == 2:
            fp.write(' %g %g %d' % (_plotyin[1], _plotyin[2], int(_plotyin[3])))
        fp.write('\n')

    if _plotref != None:
        fp.write('plotref = %d' % _plotref[0])
        if _plotref[0] == 2:
            fp.write(' %g %g %d' % (_plotref[1], _plotref[2], int(_plotref[3])))
        fp.write('\n')

    if _plotspara != None:
        fp.write('plotspara = %d' % _plotspara[0])
        if _plotspara[0] == 2:
            fp.write(' %g %g %d' % (_plotspara[1], _plotspara[2], int(_plotspara[3])))
        fp.write('\n')

    if _plotcoupling != None:
        fp.write('plotcoupling = %d' % _plotcoupling[0])
        if _plotcoupling[0] == 2:
            fp.write(' %g %g %d' % (_plotcoupling[1], _plotcoupling[2], int(_plotcoupling[3])))
        fp.write('\n')

    if _plotfar0d != None:
        fp.write('plotfar0d = %d %g %g' % (_plotfar0d[0], _plotfar0d[1], _plotfar0d[2]))
        if _plotfar0d[2] == 2:
            fp.write(' %g %g %d' % (_plotfar0d[3], _plotfar0d[4], int(_plotfar0d[5])))
        fp.write('\n')

    if _freqdiv != None:
        fp.write('freqdiv = %d\n' % int(_freqdiv))

    # 遠方界1D
    
    for p in _plotfar1d:
        fp.write('plotfar1d = %s %d' % (p[0], int(p[1])))
        if p[2] != 0:
            fp.write(' %g' % p[2])
        fp.write('\n')

    if _far1dstyle != None:
        fp.write('far1dstyle = %d\n' % _far1dstyle)

    if _far1dcomponent != None:
        fp.write('far1dcomponent = %d %d %d\n' % (int(_far1dcomponent[0]), int(_far1dcomponent[1]), int(_far1dcomponent[2])))

    if _far1ddb != None:
        fp.write('far1ddb = %d\n' % _far1ddb)

    if _far1dnorm!= None:
        fp.write('far1dnorm = %d\n' % _far1dnorm)

    if _far1dscale != None:
        fp.write('far1dscale = %g %g %d\n' % (_far1dscale[0], _far1dscale[1], int(_far1dscale[2])))

    # 遠方界2D

    if _plotfar2d != None:
        fp.write('plotfar2d = %d %d\n' % (int(_plotfar2d[0]), int(_plotfar2d[1])))

    if _far2dcomponent != None:
        fp.write('far2dcomponent = %d %d %d %d %d %d %d\n' %
            (int(_far2dcomponent[0]), int(_far2dcomponent[1]), int(_far2dcomponent[2]), int(_far2dcomponent[3]), int(_far2dcomponent[4]), int(_far2dcomponent[5]), int(_far2dcomponent[6])))

    if _far2ddb != None:
        fp.write('far2ddb = %d\n' % _far2ddb)

    if _far2dscale != None:
        fp.write('far2dscale = %g %g %d\n' % (_far2dscale[0], _far2dscale[1], int(_far2dscale[2])))

    if _far2dobj != None:
        fp.write('far2dobj = %g\n' % _far2dobj)

    # 近傍界1D

    for p in _plotnear1d:
        fp.write('plotnear1d = %s %s %g %g\n' % (p[0], p[1], p[2], p[3]))

    if _near1ddb != None:
        fp.write('near1ddb = %d\n' % _near1ddb)

    if _near1dnoinc != None:
        fp.write('near1dnoinc = %d\n' % _near1dnoinc)

    if _near1dscale != None:
        fp.write('near1dscale = %g %g %d\n' % (_near1dscale[0], _near1dscale[1], int(_near1dscale[2])))

    # 近傍界2D

    for p in _plotnear2d:
        fp.write('plotnear2d = %s %s %g\n' % (p[0], p[1], p[2]))

    if _near2ddim != None:
        fp.write('near2ddim = %d %d\n' % (_near2ddim[0], _near2ddim[1]))

    if _near2dframe != None:
        fp.write('near2dframe = %d\n' % _near2dframe)

    if _near2ddb != None:
        fp.write('near2ddb = %d\n' % _near2ddb)

    if _near2dscale != None:
        fp.write('near2dscale = %g %g\n' % (_near2dscale[0], _near2dscale[1]))

    if _near2dcontour != None:
        fp.write('near2dcontour = %d\n' % _near2dcontour)

    if _near2dobj != None:
        fp.write('near2dobj = %d\n' % _near2dobj)

    if _near2dnoinc != None:
        fp.write('near2dnoinc = %d\n' % _near2dnoinc)

    if _near2dzoom != None:
        fp.write('near2dzoom = %g %g %g %g\n' % (_near2dzoom[0], _near2dzoom[1], _near2dzoom[2], _near2dzoom[3]))

    fp.write('end\n')
    
    fp.close()

    print('output -> %s' % filename)
