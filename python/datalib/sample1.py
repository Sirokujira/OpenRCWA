# -*- coding: utf-8 -*-
"""
sample1.py
OpenRCWA データ作成ライブラリ (Python) サンプルプログラム (1)

使い方:
(1) 本ファイルを編集し、orcwa_datalib.py と同じフォルダにおく
(2) > python sample1.py
(3) OpenRCWA入力ファイル(sample1.orcwa)が出力される
"""

import orcwa_datalib as orcwa

orcwa.title('dipole antenna')

orcwa.xsection(-0.075, 30, +0.075)
orcwa.ysection(-0.075, 30, +0.075)
orcwa.zsection(-0.075, 10, -0.025, 11, 0.025, 10, 0.075)

#orcwa.material(2.0, 0.0, 1.0, 0.0)

orcwa.geometry(1, 1, 0.0, 0.0, 0.0, 0.0, -0.025, +0.025)

orcwa.feed('Z', 0.0, 0.0, 0.0, 1.0, 0.0, 50)
#orcwa.planewave(90, 0, 1)
#orcwa.pml(5, 2, 1e-5)

orcwa.frequency1(2e9, 3e9, 10)
orcwa.frequency2(3e9, 3e9, 0)

orcwa.solver(1000, 50, 1e-3)

orcwa.plotiter(1)
orcwa.plotzin(1)
orcwa.plotyin(1)
orcwa.plotref(1)

orcwa.plotfar1d('X', 90)
#orcwa.far1dstyle(2)
#orcwa.far1dcomponent(1, 1, 1)
#orcwa.far1ddb(1)
#orcwa.far1dnorm(1)
#orcwa.far1dscale(-30, 10, 4)

orcwa.plotfar2d(18, 36)
#orcwa.far2dcomponent(1, 0, 0, 0, 0, 0, 0)
#orcwa.far2ddb(1)
#orcwa.far2dscale(-20, 10, 6)
#orcwa.far2dobj(0.5)

orcwa.plotnear1d('E', 'Z', 0.03, 0.0)
#orcwa.near1ddb(1)
#orcwa.near1dnoinc(1)
#orcwa.near1dscale(-30, 20, 5)

orcwa.plotnear2d('E', 'X', 0.03)
#orcwa.near2ddim(1, 1)
#orcwa.near2dframe(20)
#orcwa.near2ddb(1)
#orcwa.near2dscale(-40, +20)
#orcwa.near2dcontour(1)
orcwa.near2dobj(1)
#orcwa.near2dnoinc(1)
#orcwa.near2dzoom(-0.1, +0.1, -0.1, +0.1)

orcwa.output('sample1.orcwa')
