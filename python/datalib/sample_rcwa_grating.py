# -*- coding: utf-8 -*-
"""
sample_rcwa_grating.py
OpenRCWA データ作成ライブラリ (Python) サンプルプログラム
可視光領域の1次元誘電体回折格子 (RCWA向け周期構造) の入力データ作成例

使い方:
(1) 本ファイルを編集し、orcwa_datalib.py と同じフォルダにおく
(2) > python sample_rcwa_grating.py
(3) OpenRCWA入力ファイル(sample_rcwa_grating.orcwa)が出力される

モデル概要:
  - 周期 Lambda = 0.5 um のガラス(epsr=2.25)製1次元矩形格子
  - デューティ比 50% (リッジ幅 0.25 um, 格子層厚さ 0.2 um)
  - x, y方向を周期境界(PBC)、z方向をPML吸収境界とする
  - 上方(空気側)から平面波を垂直入射 (theta=0, phi=0)
  - 可視光波長 400nm - 700nm に対応する周波数範囲で計算
"""

import orcwa_datalib as orcwa

# 単位: メートル [m]
um = 1.0e-6

# 周期 (格子定数)
Lambda = 0.5 * um

# 格子形状パラメータ
ridge_width = 0.25 * um   # リッジ幅 (デューティ比 50%)
grating_h = 0.2 * um      # 格子層の厚さ
substrate_h = 0.3 * um    # 基板の厚さ (計算領域内)
air_h = 0.5 * um          # 上部空気層の厚さ (計算領域内)

orcwa.title('1D dielectric grating (visible light, RCWA)')

# 計算領域: x, y は格子周期、z は基板 - 格子層 - 空気層
orcwa.xsection(-Lambda / 2, 20, +Lambda / 2)
orcwa.ysection(-Lambda / 2, 20, +Lambda / 2)
orcwa.zsection(-substrate_h, 15, 0.0, 10, grating_h, 15, grating_h + air_h)

# 材料: index 2 = ガラス (n = 1.5, epsr = 2.25)
orcwa.material(2.25, 0.0, 1.0, 0.0)

# 基板 (z < 0 の全面): ガラス
orcwa.geometry(2, 1, -Lambda / 2, +Lambda / 2, -Lambda / 2, +Lambda / 2, -substrate_h, 0.0)

# 格子リッジ (0 <= z <= grating_h, x方向に半周期幅): ガラス
orcwa.geometry(2, 1, -ridge_width / 2, +ridge_width / 2, -Lambda / 2, +Lambda / 2, 0.0, grating_h)

# 周期境界条件: x, y方向を周期境界、z方向はPML
orcwa.pbc(1, 1, 0)
orcwa.pml(8, 2, 1e-5)

# 平面波入射: 垂直入射 (theta=0, phi=0), pol=0
orcwa.planewave(0, 0, 0)

# 観測点: 上部(反射側, -Z伝搬)と下部(透過側, +Z伝搬)
orcwa.point('Z', 0.0, 0.0, grating_h + air_h * 0.5, '-Z')
orcwa.point('Z', 0.0, 0.0, -substrate_h * 0.5, '+Z')

# 周波数範囲: 可視光 400nm - 700nm
c0 = 2.99792458e8
f_700nm = c0 / (700 * 1e-9)
f_400nm = c0 / (400 * 1e-9)
orcwa.frequency1(f_700nm, f_400nm, 30)
orcwa.frequency2(f_700nm, f_400nm, 30)

orcwa.solver(3000, 50, 1e-3)

# ポスト処理
orcwa.plotiter(1)
orcwa.plotref(1)
orcwa.plotspara(1)

orcwa.output('sample_rcwa_grating.orcwa')
