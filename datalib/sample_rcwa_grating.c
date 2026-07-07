/*
sample_rcwa_grating.c
OpenRCWA データ作成ライブラリ (C) サンプルプログラム
可視光領域の1次元誘電体回折格子 (RCWA向け周期構造) の入力データ作成例
(python/datalib/sample_rcwa_grating.py のC版)

使い方:
(1) > orcwa_sample_grating
(2) OpenRCWA入力ファイル(sample_rcwa_grating.orcwa)が出力される

モデル概要:
  - 周期 Lambda = 0.5 um のガラス(epsr=2.25)製1次元矩形格子
  - デューティ比 50% (リッジ幅 0.25 um, 格子層厚さ 0.2 um)
  - x, y方向を周期境界(PBC)、z方向を吸収境界(Mur-1次)とする
    (本ソルバではPBC使用時の吸収境界はMur-1次となる)
  - 上方(空気側)から平面波を垂直入射 (theta=0, phi=0, V偏波)
  - 可視光波長 400nm - 700nm に対応する周波数範囲で計算
*/

#include "orcwa_datalib.h"

int main(void)
{
	/* 単位: メートル [m] */
	const double um = 1.0e-6;

	/* 周期 (格子定数) */
	const double Lambda = 0.5 * um;

	/* 格子形状パラメータ */
	const double ridge_width = 0.25 * um;  /* リッジ幅 (デューティ比 50%) */
	const double grating_h   = 0.2  * um;  /* 格子層の厚さ */
	const double substrate_h = 0.3  * um;  /* 基板の厚さ (計算領域内) */
	const double air_h       = 0.5  * um;  /* 上部空気層の厚さ (計算領域内) */

	orcwa_init();

	orcwa_title("1D dielectric grating (visible light, RCWA)");

	/* 計算領域: x, y は格子周期、z は基板 - 格子層 - 空気層 */
	orcwa_xsection(2, -Lambda / 2, +Lambda / 2);
	orcwa_xdivision(1, 20);
	orcwa_ysection(2, -Lambda / 2, +Lambda / 2);
	orcwa_ydivision(1, 20);
	orcwa_zsection(4, -substrate_h, 0.0, grating_h, grating_h + air_h);
	orcwa_zdivision(3, 15, 10, 15);

	/* 材料: index 2 = ガラス (n = 1.5, epsr = 2.25) */
	orcwa_material(2.25, 0.0, 1.0, 0.0);

	/* 基板 (z < 0 の全面): ガラス */
	orcwa_geometry(2, 1, -Lambda / 2, +Lambda / 2, -Lambda / 2, +Lambda / 2, -substrate_h, 0.0);

	/* 格子リッジ (0 <= z <= grating_h, x方向に半周期幅): ガラス */
	orcwa_geometry(2, 1, -ridge_width / 2, +ridge_width / 2, -Lambda / 2, +Lambda / 2, 0.0, grating_h);

	/* 周期境界条件: x, y方向を周期境界 */
	/* (z方向は吸収境界。PBC使用時はMur-1次が用いられるためPMLは指定しない) */
	orcwa_pbc(1, 1, 0);

	/* 平面波入射: 垂直入射 (theta=0, phi=0), pol=1 (V偏波) */
	orcwa_planewave(0, 0, 1);

	/* 観測点: 上部(反射側, -Z伝搬)と下部(透過側, +Z伝搬) */
	orcwa_point('Z', 0.0, 0.0, grating_h + air_h * 0.5, "-Z");
	orcwa_point('Z', 0.0, 0.0, -substrate_h * 0.5, "+Z");

	/* 周波数範囲: 可視光 400nm - 700nm */
	{
		const double c0 = 2.99792458e8;
		const double f_700nm = c0 / (700 * 1e-9);
		const double f_400nm = c0 / (400 * 1e-9);
		orcwa_frequency1(f_700nm, f_400nm, 30);
		orcwa_frequency2(f_700nm, f_400nm, 30);
	}

	orcwa_solver(3000, 50, 1e-3);

	/* ポスト処理 */
	orcwa_plotiter(1);
	orcwa_plotref(1, 0, 0, 0);
	orcwa_plotspara(1, 0, 0, 0);

	orcwa_outdata("sample_rcwa_grating.orcwa");

	return 0;
}
