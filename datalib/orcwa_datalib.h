/*
orcwa_datalib.h
OpenRCWA データ作成ライブラリ (C) (OpenFDTD の ofd_datalib に準拠)
*/

#ifndef _ORCWA_DATALIB_H_
#define _ORCWA_DATALIB_H_

#ifdef __cplusplus
extern "C" {
#endif

/* === 計算 === */

void orcwa_init(void);
void orcwa_title(const char title[]);

/* メッシュ: section は区間端の座標(double)、division は各区間の分割数(int) */
void orcwa_xsection(int n, ...);
void orcwa_ysection(int n, ...);
void orcwa_zsection(int n, ...);
void orcwa_xdivision(int n, ...);
void orcwa_ydivision(int n, ...);
void orcwa_zdivision(int n, ...);

/* 物性値 */
void orcwa_material(double epsr, double esgm, double amur, double msgm);
void orcwa_material_dispersion(double einf, double ae, double be, double ce);

/* 物体形状: shape に応じて 6 個または 8 個の double を続ける */
void orcwa_geometry(int material, int shape, ...);

/* 給電点/平面波入射/負荷/観測点 */
void orcwa_feed(char dir, double x, double y, double z, double amp, double delay, double z0);
void orcwa_planewave(double theta, double phi, int pol);
void orcwa_load(char dir, double x, double y, double z, char type, double rcl);
void orcwa_point(char dir, double x, double y, double z, const char prop[]);

/* その他 */
void orcwa_rfeed(double rfeed);
void orcwa_pulsewidth(double pulsewidth);
void orcwa_timestep(double timestep);

/* 吸収境界/周期境界 */
void orcwa_pml(int l, double m, double r0);
void orcwa_pbc(int pbcx, int pbcy, int pbcz);

/* 周波数 */
void orcwa_frequency1(double fstart, double fend, int div);
void orcwa_frequency2(double fstart, double fend, int div);

/* 計算条件 */
void orcwa_solver(int maxiter, int nout, double converg);

/* === ポスト処理 === */

/* 周波数特性 */
void orcwa_matchingloss(int i0);
void orcwa_plotiter(int i0);
void orcwa_plotfeed(int i0);
void orcwa_plotpoint(int i0);
void orcwa_plotsmith(int i0);
/* scale == 2 のときのみ smin, smax, sdiv が出力される */
void orcwa_plotzin(int scale, double smin, double smax, int sdiv);
void orcwa_plotyin(int scale, double smin, double smax, int sdiv);
void orcwa_plotref(int scale, double smin, double smax, int sdiv);
void orcwa_plotspara(int scale, double smin, double smax, int sdiv);
void orcwa_plotcoupling(int scale, double smin, double smax, int sdiv);
void orcwa_plotfar0d(double theta, double phi, int scale, double smin, double smax, int sdiv);
void orcwa_freqdiv(int freqdiv);

/* 遠方界1D: angle == 0 のときは angle は出力されない */
void orcwa_plotfar1d(char dir, int div, double angle);
void orcwa_far1dstyle(int i0);
void orcwa_far1dcomponent(int i0, int i1, int i2);
void orcwa_far1ddb(int i0);
void orcwa_far1dnorm(int i0);
void orcwa_far1dscale(double smin, double smax, int sdiv);

/* 遠方界2D */
void orcwa_plotfar2d(int divtheta, int divphi);
void orcwa_far2dcomponent(int i0, int i1, int i2, int i3, int i4, int i5, int i6);
void orcwa_far2ddb(int i0);
void orcwa_far2dscale(double smin, double smax, int sdiv);
void orcwa_far2dobj(double obj);

/* 近傍界1D */
void orcwa_plotnear1d(const char component[], char dir, double p1, double p2);
void orcwa_near1ddb(int i0);
void orcwa_near1dnoinc(int i0);
void orcwa_near1dscale(double smin, double smax, int sdiv);

/* 近傍界2D */
void orcwa_plotnear2d(const char component[], char dir, double p);
void orcwa_near2ddim(int i0, int i1);
void orcwa_near2dframe(int i0);
void orcwa_near2ddb(int i0);
void orcwa_near2dscale(double smin, double smax);
void orcwa_near2dcontour(int i0);
void orcwa_near2dobj(int i0);
void orcwa_near2dnoinc(int i0);
void orcwa_near2dzoom(double p0, double p1, double p2, double p3);

/* ファイル出力 */
void orcwa_outdata(const char filename[]);

#ifdef __cplusplus
}
#endif

#endif /* _ORCWA_DATALIB_H_ */
