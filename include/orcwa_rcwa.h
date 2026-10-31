/*
orcwa_rcwa.h

RCWA モードの共有定義。
.ofd の入力 (sol/input_data.c) とブリッジ (sol/rcwa_bridge.cpp) の橋渡し。
NRcwaLayer >= 2 のとき orcwa は FDTD を実行せず RCWA コアで解く。
*/
#ifndef ORCWA_RCWA_H
#define ORCWA_RCWA_H

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

/* 1層 = 1 周期の中に eps1 の矩形が 1 つある格子。
   eps1 : x = 0 .. fill*Λx かつ y = 0 .. fill_y*Λy の比誘電率
   eps2 : それ以外の比誘電率
   fill_y = 1 (既定) なら y 方向に一様な 1D 線格子 (従来の書式と同じ)。
   fill_y < 1 は 2D 格子 (rcwa に Ny Λy を書いたとき) でだけ使える。
   一様層は eps1 = eps2 とする

   誘電率の虚部 (eps1i/eps2i) は損失を表す。時間規約は exp(-iwt) なので
   損失は「正」の虚部。負にすると利得媒質になり T>1 が出るので注意。
   金属など実部が負の材料も指定できる。 */
typedef struct {
	double eps1;
	double eps2;
	double eps1i;      /* eps1 の虚部 (省略時 0)。正 = 損失 */
	double eps2i;      /* eps2 の虚部 (省略時 0)。正 = 損失 */
	double fill;       /* eps1 の x 方向の占有率 (1D: 0 < fill < 1, 2D: 0 < fill <= 1) */
	double fill_y;     /* eps1 の y 方向の占有率 (0 < fill_y <= 1, 省略時 1) */
	double thickness;  /* 層厚 [m] (先頭/末尾の半無限層では無視) */
} rcwalayer_t;

extern int          NRcwaHarmonics;  /* x 方向の空間高調波次数 Nx (計算次数は 2Nx+1) */
extern double       RcwaPeriod;      /* x 方向の格子周期 Λx [m] */
/* y 方向。rcwa = Nx Λx Ny Λy と書いたときだけ 2D 格子になる。
   省略時は NRcwaHarmonicsY = 0 (y 方向は一様) / RcwaPeriodY = RcwaPeriod */
extern int          NRcwaHarmonicsY;
extern double       RcwaPeriodY;
extern int          NRcwaLayer;      /* 層数 (0 = FDTD 経路) */
extern rcwalayer_t  *RcwaLayer;      /* 層スタック (入射側 → 透過側) */
/* 入射方向 [deg]。planewave キーから複写する (無指定なら 0 = 法線入射)。
   planewave の pol は使わない — RCWA モードは常に TE/TM 両方を出力する。 */
extern double       RcwaTheta;
extern double       RcwaPhi;

/* frequency1 掃引の回折効率を計算し rcwa_efficiency.csv に出力する。
   正常終了で 0 */
int rcwa_run(FILE *fp_log);

#ifdef __cplusplus
}
#endif

#endif /* ORCWA_RCWA_H */
