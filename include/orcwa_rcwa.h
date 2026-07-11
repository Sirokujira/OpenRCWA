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

/* 1層 = x 方向 2 セルの 1D 格子 (y 方向は一様)。
   eps1 : x = 0 .. fill*period の比誘電率
   eps2 : x = fill*period .. period の比誘電率
   一様層は eps1 = eps2 とする */
typedef struct {
	double eps1;
	double eps2;
	double fill;       /* eps1 の占有率 (0 < fill < 1) */
	double thickness;  /* 層厚 [m] (先頭/末尾の半無限層では無視) */
} rcwalayer_t;

extern int          NRcwaHarmonics;  /* 空間高調波次数 N (計算次数は 2N+1) */
extern double       RcwaPeriod;      /* 格子周期 [m] */
extern int          NRcwaLayer;      /* 層数 (0 = FDTD 経路) */
extern rcwalayer_t  *RcwaLayer;      /* 層スタック (入射側 → 透過側) */

/* frequency1 掃引の回折効率を計算し rcwa_efficiency.csv に出力する。
   正常終了で 0 */
int rcwa_run(FILE *fp_log);

#ifdef __cplusplus
}
#endif

#endif /* ORCWA_RCWA_H */
