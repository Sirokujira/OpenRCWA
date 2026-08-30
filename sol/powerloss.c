/*
powerloss.c (OpenRCWA)

DFT 後の周波数領域の場から発熱量密度 [W/m^3] を求める。
直列 / MPI / CUDA / CUDA+MPI のどのビルドでも同じ結果になるよう
1 箇所に集約する。
*/

#include "orcwa.h"
#include "orcwa_prototype.h"

/*
時間平均の消費電力密度:

    P = 1/2 * sigma_e * |E|^2 + 1/2 * sigma_m * |H|^2

Yee 格子では E/H の各成分が別の位置にあり、成分ごとに材料 ID
(iEx/iEy/iEz, iHx/iHy/iHz) が異なる。したがって導電率も成分ごとに引く。

以前の実装は material_id = 0 (= 空気) の導電率を全セルに使い、さらに磁気損失を
mu'' = 1e-3 のマジック定数にしていた。そのため PEC + 空気だけの無損失問題
(data/sample/dipole.ofd) でも P_loss が最大 2.5e11 W/m^3 を返していた
(正しくは全セル厳密に 0)。

PEC (材料 ID 1) は C1 = C2 = 0 で E が 0 に固定されるため、導電率は 0 のままで
よい (完全導体は電力を消費しないモデル)。空気 (材料 ID 0) も esgm = msgm = 0。

振幅の規約: 係数 1/2 は cEx_r 等が波高値 (peak) の位相子であることを前提にする
(従来の実装と同じ)。

注意: DFT 配列 (cEx_r 等) は float。double * で受けると読み越しになり
Windows ではアクセス違反になる (glibc では偶然動作していた)。
*/
void calculatePowerLoss(double *P_loss)
{
    for (int ifreq = 0; ifreq < NFreq2; ifreq++) {
        const int64_t n0 = (int64_t)ifreq * NN;
        for (int64_t nn = 0; nn < NN; nn++) {
            const int64_t n = n0 + nn;

            const double ex2 = (double)cEx_r[n] * cEx_r[n] + (double)cEx_i[n] * cEx_i[n];
            const double ey2 = (double)cEy_r[n] * cEy_r[n] + (double)cEy_i[n] * cEy_i[n];
            const double ez2 = (double)cEz_r[n] * cEz_r[n] + (double)cEz_i[n] * cEz_i[n];

            const double hx2 = (double)cHx_r[n] * cHx_r[n] + (double)cHx_i[n] * cHx_i[n];
            const double hy2 = (double)cHy_r[n] * cHy_r[n] + (double)cHy_i[n] * cHy_i[n];
            const double hz2 = (double)cHz_r[n] * cHz_r[n] + (double)cHz_i[n] * cHz_i[n];

            P_loss[n] =
                0.5 * (Material[iEx[nn]].esgm * ex2
                     + Material[iEy[nn]].esgm * ey2
                     + Material[iEz[nn]].esgm * ez2)
              + 0.5 * (Material[iHx[nn]].msgm * hx2
                     + Material[iHy[nn]].msgm * hy2
                     + Material[iHz[nn]].msgm * hz2);
        }
    }
}
