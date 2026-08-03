/*
rcwa_hdf5.h

RCWA の計算結果 (波長掃引スペクトル) を HDF5 に書き出す共有ライター。
.ofd 経路 (sol/rcwa_bridge.cpp) と .orcwa 経路 (src/rcwa_Main.cpp) の
両方から呼ばれ、同一のレイアウトを出力する。

既存の CSV 出力 (rcwa_efficiency.csv / -o の出力) は残したまま並行して
書き出す。CSV を読んでいる既存の流れを壊さないため。

出力レイアウト (GUI 側の読み込み仕様):

  /metadata/
      solver_mode     文字列 "RCWA"   ← FDTD 出力と区別する目印
      Title           char[256]
      NFreq           int             スペクトルの点数
      npol            int             偏波の数 (1 or 2)
      pol_labels      char[npol][16]  各偏波の名前 ("TE" / "TM" など)
      rcwa_harmonics  int             空間高調波次数 N
      rcwa_period     double          格子周期 [m] (.orcwa 経路では 0)
      rcwa_nlayer     int             層数
      theta, phi      double          入射方向 [deg]
  /rcwa/
      spectrum        compound[npol][NFreq]   rcwa_spectrum_row_t と同じ並び

spectrum を [npol][NFreq] の 2 次元にしてあるので、偏波ごとに 1 本の
曲線として素直に読める:
  - .ofd 経路   : npol=2, pol_labels = {"TE", "TM"}
  - .orcwa 経路 : npol=1, pol_labels = {入力の planewave で指定した偏波}

FDTD 経路が出力する time_series_data.h5 には /metadata/solver_mode が
無いので、GUI 側はこのデータセットの有無で RCWA / FDTD を判別できる。
*/
#ifndef RCWA_HDF5_H
#define RCWA_HDF5_H

#ifdef __cplusplus
extern "C" {
#endif

/* 偏波ラベルの最大長 (終端含む) */
#define RCWA_POL_LABEL_LEN 16

/* スペクトル 1 点。吸収 A = 1 - R - T (損失媒質で正) */
typedef struct {
	double frequency;  /* [Hz] */
	double lambda;     /* [m]  */
	double R;
	double T;
	double A;
} rcwa_spectrum_row_t;

/* 計算条件 */
typedef struct {
	const char *title;
	int    nharmonics;
	double period;     /* [m]。.orcwa 経路など不定なら 0 */
	int    nlayer;
	double theta;      /* [deg] */
	double phi;        /* [deg] */
} rcwa_meta_t;

/* rows は [npol][nfreq] の行優先。pol_labels は npol 個の文字列。
   成功で 0、失敗で 0 以外。失敗しても呼び出し側は続行してよい
   (HDF5 は CSV の補助出力なので、書けないことは致命的ではない) */
int rcwa_write_hdf5(const char *filename,
                    const rcwa_meta_t *meta,
                    const rcwa_spectrum_row_t *rows,
                    int npol,
                    int nfreq,
                    const char *const *pol_labels);

#ifdef __cplusplus
}
#endif

#endif /* RCWA_HDF5_H */
