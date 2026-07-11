# OpenRCWA

周期構造の解析を目的とした RCWA (厳密結合波解析) プロジェクト。
[OpenFDTD](https://github.com/Sirokujira/OpenFDTD) のコード構造
(入力データ・メッシュ・材料・後処理) をベースに、Eigen/MKL による
RCWA コア (`rcwa/`) を持ちます。GUI フロントエンド
[OpenFDTD-X](https://github.com/Sirokujira/OpenFDTD-X) から
subprocess として起動されます。

`orcwa` は入力に応じて 2 つの経路で動作します:

- **RCWA モード** — `.ofd` に `rcwa` / `rcwalayer` キーがあると、
  `rcwa/RCWASolver` (平面波回折) で層スタックを解き、周波数掃引の
  回折効率を `rcwa_efficiency.csv` に出力します (垂直入射、TE/TM 両偏波)。
  ブリッジは `sol/rcwa_bridge.cpp`。単一界面・薄膜干渉・1/4 波長 AR コートで
  フレネル解析解と一致することを検証済み、CI でエネルギー保存
  (R+T=1) を常時チェックします。
- **FDTD モード** — それ以外の入力は OpenFDTD 由来の FDTD 経路で
  動作します (従来通り)。

なお `core/` の積分器 (`RCWAIntegrator`, `RedhefferIntegrator`,
`DifferentialIntegrator`) は導波モード基底の別系統で、`tests/` の
テストバイナリからのみ利用できます。

## 処理部の構成

| ディレクトリ | 役割 | 状態 |
|---|---|---|
| `src/sol_Main.c` | ソルバー `orcwa` のエントリ。`rcwalayer` があれば RCWA モードへ分岐 | 動作 |
| `sol/rcwa_bridge.cpp` | RCWA ブリッジ: 入力キー → `RCWASolver` → `rcwa_efficiency.csv` | 動作 |
| `sol/` (その他) | FDTD 処理本体 (Yee 更新/境界条件/DFT/出力)。OpenFDTD と同系 + HDF5/熱解析レイヤ (実験的) | 動作 |
| `rcwa/` | RCWA コア: 層状周期構造の固有値問題 (`FourierSolver1D/2D`, `Layer`, `RCWASolver`)。MKL があれば LAPACKE、無ければ Eigen ソルバへフォールバック | `orcwa` に結線済み |
| `core/` | RCWA の積分器 (`RCWAIntegrator` / `RedhefferIntegrator` / `DifferentialIntegrator`)。導波モード基底の別系統 | テストからのみ到達可 |
| `gdstk/` | GDS 形状の取り扱い (gdstk 由来) | rcwa 用 |
| `tests/` | RCWA コアのテストバイナリ (`test_cpu`, `test_gds_cpu` ほか)。`-DWITH_RCWA_TESTS=ON` | ビルドのみ (CTest 未登録・アサーション整備は今後) |
| `python/` | 入力生成スクリプト (`datalib/sample_rcwa_grating.py` など) | FDTD 前提 |
| `post/` | ポスト処理 (ev2/ev3、周波数/遠方界/近傍界プロット、HDF5 読込) | 動作 |

### 依存関係 (CPU 版)

- 必須: C99/C++17 コンパイラ、CMake 3.18+、OpenMP、libhdf5、LAPACK、Eigen3
- オプション: MKL (あれば LAPACKE 経由の固有値ソルバ)、TBB、CUDA
- macOS: LAPACK は Accelerate framework、OpenMP は Homebrew libomp。
  **RCWA モードでは Homebrew `lapack` (LAPACKE) の導入を推奨** —
  Eigen の固有値ソルバへのフォールバックは ±m 次数の縮退固有値で
  基底が不良条件になることがあり、回折効率が不安定になる
  (CI の Configure ステップの `-DLAPACKE_*` 指定を参照)

## ビルド

```sh
cmake -B build -DCMAKE_BUILD_TYPE=Release -DWITH_CUDA=OFF -DWITH_MPI=OFF
cmake --build build -j       # → bin/orcwa, bin/orcwa_post

# RCWA コアのテストバイナリも作る場合
cmake -B build -DWITH_RCWA_TESTS=ON
```

## 実行

入力は OpenFDTD 形式 (`.ofd`、ヘッダ行は `OpenRCWA 4 2`)。

### RCWA モード (回折効率)

```sh
cp data/sample/grating.ofd /tmp && cd /tmp
/path/to/bin/orcwa -n 4 grating.ofd
cat rcwa_efficiency.csv    # frequency, lambda, R_TE, T_TE, R_TM, T_TM
```

入力キー:

```
# rcwa = <空間高調波次数 N> <格子周期[m]>       (計算次数は 2N+1)
rcwa = 5 3.0e-7
# rcwalayer = <eps1> <eps2> <fill> <厚さ[m]>    (入射側 → 透過側の順)
#   eps1: x = 0 .. fill*period, eps2: 残り。一様層は eps1 = eps2。
#   先頭と末尾は半無限層 (厚さは無視)。
rcwalayer = 1.0  1.0  0.5 0
rcwalayer = 4.0  1.0  0.5 2.0e-7
rcwalayer = 2.25 2.25 0.5 0
# 波長は frequency1 から lambda = c/f で決まる
frequency1 = 4.3e14 7.5e14 16
```

制約: 垂直入射のみ、y 方向一様 (1D 格子)、損失なし誘電体。
斜入射・2D パターン・GDS 形状は `tests/` の別系統が対象。

### FDTD モード

```sh
cp data/sample/dipole.ofd /tmp && cd /tmp
/path/to/bin/orcwa -n 4 dipole.ofd
grep "normal end" orcwa.log
/path/to/bin/orcwa_post -n 4 dipole.ofd   # ev.ev2/ev.ev3 等を生成
```

## CI / Release

- push / PR ごとに Linux (gcc + LAPACKE/Eigen) と macOS
  (AppleClang + Accelerate + Eigen 5.x) で CPU ビルド + スモーク実行
  (dipole = FDTD 経路、grating = RCWA 経路 + エネルギー保存チェック)
- artifact (`orcwa-linux-x64` / `orcwa-macos-arm64`) 保存、
  `v*` タグ push で Release にバイナリ自動添付

## 今後の課題

1. RCWA モードの拡張: 斜入射 (入射角掃引)、2D パターン層、
   損失性材料 (複素誘電率)、次数別効率の出力
2. `tests/` へ基準データ同梱 + 許容閾値 + CTest 登録

## Reference

- OpenFDTD — http://www.e-em.co.jp/OpenFDTD/
