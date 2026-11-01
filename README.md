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
  回折効率を `rcwa_efficiency.csv` (TE/TM 両偏波の R, T) と
  `rcwa_orders.csv` (回折次数ごと) に出力します。斜め入射・円錐回折、
  1D / 2D 格子 (層ごとに矩形 1 つ)、複素誘電率と Lorentz 分散に対応します。
  ブリッジは `sol/rcwa_bridge.cpp`。Fresnel・Brewster・1/4 波長 AR コート・
  金属薄膜の特性行列法 (解析解) と一致し、`.orcwa` 経路とも一致することを
  CI で常時検証しています (無損失なら R+T=1)。
- **FDTD モード** — それ以外の入力は OpenFDTD 由来の FDTD 経路で
  動作します (従来通り)。

`orcwa` とは別に、RCWA 専用のドライバ **`orcwa_rcwa`** があります。
こちらは `.orcwa` 形式 (テキスト、`.ofd` と共通の書式) を直接受け取り、
`orcwa` の RCWA モード (`.ofd`) では扱えない以下に対応します:

- 任意の直線偏波・円偏波 (`planewave = θ φ pol [psi]`。`.ofd` は TE/TM 固定)
- 円柱・球などの形状レイヤ、GDSII 由来パターン (`.ofd` は層ごとに矩形 1 つ)
- 磁性材料と導電率による損失 (`material` の amur/msgm/esgm)
- 断面場 (Ex..Hz) の CSV / 画像出力

```sh
./bin/orcwa_rcwa [-o out.csv] [-v] \
    [-field Ez -slice xz -scoord 0 -sopt mod] [-device] input.orcwa
```

パーサは `rcwa/RCWAInput`、ドライバは `rcwa/RCWADriver`。単体テストは
`bin/tests/test_rcwa_input` (Fresnel / Brewster / エネルギー保存の
解析解比較を含む) で、CI の Linux ジョブが常時実行します。

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
| `tests/` | RCWA の単体テスト `test_rcwa_input` (解析解比較つき) と旧テストバイナリ (`test_cpu`, `test_gds_cpu` ほか)。`WITH_RCWA` (既定 ON) | `test_rcwa_input` は CTest 登録済み。旧テストはビルドのみ |
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
cmake --build build -j       # → bin/orcwa, bin/orcwa_post, bin/orcwa_rcwa

# CUDA / MPI 版を作る場合 (並列 HDF5 が要る: libhdf5-openmpi-dev)
cmake -B build-mpi -DCMAKE_BUILD_TYPE=Release -DWITH_MPI=ON -DHDF5_PREFER_PARALLEL=ON
cmake --build build-mpi -j   # → bin/orcwa, bin/orcwa_mpi, bin/orcwa_post
```

- **`orcwa` (CPU) はどの構成でも必ずビルドされる。** CUDA / MPI 版
  (`orcwa_cuda` / `orcwa_mpi` / `orcwa_cuda_mpi`) は FDTD 専用で、RCWA 入力
  (`rcwalayer`) は「orcwa を使え」と案内して終了する。案内先の `orcwa` は
  同じビルドの `bin/` に並んでできる。
- CUDA 13 以降は `-DCMAKE_CUDA_ARCHITECTURES=<実機の値>` を渡す
  (既定は nvcc 13 以降 75)。Windows (MSVC) での手順は OpenFDTD-X の
  `docs/windows-cuda-mpi-build.md`。
- RCWA ドライバ `orcwa_rcwa` と単体テスト `test_rcwa_input` は `WITH_RCWA`
  (既定 ON) で作られる。不要なら `-DWITH_RCWA=OFF`。

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
# rcwa = <Nx> <周期 Λx[m]> [<Ny> <周期 Λy[m]>]   (計算次数は (2Nx+1)(2Ny+1))
#   Ny Λy を書くと 2D 格子、省略すると 1D (y 方向一様)
rcwa = 5 3.0e-7
# rcwalayer = <eps1r> <eps2r> <fill> <厚さ[m]> [<eps1i> <eps2i> [<fill_y>]]
#   1 周期内で eps1 が矩形 [0, fill*Λx] x [0, fill_y*Λy] を占め、残りが eps2。
#   fill_y は省略時 1 (1D 線格子)。一様層は eps1 = eps2。層ごとに fill が違ってよい。
#   虚部は損失 (正の値、exp(-iωt))。実部が負の金属も可。
#   入射側 → 透過側の順。先頭と末尾は半無限層 (厚さは無視、均質であること)。
rcwalayer = 1.0  1.0  0.5 0
rcwalayer = 4.0  1.0  0.5 2.0e-7
rcwalayer = 2.25 2.25 0.5 0
# 分散材料 (任意): 層 (1 始まり) の eps1 (1) / eps2 (2) に Lorentz 極を与える
#   rcwadispersion = <layer> <1|2> <einf> <ae> <be> <ce>   [rad/s]
#   eps(ω) = einf + Σ ae²/(ce² − ω² − i·be·ω)   (複数行で極が加算される)
# 入射方向 (任意): θ, φ [deg]。pol は無視され、TE/TM の両方を出力する
planewave = 0 0 1
# 波長は frequency1 から lambda = c/f で決まる
frequency1 = 4.3e14 7.5e14 16
```

出力:

| ファイル | 内容 |
|---|---|
| `rcwa_efficiency.csv` | `frequency[Hz],lambda[m],R_TE,T_TE,R_TM,T_TM` (全次数の合計) |
| `rcwa_orders.csv` | `frequency[Hz],lambda[m],pol,m,n,R,T` (伝搬する回折次数ごと) |
| `time_series_data.h5` | `/rcwa/spectrum` と `/rcwa/R_orders`・`T_orders` ほか (GUI 表示用) |

サンプルは `data/sample/` (`grating.ofd`, `rcwa_oblique.ofd`, `rcwa_metal.ofd`,
`rcwa_rect2d.ofd`, `rcwa_drude_film.ofd` ほか。期待値は `AGENTS.md`)。
円柱など任意形状の層・磁性材料・任意偏波・断面場の出力は `.orcwa` +
`orcwa_rcwa` を使います。

RCWA の回帰テストは CTest に登録してあります:

```sh
ctest --test-dir build    # test_rcwa_input + .ofd 経路 (2D / 分散 / 次数別効率)
```

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

1. RCWA モード: 入射角の掃引はカーネル 1 回につき 1 角度 (GUI が角度ごとに
   orcwa を実行する)。円柱など任意形状の層は `.ofd` 経路では書けない
   (`.orcwa` 経路で対応)
2. 旧テストバイナリ (`tests/` の `test_cpu` / `test_gds_cpu` ほか) の基準データ
   同梱 + 許容閾値 + CTest 登録

## Reference

- OpenFDTD — http://www.e-em.co.jp/OpenFDTD/
