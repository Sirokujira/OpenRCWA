# OpenRCWA

RCWA (厳密結合波解析) + FDTD 互換入力の光学ソルバー (C++)。
OpenFDTD-X (GUI) から QProcess で起動される処理カーネル。

## ビルド / テスト

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release -DWITH_CUDA=OFF -DWITH_MPI=OFF
cmake --build build -j"$(nproc)"

# dipole スモーク
mkdir -p /tmp/smoke && cp data/sample/dipole.ofd /tmp/smoke/ && cd /tmp/smoke
$OLDPWD/bin/orcwa -n 2 dipole.ofd && grep "normal end" orcwa.log

# RCWA コア (grating): 全周波数・両偏波で |R+T-1| <= 3e-3 (エネルギー保存)
cp data/sample/grating.ofd /tmp/smoke/ && cd /tmp/smoke
$OLDPWD/bin/orcwa -n 2 grating.ofd && cat rcwa_efficiency.csv
```

## RCWA コアの規則 (数値安定性 — 実際に踏んだもの)

- 入力キー: `rcwa = <N> <period[m]>` / `rcwalayer = <eps1> <eps2> <fill>
  <thickness[m]>`。**厚みは物理長** (γ = k0·neff スケーリング前提)。
- 固有値ソルバーは `rcwa/MKLEigenSolver.hpp` の 3 分岐
  (MKL / LAPACKE / Eigen フォールバック)。**±m 次数の縮退固有値**で
  固有ベクトル基底が不良条件になり回折効率が発散する事象があるため、
  Eigen 経路には対角への非一様微小摂動 (~1e-11·‖A‖)、LAPACKE 経路には
  残差検算 ‖A·V−V·D‖ + Eigen 解き直しが入っている。**外さない**。
- `EIGEN_DONT_PARALLELIZE` は OpenMP との競合回避。外さない。
- エネルギー保存 (R+T=1) チェックを新機能でも維持する。

## 移植性

- C99 VLA 禁止 (MSVC)。libm は `MATH_LIB` 変数経由。
- **`WITH_RCWA` が制御するのは `tests/` のテストバイナリだけ**で、orcwa 本体の
  RCWA コア (`rcwa/*.cpp`) は常にビルドされる。Windows CI は
  `-DWITH_RCWA=OFF` だが、**grating スモーク (エネルギー保存 |R+T-1| <= 3e-3)
  は Windows でも実行・合格している** — RCWA は Windows でも動作する。
  OFF にしているのは `gdstk/utils.cpp` が LAPACK の `dgesv_` を直接呼び、
  Windows に標準 LAPACK が無いため (gdstk はテストバイナリのみが使う)。
  Windows のテストバイナリを有効化するには vcpkg で LAPACK を導入する。
- Windows の固有値計算は Eigen フォールバック経路 (縮退対策の対角摂動あり)。
  **vcpkg での LAPACKE 導入は 2 回試行してどちらも不成立** (2026-07 実測):
  `lapack` ポートは lapacke.h を提供せず、`lapack-reference` ポートは
  VS 同梱 flang が LAPACK 3.12.1 の sgedmd.f90 をコンパイルできない
  (flang のバグ)。vcpkg/LLVM 側の修正が入るまで再挑戦しないこと。
- macOS は Homebrew lapack (keg-only) の LAPACKE を明示指定
  (Eigen フォールバックは縮退問題があるため LAPACKE を使う)。

## CI

`.github/workflows/ci.yml`: Linux / macOS / Windows。
タグ `v*` push で Release にバイナリ添付。
