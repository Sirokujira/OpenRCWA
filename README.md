# OpenRCWA

周期構造の解析を目的とした RCWA (厳密結合波解析) プロジェクト。
[OpenFDTD](https://github.com/Sirokujira/OpenFDTD) のコード構造
(入力データ・メッシュ・材料・後処理) をベースに、Eigen/MKL による
RCWA コア (`rcwa/`) を持ちます。GUI フロントエンド
[OpenFDTD-X](https://github.com/Sirokujira/OpenFDTD-X) から
subprocess として起動されます。

> ⚠ **現状の重要な注意**: CPU 版バイナリ `orcwa` / `orcwa_post` は
> OpenFDTD 由来の **FDTD 経路で動作しており、RCWA コアはまだ結線されて
> いません**。RCWA の計算実体 (`core/RCWAIntegrator`,
> `RedhefferIntegrator`, `DifferentialIntegrator`) は `tests/` の
> テストバイナリからのみ利用できます。`orcwa` への結線 (solve 経路の
> RCWA モード分岐と、回折次数・格子周期などの入力キー追加) は今後の
> 課題です。

## 処理部の構成

| ディレクトリ | 役割 | 状態 |
|---|---|---|
| `src/sol_Main.c` | ソルバー `orcwa` のエントリ | FDTD 経路で動作 |
| `sol/` | FDTD 処理本体 (Yee 更新/境界条件/DFT/出力)。OpenFDTD と同系 + HDF5/熱解析レイヤ (実験的) | 動作 |
| `rcwa/` | RCWA コア: 層状周期構造の固有値問題 (`FourierSolver1D/2D`, `Layer`, `RCWASolver`)。MKL があれば LAPACKE、無ければ Eigen ソルバへフォールバック | テストからのみ到達可 |
| `core/` | RCWA の積分器 (`RCWAIntegrator` / `RedhefferIntegrator` / `DifferentialIntegrator`) | テストからのみ到達可 |
| `gdstk/` | GDS 形状の取り扱い (gdstk 由来) | rcwa 用 |
| `tests/` | RCWA コアのテストバイナリ (`test_cpu`, `test_gds_cpu` ほか)。`-DWITH_RCWA_TESTS=ON` | ビルドのみ (CTest 未登録・アサーション整備は今後) |
| `python/` | 入力生成スクリプト (`datalib/sample_rcwa_grating.py` など) | FDTD 前提 |
| `post/` | ポスト処理 (ev2/ev3、周波数/遠方界/近傍界プロット、HDF5 読込) | 動作 |

### 依存関係 (CPU 版)

- 必須: C99/C++17 コンパイラ、CMake 3.18+、OpenMP、libhdf5、LAPACK、Eigen3
- オプション: MKL (あれば LAPACKE 経由の固有値ソルバ)、TBB、CUDA
- macOS: LAPACK は Accelerate framework、OpenMP は Homebrew libomp
  (LAPACKE が無いため Eigen ソルバへ自動フォールバック)

## ビルド

```sh
cmake -B build -DCMAKE_BUILD_TYPE=Release -DWITH_CUDA=OFF -DWITH_MPI=OFF
cmake --build build -j       # → bin/orcwa, bin/orcwa_post

# RCWA コアのテストバイナリも作る場合
cmake -B build -DWITH_RCWA_TESTS=ON
```

## 実行

入力は OpenFDTD 形式 (`.ofd`、ヘッダ行は `OpenRCWA 4 2`):

```sh
cp data/sample/dipole.ofd /tmp && cd /tmp
/path/to/bin/orcwa -n 4 dipole.ofd
grep "normal end" orcwa.log
/path/to/bin/orcwa_post -n 4 dipole.ofd   # ev.ev2/ev.ev3 等を生成
```

## CI / Release

- push / PR ごとに Linux (gcc + LAPACKE/Eigen) と macOS
  (AppleClang + Accelerate + Eigen 5.x) で CPU ビルド + dipole スモーク
- artifact (`orcwa-linux-x64` / `orcwa-macos-arm64`) 保存、
  `v*` タグ push で Release にバイナリ自動添付

## 今後の課題 (RCWA 結線)

1. `solve()` 経路に RCWA モード分岐を追加し、`core/RCWAIntegrator` を
   C++ ブリッジ経由で呼び出す (CMake では `RCWA_SRC_FILES` を
   `orcwa` にリンクする)
2. `.ofd` 入力に RCWA 固有キー (回折次数 / 空間高調波数 / 格子周期 /
   入射角掃引) を追加する — GUI 側 (OpenFDTD-X の OpticalTab) には
   既に対応する設定 UI がある
3. `tests/` へ基準データ同梱 + 許容閾値 + CTest 登録

## Reference

- OpenFDTD — http://www.e-em.co.jp/OpenFDTD/
