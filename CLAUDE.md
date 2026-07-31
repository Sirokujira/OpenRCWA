# CLAUDE.md — OpenRCWA

RCWA (厳密結合波解析) + FDTD 互換入力の光学ソルバー (C++/C)。
OpenFDTD-X (GUI) から QProcess で起動される処理カーネル。

実行ファイルは 2 系統ある:

- **`orcwa`** — FDTD ソルバ。`.ofd` 入力。`rcwa` / `rcwalayer` キーがあれば
  `sol/rcwa_bridge.cpp` 経由で RCWA コアに分岐し `rcwa_efficiency.csv` を出す。
- **`orcwa_rcwa`** — RCWA 専用ドライバ。`.orcwa` 入力 (`rcwa/RCWAInput` /
  `rcwa/RCWADriver`)。斜め入射・2D 格子・任意偏波・場出力に対応。

## ビルド / テスト

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release -DWITH_CUDA=OFF -DWITH_MPI=OFF
cmake --build build -j"$(nproc)"

# dipole スモーク (FDTD)
mkdir -p /tmp/smoke && cp data/sample/dipole.ofd /tmp/smoke/ && cd /tmp/smoke
$OLDPWD/bin/orcwa -n 2 dipole.ofd && grep "normal end" orcwa.log

# RCWA コア (grating): 全周波数・両偏波で |R+T-1| <= 3e-3 (エネルギー保存)
cp data/sample/grating.ofd /tmp/smoke/ && cd /tmp/smoke
$OLDPWD/bin/orcwa -n 2 grating.ofd && cat rcwa_efficiency.csv

# RCWA 単体テスト (パーサ/ドライバ/場出力)
./bin/tests/test_rcwa_input
```

- `test_rcwa_input` の終了コード = 失敗数。`All tests PASSED` が成功の目印。
  テストは `CHECK` マクロで失敗を数え、`main` の戻り値 = 失敗数。
- ソルバの進捗ログが大量に出るため、確認時は
  `grep -v "compute interface\|begin to solve\|scatter matrices\|layer state"` で除去する。
- **物理変更をしたら必ずエネルギー保存 (無損失で R+T≈1) のテストを通すこと。**
  解析解 (Fresnel / Brewster / エネルギー保存) との比較を
  `tests/test_rcwa_input.cpp` に追加する。
- 新キーワードはパース単体テスト + solve までの結合テストの両方を書く。

## 実行 (orcwa_rcwa)

```bash
./bin/orcwa_rcwa [-o out.csv] [-v] \
    [-field Ez -slice xz -scoord 0 -sopt mod] [-device] input.orcwa
```

## ディレクトリ構成

| パス | 内容 |
|---|---|
| `rcwa/` | RCWA エンジン本体 (`RCWASolver`, `Layer`, Fourier ソルバ, `MKLEigenSolver.hpp`) と `.orcwa` パーサ (`RCWAInput`) / ドライバ (`RCWADriver`) |
| `sol/` | FDTD 処理本体 (C) + `rcwa_bridge.cpp` (`.ofd` → RCWA コア) |
| `src/` | 実行ファイルの main (`sol_Main.c` = orcwa, `rcwa_Main.cpp` = orcwa_rcwa) |
| `core/` | 積分器・レイヤサンプラ (導波モード基底の別系統。テストからのみ到達) |
| `include/` | FDTD 共通ヘッダ |
| `python/` | FDTD ソルバ (Python + Numba)・入力生成スクリプト |
| `cuda/`, `cuda_mpi/`, `gpu/` | GPU 実装 |
| `gdstk/` | GDSII 読み込みライブラリ (ベンダリング) |
| `post/` | ポスト処理 (`orcwa_post`) |
| `tests/` | 単体テスト (`test_rcwa_input.cpp` が RCWA 系の主テスト) |
| `ci/`, `data/sample/` | CI 用スモークテスト入力 |

## RCWA コアの規則 (数値安定性 — 実際に踏んだもの)

- `.ofd` 入力キー: `rcwa = <N> <period[m]>` / `rcwalayer = <eps1> <eps2> <fill>
  <thickness[m]>`。**厚みは物理長** (γ = k0·neff スケーリング前提)。
- 固有値ソルバーは `rcwa/MKLEigenSolver.hpp` の 3 分岐
  (MKL / LAPACKE / Eigen フォールバック)。**±m 次数の縮退固有値**で
  固有ベクトル基底が不良条件になり回折効率が発散する事象があるため、
  Eigen 経路には対角への非一様微小摂動 (~1e-11·‖A‖)、LAPACKE 経路には
  残差検算 ‖A·V−V·D‖ + Eigen 解き直しが入っている。**外さない**。
- `EIGEN_DONT_PARALLELIZE` は OpenMP との競合回避。外さない。
- エネルギー保存 (R+T=1) チェックを新機能でも維持する。

## 重要な物理・実装規約

詳細は `.claude/rules/rcwa-conventions.md` を参照 (rcwa/ 配下を触ると自動ロード)。
エージェント向けの同内容は `AGENTS.md` にもある。要点:

- **時間規約は exp(−iωt)**: 損失媒質は誘電率の**正**の虚部。負の虚部は利得になる。
- **内部単位は μm**: `.orcwa` はメートル指定で、パース時に μm へ変換される。
- **エネルギー正規化は `1/kzinc`**: 偏波ベクトルと対で決まっており、片方だけ
  変えると斜め入射で R+T≠1 になる。
- **材料インデックス**: 0=空気, 1=PEC, 2以降がユーザ定義 (`material` 行の順)。
  複素誘電率は `material_eps = m epsr epsi` / `material_index = m n k` で直接指定できる。
- **RCWA は周期境界のみ**: `pbc=0` は警告の上で無視される。
- **偏波指定**: `planewave = θ φ pol [psi]` — pol は 1=TM(p), 2=TE(s),
  3=psi[deg] の直線偏波 (0=TM, 90=TE), 4=右円偏波, 5=左円偏波。
- **磁性材料**: `material = type epsr esgm amur msgm` の amur/msgm が効く。
  `material_mu = m mur [mui]` でも直接指定できる。μ=1 の問題は高速経路を通る。
- **多極分散**: `material_dispersion` を同じ材料に複数行書くと極が**加算**される
  (`eps = einf + Σ_p ae_p²/(ce_p²−ω²−i·be_p·ω)`)。導電率項とも加算される。

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
  LAPACKE 経路 (残差検算つき) の方が堅牢なので、vcpkg で LAPACKE を
  入れられるなら移行が望ましい。
- macOS は Homebrew lapack (keg-only) の LAPACKE を明示指定
  (Eigen フォールバックは縮退問題があるため LAPACKE を使う)。

## CI

`.github/workflows/ci.yml`: Linux / macOS / Windows。
Linux ジョブは FDTD スモーク・RCWA grating スモークに加えて
`test_rcwa_input` (単体) と `orcwa_rcwa` の Fresnel スモークを実行する。
タグ `v*` push で Release にバイナリ添付。

## Git 運用

- 指定された作業ブランチで開発し `git push -u origin <branch>` でプッシュする。
- 明示的に依頼されない限り Pull Request は作成しない。
- コミットメッセージは変更内容を具体的に (物理的な修正は理由も書く)。
