# AGENTS.md — OpenRCWA

OpenRCWA は周期構造の電磁界解析ソルバ集。RCWA (厳密結合波解析) エンジンと
OpenFDTD 由来の FDTD ソルバ (C / Python / CUDA) を含む。

実行ファイルは 2 系統ある:

- **`orcwa`** — FDTD ソルバ。`.ofd` 入力。`rcwa` / `rcwalayer` キーがあれば
  `sol/rcwa_bridge.cpp` 経由で RCWA コアに分岐し `rcwa_efficiency.csv` を出す。
- **`orcwa_rcwa`** — RCWA 専用ドライバ。`.orcwa` 入力 (`rcwa/RCWAInput` /
  `rcwa/RCWADriver`)。斜め入射・2D 格子・任意偏波・場出力に対応。

このファイルは AI エージェント向けの作業規約。**「物理規約」「落とし穴」の節は、
知らずに触ると結果が静かに壊れる箇所**なので、`rcwa/` 配下を変更する前に必ず読むこと。

## ビルド

```bash
# 依存パッケージ (Ubuntu)
sudo apt-get install -y cmake gcc g++ libhdf5-dev \
                        liblapack-dev liblapacke-dev libeigen3-dev

# 構成 (MKL/TBB/CUDA はなくてもよい — 自動で Eigen ソルバにフォールバック)
cmake -B build -DCMAKE_BUILD_TYPE=Release -DWITH_CUDA=OFF -DWITH_MPI=OFF
cmake --build build -j"$(nproc)"
```

生成物は `bin/` (実行ファイル) と `bin/tests/` (テスト) に置かれる。
`rcwa/*.cpp` の RCWA コアは `orcwa` に組み込まれるため常にビルドされる。
ビルドオプションは 2 段階:

| オプション | 既定 | 制御対象 |
|---|---|---|
| `WITH_RCWA` | ON | `orcwa_rcwa` ドライバ + `test_rcwa_input` |
| `WITH_RCWA_LEGACY_TESTS` | ON (MSVC は OFF) | `test_gds_cpu` / `test_cpu` / `test_metasurface_cpu` / `test_cpu_diff` |

`WITH_RCWA_LEGACY_TESTS` が MSVC で既定 OFF なのは、これらが依存する
`gdstk/utils.cpp` が LAPACK の `dgesv_` を直接呼び、Windows に標準 LAPACK が
無いため。**`orcwa_rcwa` と `test_rcwa_input` は `rcwa/*.cpp` にしか依存しない
ので Windows でもビルド・実行できる** (CI で実行している)。

## テスト

```bash
./bin/tests/test_rcwa_input   # RCWA 単体テスト (パーサ/ドライバ/場出力)

# .ofd 経路の RCWA スモーク (全周波数・両偏波で |R+T-1| <= 3e-3)
mkdir -p /tmp/smoke && cp data/sample/grating.ofd /tmp/smoke/ && cd /tmp/smoke
$OLDPWD/bin/orcwa -n 2 grating.ofd && cat rcwa_efficiency.csv
```

- 終了コード = 失敗数。`All tests PASSED` が成功の目印。
- テストは `CHECK` マクロで失敗を数え、`main` の戻り値 = 失敗数。
- ソルバの進捗ログが大量に出るため、確認時は
  `grep -v "compute interface\|begin to solve\|scatter matrices\|layer state"` で除去する。
- **物理変更をしたら必ずエネルギー保存 (無損失で R+T≈1) のテストを通すこと。**
  解析解 (Fresnel / Brewster / エネルギー保存) との比較テストを
  `tests/test_rcwa_input.cpp` に追加する。
- 新キーワードはパース単体テスト + solve までの結合テストの両方を書く。
- CI は `.github/workflows/ci.yml` (ビルド + スモーク + 単体テスト)。
  `test_rcwa_input` と `orcwa_rcwa` の Fresnel スモークは **Linux / macOS /
  Windows の 3 ジョブすべて**で実行する (固有値の経路が Linux/macOS =
  LAPACKE、Windows = Eigen フォールバックと異なるため)。

## 実行

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

## 物理規約 (違反すると結果が静かに壊れる)

- **時間規約 exp(−iωt)**: 損失媒質は誘電率の**正**の虚部。負の虚部は利得になる。
  - Lorentz 分散 `eps = einf + ae²/(ce² − ω² − i·be·ω)` は正虚部を返す (正しい)。
  - 導電率項は `e += i·σ/(ωε₀)` (正)。符号を落とすと利得媒質になり T>1 が出る。
- **内部単位は μm**: `.orcwa` の座標はメートル指定で、`RCWAInput.cpp` が
  `M_TO_UM` で μm へ変換する。波長も μm。
- **エネルギー正規化**: `scatterPlaneWave` の REF/TRN 正規化は、単位振幅
  TE/TM 基底 (|E_inc_3D|²=1) を使う前提で `1/kzinc`。**偏波ベクトルと正規化は
  対で決まっており、片方だけ変えると斜め入射で R+T≠1 になる。**
- **偏波ベクトル** (RCWADriver): TM `(cosθcosφ, cosθsinφ)`, TE `(−sinφ, cosφ)`。
  法線入射 (θ≈0) では縮退するため x/y 方向を直接使う。
  任意偏波は複素係数 `(a_TM, a_TE)` で両基底を合成する
  (`|a_TM|²+|a_TE|²=1` を保てば単位振幅規格化が維持される)。
- **場の再構成には Bloch 因子が要る**: `f(r) = Σ f_mn·exp(i(kx0+2πm/Lx)x + …)`。
  `kx0_`/`ky0_` を落とすと |f| は正しいまま Re/Im だけが壊れる (見つけにくい)。
- **縦成分の導出** (コード内 H 規格化 `H_code = i/(ωε₀)·H_phys`):
  - `Hz = (i/k0²)(kx·Ey − ky·Ex)` (Faraday 則, μ=1)
  - `Ez = i·[[ε]]⁻¹(kx·Hy − ky·Hx)` (`[[ε]]` = ContinuousXY 畳み込み行列)
- **μ は P/Q に ε と対称に入る**:
  `P = [[Kx εzz⁻¹Ky, k0²μyy − Kx εzz⁻¹Kx], [Ky εzz⁻¹Ky − k0²μxx, −Ky εzz⁻¹Kx]]`、
  `Q` は ε↔μ を入れ替えた形 (全体の符号が反転)。半無限媒質は `k² = k0²εμ`、
  透過率の規格化には `μ_ref/μ_trn` が要る (`S_z ∝ kz/μ`)。
  `Layer::isMagnetic_` が false なら追加の Fourier 変換と SVD を省く高速経路。

## 構造上の前提

- **RCWA は周期境界のみ**: `pbc=0` は警告の上で無視される。
- 層スタックは「最上段 = 入射側半無限」「最下段 = 透過側半無限」。
  **両端のスラブは均質でなければならない** (格子層が端に来ると R+T が壊れる)。
- **材料インデックス**: 0=空気, 1=PEC (eps=1+1e8i), 2以降がユーザ定義
  (`material` 行の順)。複素誘電率は `material_eps = m epsr epsi` /
  `material_index = m n k` で直接指定できる。
- **一様層の固有モードは縮退しており、列インデックス ≠ ハーモニクス次数。**
  平面波励起は必ず `generateHorizontalPlaneWave` で係数化する
  (固有モード添字の直接指定は ±Kx が混ざり平面波にならない)。

## `.orcwa` 入力の要点

- **偏波指定**: `planewave = θ φ pol [psi]` — pol は 1=TM(p), 2=TE(s),
  3=psi[deg] の直線偏波 (0=TM, 90=TE), 4=右円偏波, 5=左円偏波。
- **磁性材料**: `material = type epsr esgm amur msgm` の amur/msgm が効く。
  `material_mu = m mur [mui]` でも直接指定できる。μ=1 の問題は高速経路を通る。
- **多極分散**: `material_dispersion` を同じ材料に複数行書くと極が**加算**される
  (`eps = einf + Σ_p ae_p²/(ce_p²−ω²−i·be_p·ω)`)。導電率項とも加算される。
- `.orcwa` は FDTD と共有する形式。未知キーワードは警告するが、`plot*` /
  `far1d*` / `far2d*` / `near*` と時間領域固有の設定は FDTD 専用として黙殺する
  (`isFDTDOnlyKeyword`)。FDTD 側にキーワードが増えても誤警告しないよう
  接頭辞判定にしてある。

## 数値安定性 (実際に踏んだもの — 外さないこと)

- 固有値ソルバーは `rcwa/MKLEigenSolver.hpp` の 3 分岐
  (MKL / LAPACKE / Eigen フォールバック)。**±m 次数の縮退固有値**で
  固有ベクトル基底が不良条件になり回折効率が発散する事象があるため、
  Eigen 経路には対角への非一様微小摂動 (~1e-11·‖A‖)、LAPACKE 経路には
  残差検算 ‖A·V−V·D‖ + Eigen 解き直しが入っている。**削除しない**。
- `EIGEN_DONT_PARALLELIZE` は OpenMP との競合回避。外さない。
- C99 VLA 禁止 (MSVC 対応)。libm は `MATH_LIB` 変数経由 (MSVC では空)。
  C++ ランタイムは CMake が自動でリンクするので `stdc++` を明示しない。
- MSVC は `/bigobj` 必須。Eigen のテンプレートを多用する翻訳単位
  (`rcwa/`, `tests/test_rcwa_input.cpp`) は obj のセクション数が既定上限を
  超えて C1128 になる。
- Windows CI は `-DWITH_RCWA_LEGACY_TESTS=OFF` (gdstk/core 依存の旧テストのみ
  無効)。`orcwa_rcwa` と `test_rcwa_input` は Windows でもビルド・実行される。
- **テストが作るテンポラリファイルは `tmpPath()` を通す** (`test_rcwa_input.cpp`)。
  Windows に `/tmp` は無い。`TMPDIR`/`TEMP`/`TMP` を見て `/tmp` にフォールバックする。
- **libc++ でしか出ないコンパイルエラーがある** (実績: Eigen の 1×1 `Product`
  → `std::complex` の暗黙変換は libstdc++ では通り libc++ では通らない)。
  Linux CI に `clang++ -stdlib=libc++ -fsyntax-only` の先行チェックを入れてある。
  ローカルで再現するには:
  ```bash
  sudo apt-get install -y clang libc++-dev libc++abi-dev
  clang++ -std=c++17 -stdlib=libc++ -fsyntax-only -I/usr/include/eigen3 -I. rcwa/*.cpp
  ```

## 落とし穴 (Gotchas)

- **Wood アノマリー**: λ が周期 Lx に一致すると回折次数がすれすれ (kz≈0)
  になり数値破綻する。テストの波長範囲は回折次数が確実に evanescent か
  伝搬になるよう選ぶ (例: Lx=0.5 μm なら λ=0.6–0.9 μm)。
- `SET_HARMONICS_2D` 相当のラムダは `RCWASolver.cpp` 内に**2 箇所**ほぼ同文で
  存在する。片方だけ直すと `saveFieldImage` のオーバーロード間で挙動が割れる。
  場の再構成ブロック (wavelet 合成) も同じく 2 箇所ある。
- **Eigen の `dot()` は第一引数を共役する**。ハーモニクス合成に使うと位相が
  反転するので `transpose() *` を用いる。
- **`harmonics2D` は (nx_ × ny_)**。YZ 断面で y ハーモニクスへ縮約するには
  `transpose()` が必要 (nx_≠ny_ だと次元不整合で落ちる)。
- 場 CSV は 8 桁精度 (`setprecision(8)`)。数値比較の許容誤差は 1e-6 程度に。
- テストで `.orcwa` を組み立てるとき、`ostringstream` の既定精度は 6 桁。
  厳密一致を見るケースでは `std::setprecision(17)` を明示する。

## Git 運用

- 指定された作業ブランチで開発し `git push -u origin <branch>` でプッシュする。
- 明示的に依頼されない限り Pull Request は作成しない。
- コミットメッセージは変更内容を具体的に (物理的な修正は理由も書く)。
