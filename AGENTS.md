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
./bin/tests/test_core_integrators  # core/ の積分器 (解析解・収束・積分器の一致, 約 40 秒)
ctest --test-dir build        # 上の 2 本 + .ofd 経路 (2D / 分散 / 次数別効率) の一式

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
- `core/` の積分器 (`RCWAIntegrator` / `DifferentialIntegrator` /
  `RedhefferIntegrator`) は `tests/test_core_integrators.cpp` で判定する:
  一様区間の S 行列 = `diag(exp(i·γ·L/k0))` (γ は k0² スケールの PQ の固有値の
  平方根、伝搬長は L/k0 に正規化 — `RCWAScatterMatrix::normalizedLz_`)、
  境界が動く区間での RCWA の 2 次収束 (区間数 2 倍で誤差約 1/4)、
  Differential / Redheffer と細かい RCWA の一致。`test_cpu` などの旧テストは
  正解データ作成に 1 日かかるベンチマークで判定を持たないので CI では回さない。
  旧テストと同じく gdstk/core 依存なので `WITH_RCWA_LEGACY_TESTS` (MSVC は OFF)。

## 実行

```bash
./bin/orcwa_rcwa [-o out.csv] [-v] \
    [-field Ez -slice xz -scoord 0 -sopt mod] [-device] input.orcwa
```

## サンプル入力 (`data/sample/`)

解析解が分かっているものを揃えてある。物理を変更したら**まずこれで確認**する。
CI (Linux ジョブ) が全件を実行し、エネルギー保存と解析解との一致を検証する。

| ファイル | 経路 | 内容 / 期待値 |
|---|---|---|
| `rcwa_fresnel.orcwa` | `orcwa_rcwa` | 空気/ガラス界面。R=0.04 (波長によらず一定) |
| `rcwa_brewster.orcwa` | `orcwa_rcwa` | ブリュースター角 56.30993° の TM 入射。R≈0 |
| `rcwa_pillar2d.orcwa` | `orcwa_rcwa` | 2D 正方格子の円柱 (メタサーフェス)。R+T=1 |
| `rcwa_metal_drude.orcwa` | `orcwa_rcwa` | Drude 金属薄膜。**A>0** (負なら符号規約破れ) |
| `ar_coating.ofd` | `orcwa` | 1/4 波長 AR コート。設計波長 550nm で R≈0 |
| `rcwa_oblique.ofd` | `orcwa` | 斜め入射 (ブリュースター角)。R_TM≈0, R_TE=0.147929 |
| `rcwa_metal.ofd` | `orcwa` | 複素誘電率の金属薄膜。R=0.909420, T=0.084430, A>0 |
| `rcwa_rect2d.ofd` | `orcwa` | 2D 格子 (正方格子の正方形柱)。法線入射で TE = TM、R+T=1 |
| `rcwa_drude_film.ofd` | `orcwa` | 分散材料 (`rcwadispersion`) の Drude 金属薄膜。R+T<1, A>0。`rcwa_metal_drude.orcwa` と同じ構造 |
| `grating.ofd` | `orcwa` | 周期格子の波長掃引。全点で R+T=1 |
| `dipole.ofd` | `orcwa` | FDTD (RCWA ではない)。無損失 (PEC+空気) なので **P_loss は全セル厳密に 0** |
| `lossy_block.ofd` | `orcwa` | FDTD。電気損失材と磁気損失材のブロック。P_loss がその場所にだけ立つ |
| `pec_sphere_rcs.ofd` | `orcwa` | FDTD 平面波入射。半径 0.1 m の PEC 球の後方散乱断面積が Mie 級数と一致 (`ci/check_mie.py`) |

**FDTD 平面波 (`planewave`) の確認は `pec_sphere_rcs.ofd` で行う。** 球の
後方散乱は入射方向・偏波によらないので、`planewave` の θ/φ/pol を変えても
同じ解析解と比べられる (CI は直列で法線入射 V 偏波、MPI n=3 で
θ=60°, φ=30° の H 偏波を回す)。直交格子の階段近似で球は実効的に約 3% 太るため、
`check_mie.py` は半径を ±5% 内で合わせ込んでからピーク比 6% で判定する。
実測: 正しい実装 3.0%、平面波の基準位置 r0 の z を 0.5 μm に固定上書きしていた
旧実装 11.9% (しかも収束せず max steps まで回る)。r0 は「領域中心から上流側へ
外接球半径 + λ0/2」で、t=0 で入射パルスが領域に入っていないことが前提
(散乱界定式化の初期条件)。触るときはこのサンプルで確認すること。

**格子では TM の収束が遅い。** 誘電率境界で法線 E が不連続になるためで、
RCWA の既知の性質。`grating.ofd` の実測 (λ=461nm):

| N | 次数 | R_TE | R_TM |
|---|---|---|---|
| 5 | 11 | 0.138130 | 0.204841 |
| 12 | 25 | 0.138187 | 0.206747 |
| 20 | 41 | 0.138193 | 0.206805 |

TE は N=5 でも 1e-4 精度だが、TM は N=5 だと収束値から約 1% ずれる。
このため `grating.ofd` の既定は **N=12** にしてある (N=20 との差は
R_TE 9e-6 / R_TM 8.7e-5)。**新しい格子サンプルを足すときは、TE だけで
なく TM でも N を振って収束を確認すること。** 一様層 (AR コート等) は
回折が無いので低次で足りる。

**RCWA モードの入力に `orcwa_post` は使えない。** `rcwalayer` を含む `.ofd` を
与えると `orcwa` は RCWA コアに分岐して `rcwa_efficiency.csv` だけを出力し、
時間領域の `orcwa.out` を作らない。`orcwa_post` は `rcwalayer` を検出したら
「ポスト処理なし」と表示して**正常終了 (exit 0)** する — ソルバ後に必ず post を
呼ぶ GUI/スクリプトの流れを壊さないための仕様。エラーにしてはいけない。

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

## HDF5 出力 (GUI 表示用)

`time_series_data.h5` は **FDTD と RCWA で中身が違う**。GUI 側は
`/metadata/solver_mode` の**有無**で判別する (RCWA のときだけ存在し、値は
`"RCWA"`)。CSV / `orcwa.out` は従来どおり出力されるので、既存の読み込みは
壊れない。

**FDTD 経路** (`sol/solve.c` が出力):

| パス | 内容 |
|---|---|
| `/metadata/` | `Dt`, `VFeed`, `IFeed`, `VPoint`, `Eiter`, `Hiter`, メッシュ (`Xn`/`Xc` 等), `input_impedance`, `s_parameters` ほか |
| `/data%06d/` | 時間ステップごとの `E`, `H`, `P`, `P_loss` |

**RCWA 経路** (`sol/rcwa_hdf5.c` が出力。`.ofd` / `.orcwa` 共通):

| パス | 内容 |
|---|---|
| `/metadata/solver_mode` | `"RCWA"` — FDTD と区別する目印 |
| `/metadata/npol`, `pol_labels` | 偏波数と名前 (`.ofd` は 2 = TE/TM、`.orcwa` は 1) |
| `/metadata/NFreq`, `theta`, `phi`, `rcwa_harmonics`, `rcwa_period`, `rcwa_harmonics_y`, `rcwa_period_y`, `rcwa_nlayer`, `Title` | 計算条件 (`_y` は 1D 格子なら 0) |
| `/rcwa/spectrum` | compound **[npol][NFreq]** = `{frequency[Hz], lambda[m], R, T, A}` |
| `/rcwa/R_orders`, `/rcwa/T_orders` | double **[npol][NFreq][2Nx+1][2Ny+1]** = 回折次数ごとの効率 (`.ofd` 経路のみ。エバネッセント次数は 0) |
| `/rcwa/order_m`, `/rcwa/order_n` | 各軸の次数 (−Nx..Nx / −Ny..Ny) |

- **`spectrum` は 2 次元**なので、偏波ごとに 1 本の曲線として素直に読める。
- **波長は [m]**。`.orcwa` の内部単位は μm だが、書き出し時に m へ揃えてある。
- HDF5 は CSV より精度が高い (CSV は `%.8e` で丸めている)。
- 書き出しに失敗しても警告のみで続行する (結果は CSV にあるため)。
- 出力ファイル名: `.ofd` 経路は `time_series_data.h5`、`.orcwa` 経路は
  `-o` で指定した CSV の拡張子を `.h5` に置換したもの。

### ビルド構成による出力の違い (揃えること)

GUI は `time_series_data.h5` を読むので、**ビルド構成でデータセットが欠けると
機能が黙って落ちる**。4 つのソルバ (`sol/solve.c` / `mpi/solve.c` /
`cuda/solve.cu` / `cuda_mpi/solve.cu`) は同じデータセットを出す:

| グループ | データセット |
|---|---|
| `/data%06d` | `E`, `H`, `P`, `P_loss` |
| `/metadata` | 格子・時刻・波形などのスカラ/配列, `Reflection` |

`Reflection` はセルごとの屈折率 `sqrt(Material[iEx[nn]].epsr)` を `[NN, 1]` で
出したもの。以前は CUDA 版だけが持っていた。

`/data%06d` の E/H/P/P_loss の書き出しは、直列版と CUDA 版 (どちらもホストに
全領域の配列がある) で `sol/h5cells.c` の `write_cell_fields_hdf5()` を共有する。
**周波数ごとに 1 回の `H5Dwrite`** で書くこと。以前は 1 セルずつハイパースラブを
選んで書いており (出力 1 回あたり NFreq2 × NN × 4 回)、`pec_sphere_rcs.ofd`
では計算時間の 85% を占めていた (直列 96 秒 → 15 秒)。作業バッファは 1 周波数
ぶん (NN × 6 個の double) に留め、全周波数を一度に詰めない (大規模問題で
DFT 配列と同程度のメモリを余計に食うため)。MPI 版は rank ごとの集団書き込みで別実装。

**遠方界用の近傍界サーフェス (`Surface`) は HDF5 に出さない**。以前は
`/data%06d/Surface` と `/metadata/NSurface`・`/metadata/Surface` を書いていたが、
`NSurface` と `SurfaceEx` 等は `solve()` の**後**の `outputChars()`
(`alloc_farfield()` / `setup_farfield()`) で初めて設定されるため、HDF5 を書く
`solve()` の中では常に `NSurface = 0` で、全ビルドで中身が全要素 0 だった
(per-step 側は次元も `NSurface` ではなく `NN`)。GUI (OpenFDTD-X) も
`orcwa_post` も読んでいない。並列 HDF5 は領域を先に確保するため、MPI 版では
このゼロ領域が dipole でファイルの 43% (42.6 MB / 99.8 MB) を占めていた。
遠方界は従来どおり `orcwa.out` (`writeout`) 経由で `orcwa_post` が計算する。
HDF5 から遠方界を出したくなったら、`outputChars()` の後に最終 DFT から
`[NFreq2][NSurface]` で書くこと (時間ステップごとに書く意味は無い)。

**新しいデータセットを足すときは 4 本すべてに入れること。** 直列と MPI の
一致は `ci/compare_h5.py` が CI で検証する (CUDA 系は GPU が要るため CI では
ビルドのみ)。

### 出力先のパス (現状仕様 — 変更しないこと)

`.ofd` 経路の出力 (`rcwa_efficiency.csv` / `time_series_data.h5` /
`orcwa.log`) は**カレントディレクトリ基準の相対パス**で、入力ファイルの
場所からは導出しない。FDTD 経路の `orcwa.out` も同じ規則なので、両経路で
挙動は揃っている。

GUI (OpenFDTD-X) は入力ファイルのあるディレクトリへ `cd` してから
絶対パスで入力を渡すため、結果として出力は入力の隣に出る。

- **`orcwa` の `-out` は RCWA モードでは無視される。** `rcwa_run(fp_log)` が
  出力先を受け取らず、ファイル名がハードコードされているため
  (`sol/rcwa_bridge.cpp`)。FDTD では `-out` が `orcwa.out` に効くので、
  経路によって挙動が違う点に注意。
- これは既知かつ意図的な現状仕様 (2026-08 時点でユーザー判断により維持)。
  GUI から出力先を制御する必要が出たときに、`-out` を RCWA にも
  効かせるかを改めて判断する。

## 発熱量密度 `P_loss` (FDTD)

時間平均の消費電力密度 [W/m^3]:

    P = 1/2 sigma_e |E|^2 + 1/2 sigma_m |H|^2

- **導電率はセルごとの材料から引く**。Yee 格子では E/H の各成分が別の位置にあり、
  成分ごとに材料 ID (`iEx`/`iEy`/`iEz`, `iHx`/`iHy`/`iHz`) が違うので、
  成分ごとに `Material[id].esgm` / `Material[id].msgm` を掛ける。
  実装は `sol/powerloss.c` の `calculatePowerLoss()` に 1 本化してある。
- 空気 (ID 0) も PEC (ID 1) も `esgm = msgm = 0`。PEC は `C1 = C2 = 0` で E が 0 に
  固定されるため、完全導体は電力を消費しないモデルになる。
- 係数 1/2 は `cEx_r` 等が波高値 (peak) の位相子であることを前提にする。
- **旧実装の壊れ方** (2026-08 に修正): `material_id = 0` (= 空気) の導電率を全セルに
  使い、磁気損失を `mu'' = 1e-3` のマジック定数にしていた。そのため無損失の
  `dipole.ofd` でも最大 2.5e11 W/m^3 を返していた。さらに結果を捨てるだけの
  `updateTemperature()` (温度 T はどこにも出力されない) を毎ステップ回していた。
- 検証は `ci/check_ploss.py`:

```bash
# 無損失入力 (dipole.ofd): 全セル厳密に 0
python3 ci/check_ploss.py lossless time_series_data.h5
# 損失材あり (lossy_block.ofd): 材料のある場所にだけ立つ / E・H から決まる上限を超えない
python3 ci/check_ploss.py lossy --sigma-e 0.5 --sigma-m 0.3 time_series_data.h5
```

- `P_loss` は**直列 / MPI / CUDA / CUDA+MPI の 4 ビルドすべてが出力する**。
  実装は `sol/powerloss.c` の 1 本で、MPI 系は各 rank が自分の担当範囲を
  計算して集団書き込みする。

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

## `.ofd` の RCWA モード (GUI 経路) の要点

`.orcwa` と `.ofd` は**別形式**。拡張子を変えても変換されない。GUI
(OpenFDTD-X) が開けるのは `.ofd` だけで、RCWA として解くには
`rcwa` / `rcwalayer` キーが要る (`.orcwa` の `geometry` 系は読まれない)。

```
rcwa      = <Nx> <period_x[m]> [<Ny> <period_y[m]>]
rcwalayer = <eps1r> <eps2r> <fill> <thickness[m]> [<eps1i> <eps2i> [<fill_y>]]
rcwaincidence = <theta[deg]> <phi[deg]>            # RCWA の入射角 (planewave より優先)
planewave = <theta[deg]> <phi[deg]> <pol>           # 従来の書き方 (rcwaincidence が無いとき)
```

- **2D 格子**: `rcwa` に `Ny period_y` を足すと 2D になる (省略時は従来の 1D)。
  各層は「1 周期内に eps1 の矩形 `[0, fill·Λx] × [0, fill_y·Λy]` が 1 つ、残りが
  eps2」。`fill_y` の既定は 1 (y 方向に一様 = 1D 線格子) で、書くときは虚部の組
  (無損失なら `0 0`) を先に書く。1D で `fill_y < 1` は入力エラー。
  `fill` / `fill_y` は **[0, 1]** (端の値は一様層: 0 → eps2 だけ、1 → eps1 だけ)。
  GUI (OpenFDTD-X) は一様層を fill 0 / 1 で書くことがあり、以前の 1D は
  fill ≤ 0 / ≥ 1 を入力エラーにしていたため GUI の入力がカーネルで止まっていた。
- **層ごとに fill が違ってよい**。ブリッジは全層の fill の境界を合わせた共通の
  セル格子で各層を作る (`commonEdges()`)。`RCWASolver::addLayer()` は全層の
  セル境界が一致しないと層を黙って捨てるため、以前は格子層の fill を一様層と
  違う値にするだけで orcwa が Segmentation fault で落ちていた (既存サンプルは
  全層 fill=0.5 だったので見つからなかった)。同じ値のセルへの細分は
  Fourier 係数を変えないので結果は変わらない。
- 2D の検証は `ci/check_rcwa_ofd2d.py` (1D との一致・C4 対称・90° 回転・
  `.orcwa` 経路との一致・fill 混在)。CI の Linux / macOS / Windows で実行する。
- **回折次数ごとの効率**は `rcwa_orders.csv` (`frequency[Hz],lambda[m],pol,m,n,R,T`、
  伝搬する次数だけ) と HDF5 の `/rcwa/R_orders`・`T_orders` に出る。
  `rcwa_efficiency.csv` (GUI が読む 6 列) は変えないこと。次数の並びは
  `RCWASolver::evaluateAlphaBeta()` の `k = (2Ny+1)·(m+Nx) + (n+Ny)`
  (kx = kx0 + 2πm/Λx)。伝搬判定は**格子方程式** (kx²+ky² < k0²·Re(ε) を入射側 /
  透過側で別々に) で行い、伝搬しない側は 0 に揃える。**値が 0 かどうかで判定しない
  こと**: Linux ではエバネッセント次数が厳密に 0 になるが、macOS (libc++) は複素
  平方根の丸めで 1e-33 程度が残り、全次数が「伝搬」扱いになって CI で落ちた。検証は `ci/check_rcwa_ofd_orders.py`
  (和 = 合計、格子方程式、法線入射の ±m 対称、斜め入射での m / n の符号)。

- **誘電率の虚部は省略可** (省略時 0 = 無損失、従来の書式と後方互換)。
  損失は **正**の虚部 (exp(−iωt))。負は利得なので入力チェックで弾く。
  金属のように**実部が負**の材料も指定できる。
- **入射角は `rcwaincidence = θ φ`、無ければ `planewave` の θ φ** (どちらも無ければ
  法線入射)。θ は (−90, 90)。θ は法線からの角度で、入射側の横波数は
  `k0·n_inc·sinθ·(cosφ, sinφ)` (θ > 0, φ = 0 の反射 1 次は m = −1)。
  **`rcwaincidence` を足した理由**: `planewave` は FDTD の入射方向 (波が**来る**向き、
  横波数は −sinθ·(cosφ, sinφ)。OpenFDTD-X の既定は θ = 90°) と同じ行で、意味が
  違う。GUI で平面波を有効にしたまま RCWA を選ぶと θ = 90° で
  「theta must be within (-90, 90)」になっていた。GUI は RCWA の角度を
  `rcwaincidence` で書き、FDTD の `planewave` には触れない。rcwalayer の無い入力に
  `rcwaincidence` があれば入力エラー。検証は `ci/check_rcwa_ofd_incidence.py`
  (Fresnel / Brewster、`planewave` と一致、優先順位、コニカル入射の R+T=1、
  鏡映で m → −m、入力エラー)。
- **`planewave` の pol は無視される。** RCWA モードは常に TE/TM 両方を
  計算して CSV の 4 列に出力する (列構成を固定するため)。
- 偏波基底と Bloch 波数の決め方は `rcwa/RCWADriver.cpp` と**同一にすること**。
  正規化 (1/kzinc) と対で決まっており、片方だけ変えると斜め入射で R+T≠1 になる。
- **分散材料**: `rcwadispersion = <layer(1..)> <side(1|2)> <einf> <ae> <be> <ce>`
  [rad/s] で層の eps1 (side=1) / eps2 (side=2) に Lorentz 極を与える
  (`eps = einf + Σ ae²/(ce²−ω²−i·be·ω)`)。意味付けは `.orcwa` の
  `material_dispersion` と同じ: 同じ (layer, side) の複数行は極の**加算**、
  einf は最後の行、ae≠0 の極が 1 つでもあれば `rcwalayer` の静的な値は使わない。
  be>0 が損失、be<0 (利得) は入力エラー。一様層は eps1 / eps2 の両方に書く。
  Drude 金属は `ce=0, ae=ωp, be=γ`。検証は `ci/check_rcwa_ofd_dispersion.py`
  (Drude 薄膜を特性行列法と比べて 5e-10、`rcwa_metal_drude.orcwa` と 4e-9)。
- **`.ofd` 経路の層形状は「矩形 1 つ」まで**。円柱など任意形状の 2D 格子は
  `.orcwa` 経路 (`geometry`) が必要。

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
- LAPACKE 経路は毎回 ‖A・V−V・D‖∞/‖A‖∞ ≤ 1e-8 を検算し、外れたら Eigen で
  解き直す。**macOS/arm64 (Homebrew LAPACKE) ではこの検算に落ちる事象を
  CI で実測済み**: 相対残差 1.96〜3.60 と **O(1)** で、閾値がわずかに
  厳しいのではなく LAPACKE の返す固有ベクトルが完全に壊れている
  (結果は Eigen 解き直しで正しくなり、R+T=1 も成立する)。
  - 対策 1: `geev` は `LAPACK_COL_MAJOR` で呼ぶ (詰め替え時に自前で転置)。
    row-major ラッパの内部転置を踏まず、余分なコピーも消える。
    Linux では row-major 時代と出力がビット単位一致することを確認済み。
  - 対策 2: 検算に一度落ちたらプロセス内で LAPACKE を諦め、以後は直接
    Eigen で解く (粘着フォールバック)。壊れた zgeev を毎ソルブ払わない。
    Eigen 経路は対角摂動つきで Windows が常用しており正しさは同等。
  - ログは相対残差つきで初回のみ。O(1) なら環境の LAPACKE が壊れている、
    1e-8 をわずかに超える程度なら閾値の問題、と読み分ける。
  - **結論 (CI 実測 2026-08)**: col-major 化しても残差は row-major 時代と
    **全く同じ値** (1.96334 / 3.59536) で落ちた。壊れ方はラッパ層と無関係の
    決定的なもので、**Homebrew の LAPACKE/LAPACK 本体が macOS/arm64 で
    壊れている**。したがって macOS は実効的に Eigen 経路で動いている
    (粘着フォールバックにより無駄はプロセスあたり 1 ソルブのみ)。
    Homebrew lapack のバージョンが上がったら、build-macos ログの
    `unreliable` の有無で再評価する。CMake 構成をいじって直そうと
    しないこと — リンク構成 (Accelerate との混在) は原因ではない。
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

## MPI (`orcwa_mpi`) の並列 HDF5

```bash
sudo apt-get install -y openmpi-bin libopenmpi-dev libhdf5-openmpi-dev
cmake -B build-mpi -DCMAKE_BUILD_TYPE=Release -DWITH_CUDA=OFF -DWITH_MPI=ON \
      -DHDF5_PREFER_PARALLEL=ON -DWITH_RCWA=OFF -DWITH_RCWA_LEGACY_TESTS=OFF
cmake --build build-mpi -j"$(nproc)"   # → bin/orcwa, bin/orcwa_mpi, bin/orcwa_post

# 直列版と一致するか (不均等分割を必ず含める)
mpirun --oversubscribe -n 4 bin/orcwa_mpi dipole.ofd
python3 ci/compare_h5.py <serial>/time_series_data.h5 time_series_data.h5
```

- **`WITH_MPI=ON` には並列 HDF5 が要る**。直列 HDF5 だと `H5Pset_fapl_mpio`
  等が undefined reference になるため、configure で FATAL_ERROR にしてある。
- **`orcwa` (CPU) は CUDA / MPI の有無によらず常にビルドする。** RCWA 入力を
  解けるのは RCWA コアに結線された `orcwa` だけで、`orcwa_mpi` / `orcwa_cuda` /
  `orcwa_cuda_mpi` は「orcwa を使え」と案内して止まる。以前は `orcwa` を
  CUDA でも MPI でもない構成でしか作らなかったため、MPI/CUDA ビルドには
  案内先の `orcwa` 自体が無かった。
- **`_MPI` は MPI 版のターゲット (`orcwa_mpi` / `orcwa_cuda_libmpi` /
  `orcwa_cuda_mpi`) にだけ付ける。** 全体定義にすると、同じビルドの直列
  `orcwa` まで `_MPI` 付きでコンパイルされ、`sol/` に `#ifdef _MPI` を足した
  途端に MPI 構成の `orcwa` だけが壊れる。CI の build-mpi がフラグを直接検査する。
- **並列 HDF5 の `hdf5.h` は `mpi.h` を取り込む** (`H5public.h`) のに、
  FindHDF5 の `HDF5::HDF5` は MPI のパスを運ばない。そのため並列 HDF5 では
  `HDF5::HDF5` に `MPI::MPI_C` を INTERFACE で足している (HDF5 を使う全ターゲット
  に効く)。以前は `WITH_MPI=ON` で全ターゲットをビルドすると `orcwa_post` が
  必ず `mpi.h: No such file or directory` で落ちていた (CI は `--target orcwa_mpi`
  しかビルドしていなかったので気づけなかった)。直列プログラムは `MPI_Init` を
  呼ばないので直列のまま動く。
- **`H5Gcreate` / `H5Dcreate` は集団操作**。`if (commRank == 0)` で囲うと
  他 rank が `H5Fclose` (これも集団) で待ち続けてデッドロックする。
  生の書き込み (`H5Dwrite`) だけは INDEPENDENT 転送で rank 0 に限定してよい。
- **集団 `H5Dcreate` はデータ空間の次元が全 rank で一致していなければならない**。
  ここが実際に踏んだ最大の罠で、rank ローカルな値を次元に使うと
  **均等分割では顕在化せず不均等分割でだけハングする**。実測:
  - (現在は削除済みの) `Surface` データセットの第 3 次元に rank ローカルの `NN`
    (`sol/setupSize.c` の `iMax/jMax/kMax` 依存) を使っていた
    → n=2,3,5,6 は通り **n=4 (7,7,7,9) / n=7 (4,4,4,4,4,4,6) でハング**。
    gdb で見ると 3 rank が `H5FD_truncate`→`MPI_File_set_size`、
    1 rank が `H5AC__run_sync_point`→`MPI_Barrier` で別々の
    コミュニケータの集団操作を待っていた。→ 全体側の `g_NN` に修正。
  - `metadata/Gline` の次元に `NGline` を使っていた。`Gline`/`NGline` は
    `comm_broadcast` の対象外で rank 0 以外は 0 → データセットが 0 要素に
    なっていた。→ `arrays[]` の各 size を rank 0 から `MPI_Bcast` して揃える。
  - 新しいデータセットを足すときは **次元に使う変数が全 rank で同じか**を
    必ず確認すること。
- 場データ (E/H/P) は各 rank が自分の担当セルをグローバル添字のハイパー
  スラブに置き、`H5FD_MPIO_COLLECTIVE` で 1 回にまとめて書く
  (1 セルずつ書くと遅いうえ、呼び出し回数が rank ごとに変わりやすい)。
- **給電点・観測点の波形はそのセルを持つ rank にしか溜まらない**。
  出力前に `comm_feed()` / `comm_point()` で rank 0 へ集める
  (これを呼び忘れていて、rank 0 が給電セルを持たない n≥4 で `IFeed` が
  直列と 6.6e-3 ずれていた)。
- CI の `build-mpi` ジョブは**全ターゲットをビルド**し (`orcwa` / `orcwa_mpi` /
  `orcwa_post` が揃うこと、`_MPI` が `orcwa_mpi` にだけ付いていること)、
  MPI 構成の `orcwa` で RCWA (grating, R+T=1) を解き、その `orcwa` の FDTD 出力が
  純粋な直列構成の `orcwa` と一致することを確認する。
- CI の `build-mpi` ジョブが n=2,3,4,7 を **無損失 (`dipole.ofd`) と損失材あり
  (`lossy_block.ofd`) の両方**で回し、`ci/compare_h5.py` で直列版との一致
  (E/H/P/P_loss/Reflection とメタデータは差 0、Eiter/Hiter のみ
  総和順序ぶんの許容差) を確認する。無損失だけだと `P_loss` が全セル 0 で、
  MPI 側がまるごと壊れていても一致してしまうため両方回す。

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
