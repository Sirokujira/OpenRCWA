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
| `grating.ofd` | `orcwa` | 周期格子の波長掃引。全点で R+T=1 |
| `dipole.ofd` | `orcwa` | FDTD (RCWA ではない) |

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
| `/data%06d/` | 時間ステップごとの `E`, `H`, `P`, `P_loss`, `Surface` |

**RCWA 経路** (`sol/rcwa_hdf5.c` が出力。`.ofd` / `.orcwa` 共通):

| パス | 内容 |
|---|---|
| `/metadata/solver_mode` | `"RCWA"` — FDTD と区別する目印 |
| `/metadata/npol`, `pol_labels` | 偏波数と名前 (`.ofd` は 2 = TE/TM、`.orcwa` は 1) |
| `/metadata/NFreq`, `theta`, `phi`, `rcwa_harmonics`, `rcwa_period`, `rcwa_nlayer`, `Title` | 計算条件 |
| `/rcwa/spectrum` | compound **[npol][NFreq]** = `{frequency[Hz], lambda[m], R, T, A}` |

- **`spectrum` は 2 次元**なので、偏波ごとに 1 本の曲線として素直に読める。
- **波長は [m]**。`.orcwa` の内部単位は μm だが、書き出し時に m へ揃えてある。
- HDF5 は CSV より精度が高い (CSV は `%.8e` で丸めている)。
- 書き出しに失敗しても警告のみで続行する (結果は CSV にあるため)。
- 出力ファイル名: `.ofd` 経路は `time_series_data.h5`、`.orcwa` 経路は
  `-o` で指定した CSV の拡張子を `.h5` に置換したもの。

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
rcwa      = <N> <period[m]>
rcwalayer = <eps1r> <eps2r> <fill> <thickness[m]> [<eps1i> <eps2i>]
planewave = <theta[deg]> <phi[deg]> <pol>
```

- **誘電率の虚部は省略可** (省略時 0 = 無損失、従来の書式と後方互換)。
  損失は **正**の虚部 (exp(−iωt))。負は利得なので入力チェックで弾く。
  金属のように**実部が負**の材料も指定できる。
- **`planewave` の pol は無視される。** RCWA モードは常に TE/TM 両方を
  計算して CSV の 4 列に出力する (列構成を固定するため)。使うのは θ/φ のみ。
  `planewave` 自体を省略すると法線入射。θ は (−90, 90)。
- 偏波基底と Bloch 波数の決め方は `rcwa/RCWADriver.cpp` と**同一にすること**。
  正規化 (1/kzinc) と対で決まっており、片方だけ変えると斜め入射で R+T≠1 になる。
- **`.ofd` 経路は分散 (波長依存) を扱えない**。`rcwalayer` の誘電率は固定値
  なので、金属など分散の効く材料は設計波長近傍でしか意味を持たない。
  波長依存が必要なら `.orcwa` の `material_dispersion` + `orcwa_rcwa` を使う。
- **`.ofd` 経路は 1D・単一周期のみ**。2D 格子は `.orcwa` 経路が必要。

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
  解き直す。**環境によっては毎回この検算に落ちる** (実測: macOS で
  `LAPACKE geev unreliable (info=0)` が全ソルブで発生。結果は Eigen 経路で
  正しく、R+T=1 も成立する)。原因は LAPACK 本体と LAPACKE の提供元の
  食い違いが疑わしい (macOS の `find_package(LAPACK)` は Accelerate
  framework を拾う一方、LAPACKE は Homebrew の `liblapacke.dylib` になる)。
  切り分けには次を見る:
  ```bash
  grep -i lapacke build/CMakeCache.txt   # どのヘッダ/ライブラリを選んだか
  otool -L bin/orcwa | grep -i "lapack\|Accelerate"   # 実際のリンク先
  ```
  ログは最初の 1 回だけ相対残差つきで出る (毎ソルブ出すと大量になるため)。
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
