---
paths:
  - "rcwa/**"
  - "tests/test_rcwa_input.cpp"
  - "src/rcwa_Main.cpp"
---

# RCWA 実装規約

## 物理規約 (違反すると結果が静かに壊れる)

- **時間規約 exp(−iωt)**: 損失は誘電率の**正**の虚部。
  - Lorentz 分散 `eps = einf + ae²/(ce² − ω² − i·be·ω)` は正虚部を返す (正)。
  - 導電率項は `e += i·σ/(ωε₀)` (正)。負にすると利得媒質になり T>1 が出る。
- **内部単位は μm**: `.orcwa` の座標はメートルで、`RCWAInput.cpp` が
  `M_TO_UM` で μm へ変換する。波長も μm。
- **エネルギー正規化**: `scatterPlaneWave` の REF/TRN 正規化は、単位振幅
  TE/TM 基底 (|E_inc_3D|²=1) を使う前提で `1/kzinc`。偏波ベクトルと正規化は
  対で決まっており、片方だけ変えると斜め入射で R+T≠1 になる。
- **偏波ベクトル** (RCWADriver): TM `(cosθcosφ, cosθsinφ)`, TE `(−sinφ, cosφ)`。
  法線入射 (θ≈0) では縮退するため x/y 方向を直接使う。
  任意偏波は複素係数 `(a_TM, a_TE)` で両基底を合成する
  (`|a_TM|²+|a_TE|²=1` を保てば単位振幅規格化が維持される)。
- **場の再構成には Bloch 因子が要る**: `f(r) = Σ f_mn·exp(i(kx0+2πm/Lx)x + …)`。
  `kx0_`/`ky0_` を落とすと |f| は正しいまま Re/Im だけが壊れる (見つけにくい)。
- **縦成分の導出** (コード内 H 規格化 H_code = i/(ωε₀)·H_phys):
  - `Hz = (i/k0²)(kx·Ey − ky·Ex)` (Faraday 則, μ=1)
  - `Ez = i·[[ε]]⁻¹(kx·Hy − ky·Hx)` ([[ε]] = ContinuousXY 畳み込み行列)

- **μ は P/Q に ε と対称に入る**:
  `P = [[Kx εzz⁻¹Ky, k0²μyy − Kx εzz⁻¹Kx], [Ky εzz⁻¹Ky − k0²μxx, −Ky εzz⁻¹Kx]]`,
  `Q` は ε↔μ を入れ替えた形 (全体の符号が反転)。半無限媒質は `k² = k0²εμ`、
  透過率の規格化には `μ_ref/μ_trn` が要る (S_z ∝ kz/μ)。
  `Layer::isMagnetic_` が false なら追加の Fourier 変換と SVD を省く高速経路。

## 構造上の前提

- 層スタックは「最上段 = 入射側半無限」「最下段 = 透過側半無限」。
  **両端のスラブは均質でなければならない** (格子層が端に来ると R+T が壊れる)。
- 材料インデックス: 0=空気, 1=PEC (eps=1+1e8i), 2以降がユーザ材料。
- 一様層の固有モードは縮退しており、列インデックス ≠ ハーモニクス次数。
  平面波励起は必ず `generateHorizontalPlaneWave` で係数化する
  (固有モード添字の直接指定は ±Kx が混ざり平面波にならない)。

## 落とし穴 (Gotchas)

- **Wood アノマリー**: λ が周期 Lx に一致すると回折次数がすれすれ (kz≈0)
  になり数値破綻する。テストの波長範囲は回折次数が確実に evanescent か
  伝搬になるよう選ぶ (例: Lx=0.5 μm なら λ=0.6–0.9 μm)。
- 場 CSV は 8 桁精度 (`setprecision(8)`)。数値比較の許容誤差は 1e-6 程度に。
- `SET_HARMONICS_2D` 相当のラムダは RCWASolver.cpp 内に**2 箇所**ほぼ同文で
  存在する。片方だけ直すと saveFieldImage のオーバーロード間で挙動が割れる。
  場の再構成ブロック (wavelet 合成) も同じく 2 箇所ある。
- **Eigen の `dot()` は第一引数を共役する**。ハーモニクス合成に使うと位相が
  反転するので `transpose() *` を用いる。
- **`harmonics2D` は (nx_ × ny_)**。YZ 断面で y ハーモニクスへ縮約するには
  `transpose()` が必要 (nx_≠ny_ だと次元不整合で落ちる)。
- テストで `.orcwa` を組み立てるとき、`ostringstream` の既定精度は 6 桁。
  厳密一致を見るケースでは `std::setprecision(17)` を明示する。

- `.orcwa` は FDTD と共有する形式。未知キーワードは警告するが、`plot*` /
  `far1d*` / `far2d*` / `near*` と時間領域固有の設定は FDTD 専用として黙殺する
  (`isFDTDOnlyKeyword`)。FDTD 側にキーワードが増えても誤警告しないよう
  接頭辞判定にしてある。

## テスト規約

- 物理変更には必ず解析解 (Fresnel / Brewster / エネルギー保存) との比較テストを
  `tests/test_rcwa_input.cpp` に追加する。
- テストは `CHECK` マクロで失敗を数え、`main` の戻り値 = 失敗数。
- 新キーワードはパース単体テスト + solve までの結合テストの両方を書く。
- **テンポラリファイルのパスは `tmpPath("name.csv")` を通す**。`/tmp/...` を
  直接書かないこと (Windows に `/tmp` は無く、`test_rcwa_input` は 3
  プラットフォームすべての CI で実行される)。
