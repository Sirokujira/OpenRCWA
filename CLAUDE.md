# CLAUDE.md — OpenRCWA

OpenRCWA は周期構造の電磁界解析ソルバ集。RCWA (厳密結合波解析) エンジンと
OpenFDTD 由来の FDTD ソルバ (C / Python / CUDA) を含み、`.orcwa` 入力形式を共有する。

## ビルド

```bash
# 依存パッケージ (Ubuntu)
sudo apt-get install -y libblas-dev liblapack-dev liblapacke-dev \
                        libeigen3-dev libhdf5-dev

# 構成 (MKL/TBB/GPU はなくてもよい — 自動で Eigen ソルバにフォールバック)
cmake -S . -B build -DBUILD_GPU=OFF

# 主要ターゲット
cmake --build build --target orcwa_rcwa test_rcwa_input -j$(nproc)
```

生成物は `bin/` (実行ファイル) と `bin/tests/` (テスト) に置かれる。

## テスト

```bash
./bin/tests/test_rcwa_input   # RCWA 単体テスト (パーサ/ドライバ/場出力)
```

- 終了コード = 失敗数。`All tests PASSED` が成功の目印。
- ソルバの進捗ログが大量に出るため、確認時は
  `grep -v "compute interface\|begin to solve\|scatter matrices\|layer state"` で除去する。
- 物理変更をしたら必ずエネルギー保存 (無損失で R+T≈1) のテストを通すこと。
- CI は `.github/workflows/ci.yml` (ビルド + スモーク + 単体テスト)。

## 実行

```bash
./bin/orcwa_rcwa [-o out.csv] [-v] \
    [-field Ez -slice xz -scoord 0 -sopt mod] [-device] input.orcwa
```

## ディレクトリ構成

| パス | 内容 |
|---|---|
| `rcwa/` | RCWA エンジン本体 (`RCWASolver`, `Layer`, Fourier ソルバ) と `.orcwa` パーサ (`RCWAInput`) / ドライバ (`RCWADriver`) |
| `src/` | 実行ファイルの main (`rcwa_Main.cpp` = orcwa_rcwa, ほか FDTD 系) |
| `core/` | 積分器・レイヤサンプラ (RCWA 応用層) |
| `sol/`, `include/` | FDTD ソルバ (C) |
| `python/` | FDTD ソルバ (Python + Numba) |
| `cuda/`, `cuda_mpi/`, `gpu/` | GPU 実装 |
| `gdstk/` | GDSII 読み込みライブラリ (ベンダリング) |
| `tests/` | 単体テスト (`test_rcwa_input.cpp` が RCWA 系の主テスト) |
| `ci/` | CI 用スモークテスト入力 |

## 重要な物理・実装規約

詳細は `.claude/rules/rcwa-conventions.md` を参照 (rcwa/ 配下を触ると自動ロード)。
要点:

- **時間規約は exp(−iωt)**: 損失媒質は誘電率の**正**の虚部。負の虚部は利得になる。
- **内部単位は μm**: `.orcwa` はメートル指定で、パース時に μm へ変換される。
- **材料インデックス**: 0=空気, 1=PEC, 2以降がユーザ定義 (`material` 行の順)。
  複素誘電率は `material_eps = m epsr epsi` / `material_index = m n k` で直接指定できる。
- **RCWA は周期境界のみ**: `pbc=0` は警告の上で無視される。
- **偏波指定**: `planewave = θ φ pol [psi]` — pol は 1=TM(p), 2=TE(s),
  3=psi[deg] の直線偏波 (0=TM, 90=TE), 4=右円偏波, 5=左円偏波。
- **磁性材料**: `material = type epsr esgm amur msgm` の amur/msgm が効く。
  `material_mu = m mur [mui]` でも直接指定できる。μ=1 の問題は高速経路を通る。
- **多極分散**: `material_dispersion` を同じ材料に複数行書くと極が**加算**される
  (`eps = einf + Σ_p ae_p²/(ce_p²−ω²−i·be_p·ω)`)。導電率項とも加算される。

## Git 運用

- 指定された作業ブランチで開発し `git push -u origin <branch>` でプッシュする。
- 明示的に依頼されない限り Pull Request は作成しない。
- コミットメッセージは変更内容を具体的に (物理的な修正は理由も書く)。
