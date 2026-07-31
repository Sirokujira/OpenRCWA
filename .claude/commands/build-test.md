---
description: RCWA をビルドして単体テストを実行し、結果を要約する
---

OpenRCWA の RCWA 系ターゲットをビルドし、単体テストを実行して結果を報告して。

手順:

1. `build/` が未構成なら `cmake -S . -B build -DBUILD_GPU=OFF` を実行する。
   BLAS/LAPACK が見つからない場合は依存パッケージ
   (libblas-dev liblapack-dev liblapacke-dev libeigen3-dev libhdf5-dev) を
   インストールしてから再構成する。
2. `cmake --build build --target orcwa_rcwa test_rcwa_input -j$(nproc)` でビルド。
3. `./bin/tests/test_rcwa_input` を実行。出力は
   `grep -v "compute interface\|begin to solve\|scatter matrices\|layer state"`
   でソルバ進捗ログを除去して確認する。
4. PASS/FAIL の内訳を要約する。FAIL があれば該当テストの物理的な意味
   (エネルギー保存・解析解比較など) と原因候補を説明する。
