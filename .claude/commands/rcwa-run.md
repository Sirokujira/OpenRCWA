---
description: 指定した .orcwa 入力を RCWA ソルバで解いて結果を要約する
argument-hint: <input.orcwa> [orcwa_rcwa の追加オプション]
---

`.orcwa` 入力ファイルを RCWA ソルバで解いて結果を要約して。

対象: $ARGUMENTS

手順:

1. `bin/orcwa_rcwa` がなければ `/build-test` の手順でビルドする。
2. `./bin/orcwa_rcwa $ARGUMENTS` を実行する
   (場の断面が欲しいときは `-field Ez -slice xz -scoord 0 -sopt mod` などを付ける)。
3. 出力 CSV (`<input>_rcwa.csv`) を読み、波長ごとの R / T / A を表で要約する。
4. 妥当性を確認して報告する:
   - 無損失構造なら R+T≈1 (ずれていたら警告)
   - A<0 や T>1 は物理的に不正 — 入力の材料定義 (符号規約) を疑う
   - 掃引時はピーク/ディップの波長位置を挙げる
