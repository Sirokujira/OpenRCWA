#!/usr/bin/env python3
"""time_series_data.h5 の発熱量密度 P_loss を検証する (CI の回帰チェック用)。

P_loss は時間平均の消費電力密度 [W/m^3]:

    P = 1/2 sigma_e |E|^2 + 1/2 sigma_m |H|^2

導電率はセルごとの材料から引かなければならない。以前の実装は material_id = 0
(= 空気) の導電率を全セルに使い、磁気損失を mu'' = 1e-3 のマジック定数にして
いたため、PEC + 空気だけの無損失問題でも最大 2.5e11 W/m^3 を返していた。

使い方:

    # 無損失入力 (dipole.ofd): P_loss は全セル厳密に 0
    python3 ci/check_ploss.py lossless time_series_data.h5

    # 損失材あり (lossy_block.ofd): 損失が材料のある場所にだけ立つ
    python3 ci/check_ploss.py lossy --sigma-e 0.5 --sigma-m 0.3 time_series_data.h5

一致すれば exit 0、外れたら内容を表示して exit 1。
"""
import argparse
import sys

import h5py
import numpy as np

# 損失が「材料のある場所に限られている」と見なす非ゼロセルの割合の上限。
# 修正前の実装は場のある全セルが非ゼロになった (実測 91%)。
# 正しい実装では損失ブロックの体積ぶんだけ (同じ入力で実測 2.6%)。
LOCALIZED_MAX_FRACTION = 0.5


def data_groups(f):
    return sorted(k for k in f if k.startswith("data"))


def check_lossless(f):
    problems = []
    worst = 0.0
    groups = data_groups(f)
    if not groups:
        # ここを見ないと、data グループが 1 つも無いファイルが
        # 「max = 0」で素通りしてしまう (空振りで合格する)
        return ["data グループが 1 つも無い"], ""
    checked = 0
    for g in groups:
        name = g + "/P_loss"
        if name not in f:
            problems.append(f"{name} が無い")
            continue
        p = f[name][()]
        if p.size == 0:
            problems.append(f"{name} が空 (要素数 0)")
            continue
        checked += 1
        worst = max(worst, float(np.max(np.abs(p))))
    if checked == 0:
        problems.append("中身のある P_loss が 1 つも無い")
    if worst != 0.0:
        problems.append(
            f"無損失入力なのに P_loss が 0 でない (max|P_loss| = {worst:.6g})。"
            " セルごとの材料から導電率を引けていない可能性がある"
        )
    return problems, f"max|P_loss| = {worst:.6g} ({checked} グループを検査)"


def check_lossy(f, sigma_e, sigma_m):
    problems = []
    groups = data_groups(f)
    if not groups:
        return ["data グループが無い"], ""
    g = groups[-1]  # 場が最も育っている最終ステップで見る
    for name in (g + "/P_loss", g + "/E", g + "/H"):
        if name not in f:
            return [f"{name} が無い"], ""

    p = f[g + "/P_loss"][()]
    e = f[g + "/E"][()]
    h = f[g + "/H"][()]
    # [1, NFreq2, NN, 1] / [1, NFreq2, NN, 6]
    p = p.reshape(p.shape[1], p.shape[2])
    e = e.reshape(e.shape[1], e.shape[2], 6)
    h = h.reshape(h.shape[1], h.shape[2], 6)

    nn = p.shape[1]
    if p.size == 0:
        return [f"{g}/P_loss が空 (要素数 0)"], ""
    nonzero = int(np.count_nonzero(p))
    frac = nonzero / (p.size or 1)

    if float(p.min()) < 0.0:
        problems.append(f"P_loss が負のセルがある (min = {p.min():.6g})")
    if nonzero == 0:
        problems.append("損失材があるのに P_loss が全セル 0")
    if frac > LOCALIZED_MAX_FRACTION:
        problems.append(
            f"P_loss の非ゼロ率が {100 * frac:.1f}% と高すぎる"
            f" (>{100 * LOCALIZED_MAX_FRACTION:.0f}%)。損失が材料のある場所に"
            " 限られておらず、全セルに同じ導電率を使っている可能性がある"
        )

    # 同じファイルの E/H から作る物理的な上限。セルごとの導電率は入力に現れる
    # 最大値を超えないので、これを上回ったら材料の引き方が間違っている。
    e2 = np.sum(e * e, axis=2)
    h2 = np.sum(h * h, axis=2)
    upper = 0.5 * sigma_e * e2 + 0.5 * sigma_m * h2
    excess = p - upper
    tol = 1e-9 * np.maximum(upper, 1.0)
    bad = int(np.count_nonzero(excess > tol))
    if bad:
        problems.append(
            f"P_loss が E/H から決まる上限 (sigma_e={sigma_e}, sigma_m={sigma_m})"
            f" を超えるセルが {bad} 個ある (最大超過 {float(excess.max()):.6g})"
        )

    summary = (
        f"{g}: NN={nn} 非ゼロ {nonzero} セル ({100 * frac:.1f}%)"
        f" max={float(p.max()):.6g} 上限超過 {bad} セル"
    )
    return problems, summary


def main(argv):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("mode", choices=["lossless", "lossy"])
    ap.add_argument("h5file")
    ap.add_argument("--sigma-e", type=float, default=0.0,
                    help="入力に現れる電気導電率の最大値 [S/m] (lossy のみ)")
    ap.add_argument("--sigma-m", type=float, default=0.0,
                    help="入力に現れる磁気導電率の最大値 [1/Sm] (lossy のみ)")
    args = ap.parse_args(argv)

    with h5py.File(args.h5file, "r") as f:
        if args.mode == "lossless":
            problems, summary = check_lossless(f)
        else:
            problems, summary = check_lossy(f, args.sigma_e, args.sigma_m)

    if problems:
        print(f"FAILED ({args.mode}): {args.h5file}")
        for p in problems:
            print("  " + p)
        return 1
    print(f"OK ({args.mode}): {summary}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
