#!/usr/bin/env python3
"""PEC 球の後方散乱断面積を Mie 級数と比較する (FDTD 平面波入射の回帰チェック)。

orcwa.log の "=== cross section ===" 表 (frequency[Hz] backward[m*m] forward[m*m])
の backward 列を、完全導体球の厳密解と比べる。

    python3 ci/check_mie.py --radius 0.1 orcwa.log

Mie 級数 (Bohren & Huffman, m -> ∞ の極限):
    a_n = [x j_n(x)]' / [x h_n(x)]',  b_n = j_n(x) / h_n(x)
    sigma_back = (pi / k^2) |sum_n (2n+1) (-1)^n (a_n - b_n)|^2

判定: 直交格子の階段近似では球が実効的に少し太る (a/Δ = 10 で約 +3%)。
そこで半径を幾何半径の ±fit 内で動かし、帯域内の誤差
    max_f |sigma_FDTD - sigma_Mie(a')| / max_f sigma_Mie(a')
が最小になる a' を選び、その誤差が tol 以下なら合格とする。dB 比較にしないのは、
Mie 解の谷 (ka ≈ 1.7) で dB 誤差が過大に出て判定が不安定になるため。

平面波の基準位置 r0 を誤る (t=0 で入射パルスが領域内にある) と、散乱界に
初期過渡が乗ってスペクトルの形が崩れ、どの a' でも合わない
(実測: 正しい実装 3.1%、r0 を 0.5 μm に固定した旧実装 12.1%)。

一致すれば exit 0、外れれば内容を表示して exit 1。
"""
import argparse
import math
import sys

import numpy as np

C = 2.99792458e8


def sph_jy(nmax, x):
    """球ベッセル j_n, y_n (n = 0..nmax)。j は下向き漸化式で安定に求める。"""
    y = np.zeros(nmax + 2)
    y[0] = -math.cos(x) / x
    y[1] = -math.cos(x) / x**2 - math.sin(x) / x
    for n in range(1, nmax + 1):
        y[n + 1] = (2 * n + 1) / x * y[n] - y[n - 1]
    start = nmax + 30 + int(x)
    j = np.zeros(start + 2)
    j[start + 1], j[start] = 0.0, 1e-300
    for n in range(start, 0, -1):
        j[n - 1] = (2 * n + 1) / x * j[n] - j[n + 1]
    j *= (math.sin(x) / x) / j[0]
    return j[:nmax + 2], y


def pec_backscatter(a, f):
    k = 2 * math.pi * f / C
    x = k * a
    nmax = int(x + 4 * x ** (1 / 3) + 10)
    j, y = sph_jy(nmax, x)
    h = j + 1j * y
    s = 0j
    for n in range(1, nmax + 1):
        # [x z_n(x)]' = x z_{n-1}(x) - n z_n(x)
        dxj = x * j[n - 1] - n * j[n]
        dxh = x * h[n - 1] - n * h[n]
        an = dxj / dxh
        bn = j[n] / h[n]
        s += (2 * n + 1) * (-1) ** n * (an - bn)
    return math.pi / k**2 * abs(s) ** 2


def read_cross_section(path):
    rows, on = [], False
    with open(path, encoding="utf-8", errors="replace") as fp:
        for line in fp:
            if line.startswith("=== cross section"):
                on = True
                continue
            if on:
                if line.startswith("==="):
                    break
                t = line.split()
                try:
                    rows.append((float(t[0]), float(t[1])))
                except (ValueError, IndexError):
                    pass
    return rows


def max_err(rows, a):
    ref = [pec_backscatter(a, f) for f, _ in rows]
    peak = max(ref)
    return max(abs(b - r) for (_, b), r in zip(rows, ref)) / peak, ref


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--radius", type=float, required=True, help="球の幾何半径 [m]")
    ap.add_argument("--fit", type=float, default=0.05,
                    help="実効半径を探す範囲 (幾何半径に対する比, 既定 ±5%%)")
    ap.add_argument("--tol", type=float, default=0.06,
                    help="帯域内ピークに対する最大誤差の許容 (既定 6%%)")
    ap.add_argument("log")
    args = ap.parse_args()

    rows = read_cross_section(args.log)
    if not rows:
        print(f"NG: {args.log} に cross section 表がない (planewave 入力か?)")
        return 1
    if any(not (b > 0 and math.isfinite(b)) for _, b in rows):
        print("NG: backward cross section に正の有限値でないものがある")
        return 1

    n = 200
    cands = [args.radius * (1 - args.fit + 2 * args.fit * i / n) for i in range(n + 1)]
    err, a_eff = min((max_err(rows, a)[0], a) for a in cands)
    _, ref = max_err(rows, a_eff)
    ref0 = [pec_backscatter(args.radius, f) for f, _ in rows]

    print(f"  幾何半径 {args.radius:g} m, 実効半径 {a_eff:.5g} m "
          f"({(a_eff / args.radius - 1) * 100:+.1f}%)")
    print("  frequency[Hz]   ka   FDTD[m^2]  Mie(a)[m^2]  Mie(a')[m^2]")
    for (f, b), r0, r in zip(rows, ref0, ref):
        ka = 2 * math.pi * f / C * args.radius
        print(f"  {f:13.5e}  {ka:4.2f}  {b:10.4e}  {r0:10.4e}  {r:10.4e}")
    if err > args.tol:
        print(f"NG: 最大誤差 {err * 100:.1f}% (ピーク比) > 許容 {args.tol * 100:.1f}%")
        return 1
    print(f"OK: 最大誤差 {err * 100:.1f}% (ピーク比, 許容 {args.tol * 100:.1f}%)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
