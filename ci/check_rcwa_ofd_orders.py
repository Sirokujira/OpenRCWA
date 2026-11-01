#!/usr/bin/env python3
""".ofd の RCWA の回折次数ごとの効率 (rcwa_orders.csv / HDF5) の回帰チェック。

rcwa_orders.csv: frequency[Hz],lambda[m],pol,m,n,R,T (伝搬する次数だけ)
HDF5: /rcwa/R_orders, /rcwa/T_orders [npol][nfreq][2Nx+1][2Ny+1], order_m, order_n

orcwa を実際に走らせて次を確認する:

  1. 次数の和 = rcwa_efficiency.csv の R, T (両偏波・全周波数)。
  2. 格子方程式: 出てくる次数は入射側 (R) / 透過側 (T) で伝搬する次数に含まれ、
     伝搬すると予測される ±1 次は実際に出てくる。
  3. 対称性: 法線入射の対称な格子では R(+m) = R(−m)、T(+m) = T(−m)。
  4. m の符号: x 方向に傾けた入射 (θ>0, φ=0) では、反射側で伝搬する 1 次は
     m = −1 (sinθ_m = sinθ + mλ/Λ で |sinθ_m| < 1 になる側) だけ。
  5. n の符号: 2D 格子で y 方向に傾けた入射 (φ=90°) では n = −1 だけ。
  6. HDF5 の R_orders / T_orders / order_m / order_n が CSV と一致する
     (h5py があるとき)。

    python3 ci/check_rcwa_ofd_orders.py --orcwa bin/orcwa

すべて満たせば exit 0、外れれば内容を表示して exit 1。
"""
import argparse
import csv
import math
import os
import subprocess
import sys
import tempfile

C0 = 2.99792458e8
PERIOD = 0.3e-6
EPS_SUB = 2.25
LAYERS = ["1 1 0.5 0", "4.0 1.0 0.5 2.0e-7", f"{EPS_SUB} {EPS_SUB} 0.5 0"]


def run(orcwa, work, name, rcwa, layers, theta=0.0, phi=0.0, freq_line=None):
    d = os.path.join(work, name)
    os.makedirs(d, exist_ok=True)
    lines = ["OpenRCWA 4 2", f"title = {name}", f"rcwa = {rcwa}"]
    lines += [f"rcwalayer = {l}" for l in layers]
    lines += [f"planewave = {theta} {phi} 1",
              freq_line or "frequency1 = 4.3e14 7.5e14 16", "end"]
    with open(os.path.join(d, "in.ofd"), "w", encoding="utf-8") as fp:
        fp.write("\n".join(lines) + "\n")
    r = subprocess.run([orcwa, "-n", "2", "in.ofd"], cwd=d, stdout=subprocess.PIPE,
                       stderr=subprocess.STDOUT, encoding="utf-8", errors="replace")
    if r.returncode != 0:
        print(r.stdout[-2000:])
        raise RuntimeError(f"{name}: orcwa exited {r.returncode}")
    with open(os.path.join(d, "rcwa_efficiency.csv"), encoding="utf-8") as fp:
        tot = [{k: float(v) for k, v in row.items()} for row in csv.DictReader(fp)]
    with open(os.path.join(d, "rcwa_orders.csv"), encoding="utf-8") as fp:
        ords = []
        for row in csv.DictReader(fp):
            ords.append({"f": float(row["frequency[Hz]"]), "lam": float(row["lambda[m]"]),
                         "pol": row["pol"], "m": int(row["m"]), "n": int(row["n"]),
                         "R": float(row["R"]), "T": float(row["T"])})
    return d, tot, ords


def propagating(lam, m, n, eps, theta, phi, n_inc=1.0, px=PERIOD, py=PERIOD):
    k0 = 2 * math.pi / lam
    s = n_inc * math.sin(math.radians(theta))
    kx = k0 * s * math.cos(math.radians(phi)) + 2 * math.pi * m / px
    ky = k0 * s * math.sin(math.radians(phi)) + 2 * math.pi * n / py
    return kx * kx + ky * ky < k0 * k0 * eps


def check(cond, msg, failures):
    print(("  ok  " if cond else "  NG  ") + msg)
    if not cond:
        failures.append(msg)


def main():
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    ap = argparse.ArgumentParser()
    ap.add_argument("--orcwa", required=True)
    args = ap.parse_args()
    orcwa = os.path.abspath(args.orcwa)
    failures = []

    with tempfile.TemporaryDirectory() as work:
        cases = {}
        for name, rcwa, layers, th, ph in (
                ("normal", f"8 {PERIOD}", LAYERS, 0, 0),
                ("oblique", f"8 {PERIOD}", LAYERS, 30, 0),
                ("conical2d", f"3 {PERIOD} 3 {PERIOD}",
                 ["1 1 0.5 0", "4.0 1.0 0.5 2.0e-7 0 0 0.6", f"{EPS_SUB} {EPS_SUB} 0.5 0"], 30, 90)):
            cases[name] = (th, ph) + run(orcwa, work, name, rcwa, layers, th, ph)

        print("[1] 次数の和 = rcwa_efficiency.csv の合計")
        for name, (th, ph, d, tot, ords) in cases.items():
            worst = 0.0
            for t in tot:
                for pol in ("TE", "TM"):
                    rs = [o for o in ords if o["f"] == t["frequency[Hz]"] and o["pol"] == pol]
                    worst = max(worst, abs(sum(o["R"] for o in rs) - t["R_" + pol]),
                                abs(sum(o["T"] for o in rs) - t["T_" + pol]))
            check(worst < 1e-7, f"{name}: max|Σ次数 - 合計| = {worst:.2e} (< 1e-7)", failures)

        print("[2] 格子方程式: 出てくる次数は伝搬次数で、伝搬する ±1 次は出てくる")
        for name, (th, ph, d, tot, ords) in cases.items():
            bad = [o for o in ords
                   if (o["R"] > 0 and not propagating(o["lam"], o["m"], o["n"], 1.0, th, ph)) or
                      (o["T"] > 0 and not propagating(o["lam"], o["m"], o["n"], EPS_SUB, th, ph))]
            check(not bad, f"{name}: 伝搬しないはずの次数が {len(bad)} 行", failures)
            missing = 0
            for t in tot:
                for pol in ("TE", "TM"):
                    # 1D (Ny = 0) の格子には y 方向の次数は無い
                    firsts = ((1, 0), (-1, 0)) + (((0, 1), (0, -1)) if "2d" in name else ())
                    for (m, n) in firsts:
                        if not propagating(t["lambda[m]"], m, n, EPS_SUB, th, ph):
                            continue
                        if not any(o["f"] == t["frequency[Hz]"] and o["pol"] == pol and
                                   o["m"] == m and o["n"] == n and o["T"] > 0 for o in ords):
                            missing += 1
            check(missing == 0, f"{name}: 伝搬するはずの 1 次の欠落 {missing} 件", failures)
        n_first = sum(1 for o in cases["normal"][4] if abs(o["m"]) == 1)
        check(n_first > 0, f"normal: 1 次が実際に出る ({n_first} 行)", failures)

        print("[3] 法線入射の対称な格子: R(+m) = R(-m), T(+m) = T(-m)")
        ords = cases["normal"][4]
        worst, pairs = 0.0, 0
        for o in ords:
            if o["m"] <= 0:
                continue
            mate = [p for p in ords if p["f"] == o["f"] and p["pol"] == o["pol"] and
                    p["m"] == -o["m"] and p["n"] == o["n"]]
            if not mate:
                worst = float("inf")
                continue
            pairs += 1
            worst = max(worst, abs(o["R"] - mate[0]["R"]), abs(o["T"] - mate[0]["T"]))
        check(pairs > 0 and worst < 1e-7, f"{pairs} 組、max|diff| = {worst:.2e} (< 1e-7)", failures)

        print("[4] m の符号: x 方向に傾けた入射で反射側の 1 次は m = -1 だけ")
        ords = cases["oblique"][4]
        refl1 = {o["m"] for o in ords if o["R"] > 0 and o["m"] != 0}
        check(refl1 == {-1}, f"反射側で出る m≠0 の次数 = {sorted(refl1)} (期待 [-1])", failures)

        print("[5] n の符号: 2D で y 方向に傾けた入射 (phi=90) で反射側の 1 次は n = -1 だけ")
        ords = cases["conical2d"][4]
        refl1 = {(o["m"], o["n"]) for o in ords if o["R"] > 0 and (o["m"], o["n"]) != (0, 0)}
        check(refl1 == {(0, -1)}, f"反射側で出る次数 = {sorted(refl1)} (期待 [(0, -1)])", failures)

        try:
            import h5py
        except ImportError:
            h5py = None
        if h5py is None:
            print("[6] h5py が無いので HDF5 の確認は飛ばす")
        else:
            print("[6] HDF5 の R_orders / T_orders が CSV と一致")
            for name, (th, ph, d, tot, ords) in cases.items():
                with h5py.File(os.path.join(d, "time_series_data.h5"), "r") as f:
                    R, T = f["rcwa/R_orders"][()], f["rcwa/T_orders"][()]
                    om, on = list(f["rcwa/order_m"][()]), list(f["rcwa/order_n"][()])
                    freqs = [row[0] for row in f["rcwa/spectrum"][0]]
                worst, nz = 0.0, 0
                for ip, pol in enumerate(("TE", "TM")):
                    for ifr, fr in enumerate(freqs):
                        for i, m in enumerate(om):
                            for j, n in enumerate(on):
                                rows = [o for o in ords if o["f"] == float("%.8e" % fr) and
                                        o["pol"] == pol and o["m"] == m and o["n"] == n]
                                r = rows[0]["R"] if rows else 0.0
                                t = rows[0]["T"] if rows else 0.0
                                nz += (R[ip, ifr, i, j] != 0) or (T[ip, ifr, i, j] != 0)
                                worst = max(worst, abs(R[ip, ifr, i, j] - r), abs(T[ip, ifr, i, j] - t))
                check(worst < 1e-7 and nz == len(ords),
                      f"{name}: max|HDF5 - CSV| = {worst:.2e}, 非ゼロ {nz} = CSV {len(ords)} 行",
                      failures)

    if failures:
        print(f"NG: {len(failures)} 件")
        return 1
    print("OK: .ofd の回折次数ごとの効率")
    return 0


if __name__ == "__main__":
    sys.exit(main())
