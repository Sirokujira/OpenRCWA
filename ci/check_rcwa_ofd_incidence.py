#!/usr/bin/env python3
""".ofd の RCWA の入射角キー rcwaincidence の回帰チェック。

    rcwaincidence = <theta[deg]> <phi[deg]>

RCWA の入射角は従来 planewave = θ φ pol の θ φ から取っていた。しかし planewave は
FDTD の入射方向 (波が来る向き。OpenFDTD-X の既定は θ = 90°) と同じ行なので、
GUI で平面波を有効にしたまま RCWA を選ぶと「theta must be within (-90, 90)」で
止まる。rcwaincidence は RCWA 専用の角度で、planewave より優先する。

orcwa を実際に走らせて次を確認する:

  1. 解析解: 空気 → ガラス (ε = 2.25) の単一界面で R_TE = |r_s|²、R_TM = |r_p|²
     (Fresnel)。Brewster 角 atan(1.5) で R_TM = 0。一様な界面なので φ によらない。
  2. 互換: rcwaincidence = θ φ と planewave = θ φ 1 の結果がバイト単位で一致。
  3. 優先: rcwaincidence があれば planewave の θ (FDTD の既定 90° でも) は使わない。
  4. エネルギー保存: 無損失の格子にコニカル入射 (φ ≠ 0) で R + T = 1。
  5. 鏡映対称: 左右対称な格子で (θ, φ = 0) と (−θ, 0) と (θ, 180°) が一致。
  6. 入力エラー: |θ| >= 90、値の数の誤り、rcwalayer なしの rcwaincidence。

    python3 ci/check_rcwa_ofd_incidence.py --orcwa bin/orcwa

すべて満たせば exit 0、外れれば内容を表示して exit 1。
"""
import argparse
import csv
import math
import os
import subprocess
import sys
import tempfile

PERIOD = 0.4e-6
GRATING = ["1 1 0.5 0", "4.0 1.0 0.4 2.5e-7", "2.25 2.25 0.5 0"]
FREQ = "frequency1 = 4.3e14 7.5e14 9"


def write_and_run(orcwa, work, name, lines):
    d = os.path.join(work, name)
    os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, "in.ofd"), "w", encoding="utf-8") as fp:
        fp.write("\n".join(["OpenRCWA 4 2", f"title = {name}"] + lines + ["end"]) + "\n")
    r = subprocess.run([orcwa, "-n", "2", "in.ofd"], cwd=d, stdout=subprocess.PIPE,
                       stderr=subprocess.STDOUT, encoding="utf-8", errors="replace")
    return d, r


def read_orders(work, name):
    with open(os.path.join(work, name, "rcwa_orders.csv"), encoding="utf-8") as fp:
        return {(row["frequency[Hz]"], row["pol"], int(row["m"]), int(row["n"])):
                (float(row["R"]), float(row["T"])) for row in csv.DictReader(fp)}


def solve(orcwa, work, name, rcwa, layers, extra):
    lines = [f"rcwa = {rcwa}"] + [f"rcwalayer = {l}" for l in layers] + extra + [FREQ]
    d, r = write_and_run(orcwa, work, name, lines)
    if r.returncode != 0:
        print(r.stdout[-2000:])
        raise RuntimeError(f"{name}: orcwa exited {r.returncode}")
    path = os.path.join(d, "rcwa_efficiency.csv")
    with open(path, encoding="utf-8") as fp:
        rows = [{k: float(v) for k, v in row.items()} for row in csv.DictReader(fp)]
    with open(path, "rb") as fp:
        raw = fp.read()
    return rows, raw, r.stdout


def fresnel(theta_deg, n1=1.0, n2=1.5):
    ti = math.radians(theta_deg)
    st = n1 * math.sin(ti) / n2
    ci, ct = math.cos(ti), math.sqrt(1 - st * st)
    rs = (n1 * ci - n2 * ct) / (n1 * ci + n2 * ct)
    rp = (n2 * ci - n1 * ct) / (n2 * ci + n1 * ct)
    return rs * rs, rp * rp


def check(cond, msg, failures):
    print(("  ok  " if cond else "  NG  ") + msg)
    if not cond:
        failures.append(msg)


def max_diff(a, b):
    keys = ("R_TE", "T_TE", "R_TM", "T_TM")
    return max(abs(x[k] - y[k]) for x, y in zip(a, b) for k in keys)


def main():
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    ap = argparse.ArgumentParser()
    ap.add_argument("--orcwa", required=True)
    args = ap.parse_args()
    orcwa = os.path.abspath(args.orcwa)
    failures = []

    with tempfile.TemporaryDirectory() as work:
        print("[1] 単一界面 (空気 → ガラス) の Fresnel 係数と Brewster 角")
        interface = ["1 1 0.5 0", "2.25 2.25 0.5 0"]
        brewster = math.degrees(math.atan(1.5))
        for th, ph in ((0, 0), (30, 0), (brewster, 0), (70, 0), (45, 37)):
            rows, _, log = solve(orcwa, work, f"fresnel_{th:.2f}_{ph}", f"2 {PERIOD}",
                                 interface, [f"rcwaincidence = {th!r} {ph}"])
            rs, rp = fresnel(th)
            worst = max(max(abs(r["R_TE"] - rs), abs(r["R_TM"] - rp)) for r in rows)
            cons = max(max(abs(r["R_TE"] + r["T_TE"] - 1), abs(r["R_TM"] + r["T_TM"] - 1))
                       for r in rows)
            check(worst < 1e-6 and cons < 1e-6,
                  f"θ={th:.2f} φ={ph}: max|R - Fresnel| = {worst:.1e}, max|R+T-1| = {cons:.1e}"
                  f" (R_TE={rs:.5f}, R_TM={rp:.5f})", failures)
            if abs(th - brewster) < 1e-9:
                check(max(r["R_TM"] for r in rows) < 1e-9,
                      f"Brewster 角 {th:.3f}°: R_TM = {max(r['R_TM'] for r in rows):.1e}",
                      failures)
            if th == 30:
                check("theta=30" in log, "ログに入射角 (theta=30) が出る", failures)

        print("[2] rcwaincidence = θ φ と planewave = θ φ 1 が一致")
        a, araw, _ = solve(orcwa, work, "inc", f"8 {PERIOD}", GRATING,
                           ["rcwaincidence = 30 20"])
        b, braw, _ = solve(orcwa, work, "pw", f"8 {PERIOD}", GRATING,
                           ["planewave = 30 20 1"])
        check(araw == braw, "rcwa_efficiency.csv がバイト単位で一致", failures)

        print("[3] rcwaincidence は planewave より優先 (FDTD の既定 θ = 90° でも止まらない)")
        c, craw, _ = solve(orcwa, work, "both", f"8 {PERIOD}", GRATING,
                           ["planewave = 90 0 1", "rcwaincidence = 30 20"])
        check(craw == araw, "planewave = 90 0 1 があっても rcwaincidence の角度で解く", failures)
        _, r = write_and_run(orcwa, work, "pw90",
                             [f"rcwa = 8 {PERIOD}"] + [f"rcwalayer = {l}" for l in GRATING] +
                             ["planewave = 90 0 1", FREQ])
        check(r.returncode != 0 and "within (-90, 90)" in r.stdout,
              "rcwaincidence が無ければ従来どおり planewave の θ を使う (90° は入力エラー)",
              failures)

        print("[4] 無損失の格子にコニカル入射: R + T = 1")
        for th, ph in ((40, 30), (60, 75)):
            rows, _, _ = solve(orcwa, work, f"conical_{th}_{ph}", f"8 {PERIOD}", GRATING,
                               [f"rcwaincidence = {th} {ph}"])
            cons = max(max(abs(r["R_TE"] + r["T_TE"] - 1), abs(r["R_TM"] + r["T_TM"] - 1))
                       for r in rows)
            check(cons < 1e-6, f"θ={th} φ={ph}: max|R+T-1| = {cons:.1e} (< 1e-6)", failures)

        print("[5] 左右対称な格子の鏡映: (θ, 0) = (-θ, 0) = (θ, 180)")
        p, _, _ = solve(orcwa, work, "m_p", f"8 {PERIOD}", GRATING, ["rcwaincidence = 25 0"])
        m, _, _ = solve(orcwa, work, "m_m", f"8 {PERIOD}", GRATING, ["rcwaincidence = -25 0"])
        f, _, _ = solve(orcwa, work, "m_f", f"8 {PERIOD}", GRATING, ["rcwaincidence = 25 180"])
        check(max_diff(p, m) < 1e-6, f"(25, 0) と (-25, 0): {max_diff(p, m):.1e}", failures)
        check(max_diff(p, f) < 1e-6, f"(25, 0) と (25, 180): {max_diff(p, f):.1e}", failures)
        # 次数ごとでは m が反転する: (25, 0) の m 次 = (-25, 0) の -m 次
        op, om = read_orders(work, "m_p"), read_orders(work, "m_m")
        worst = max((abs(v[0] - om[(f_, pol, -mm, nn)][0]) + abs(v[1] - om[(f_, pol, -mm, nn)][1])
                     if (f_, pol, -mm, nn) in om else float("inf"))
                    for (f_, pol, mm, nn), v in op.items())
        refl = sorted({mm for (_, _, mm, _), v in op.items() if v[0] > 0 and mm != 0})
        check(len(op) == len(om) and worst < 1e-6 and refl == [-1],
              f"次数ごとに m -> -m で一致 ({len(op)} 行、{worst:.1e})、"
              f"θ=+25 の反射 1 次は m = {refl} (期待 [-1])", failures)
        # 対称性が自明に成り立っているだけでないこと (斜入射で結果が変わること)
        n0, _, _ = solve(orcwa, work, "m_0", f"8 {PERIOD}", GRATING, [])
        check(max_diff(p, n0) > 1e-3, f"斜入射は法線入射と異なる: {max_diff(p, n0):.2e}", failures)

        print("[6] 入力エラー")
        base = [f"rcwa = 8 {PERIOD}"] + [f"rcwalayer = {l}" for l in GRATING] + [FREQ]
        for name, lines, msg in (
                ("theta90", base + ["rcwaincidence = 90 0"], "rcwaincidence theta must be within"),
                ("thetam90", base + ["rcwaincidence = -95 0"], "rcwaincidence theta must be within"),
                ("few", base + ["rcwaincidence = 30"], "invalid rcwaincidence data"),
                ("many", base + ["rcwaincidence = 30 0 1"], "invalid rcwaincidence data"),
                ("nolayer", ["rcwaincidence = 30 0", FREQ], "rcwaincidence needs rcwalayer")):
            _, r = write_and_run(orcwa, work, "err_" + name, lines)
            check(r.returncode != 0 and msg in r.stdout,
                  f"{name}: exit {r.returncode}, '{msg}'", failures)

    if failures:
        print(f"NG: {len(failures)} 件")
        return 1
    print("OK: .ofd の RCWA 入射角 (rcwaincidence)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
