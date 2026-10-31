#!/usr/bin/env python3
""".ofd の 2D 格子 (rcwa = Nx Λx Ny Λy / rcwalayer の fill_y) の回帰チェック。

orcwa を実際に走らせ、厳密に成り立つはずの関係で判定する:

  0. fill の違う層の混在: 一様層 (eps1 = eps2) の fill は結果に効かない。
     以前は層ごとの fill が違うだけで orcwa が落ちていた。
  1. 1D との一致: fill_y = 1 (y 一様) の 2D 入力は、同じ構造の 1D 入力と
     同じ R/T になる (法線入射・斜め入射の両方)。
  2. C4 対称: 正方格子の正方形柱を法線入射で解くと TE と TM が一致し、
     無損失なので R+T = 1。
  3. 90° 回転: y 方向に伸びた線格子 (fill=f, fill_y=1) と x 方向に伸びた
     線格子 (fill=1, fill_y=f) は、法線入射で TE と TM が入れ替わる。
  4. .orcwa 経路との一致 (--orcwa-rcwa を渡したとき): 同じ矩形柱を
     .orcwa (geometry の箱) で書き、orcwa_rcwa で解いた結果と一致する。

    python3 ci/check_rcwa_ofd2d.py --orcwa bin/orcwa [--orcwa-rcwa bin/orcwa_rcwa]

  5. 入力エラー: rcwa のトークン数の誤り、1D での fill_y < 1、fill_y > 1、
     値が多すぎる rcwalayer、Λy <= 0 は終了コード != 0 で止まる。

すべて満たせば exit 0、外れれば内容を表示して exit 1。
"""
import argparse
import csv
import os
import subprocess
import sys
import tempfile

C0 = 2.99792458e8
PERIOD = 0.5e-6           # 周期 [m] (x, y とも)
THICK = 0.2e-6            # 格子層の厚さ [m]
EPS_PILLAR = 2.25         # 柱 (ガラス)
EPS_SUB = 2.25            # 透過側基板
# λ = 0.6-0.9 μm (Wood アノマリーの λ = Λ = 0.5 μm と λ = nΛ = 0.75 μm を避ける)
# frequency1 は等間隔掃引 (f_min f_max 分割数)。λ = 0.9 → 0.62 μm を 5 点
F_MIN, F_MAX, NDIV = C0 / 0.9e-6, C0 / 0.62e-6, 4
FREQ_LINE = f"frequency1 = {F_MIN:.10e} {F_MAX:.10e} {NDIV}"
NH = 4                    # 高調波次数 (x, y とも)
# 厳密に一致するはずの量の許容差。LAPACKE 経路では CSV の桁 (1e-9) まで一致するが、
# Eigen フォールバック経路 (Windows) は縮退対策の非一様な対角摂動 (~1e-11·‖A‖) で
# 対称性がわずかに破れ、実測 3.5e-8 ずれる。物理的な差 (1e-2 程度) よりは十分小さい。
TOL = 1e-6


def run(binary, workdir, name, text, args=()):
    d = os.path.join(workdir, name)
    os.makedirs(d, exist_ok=True)
    path = os.path.join(d, name + os.path.splitext(text[0])[1])
    with open(path, "w", encoding="utf-8") as fp:
        fp.write(text[1])
    r = subprocess.run([binary, *args, os.path.basename(path)], cwd=d,
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                       encoding="utf-8", errors="replace")
    if r.returncode != 0:
        print(r.stdout[-2000:])
        raise RuntimeError(f"{name}: {os.path.basename(binary)} exited {r.returncode}")
    return d


def ofd(name, rcwa, layers, theta=0.0, phi=0.0):
    """.ofd を組み立てる。layers は rcwalayer 行の値の文字列のリスト。"""
    lines = ["OpenRCWA 4 2", f"title = {name}", f"rcwa = {rcwa}"]
    lines += [f"rcwalayer = {l}" for l in layers]
    lines.append(f"planewave = {theta} {phi} 1")
    lines.append(FREQ_LINE)
    lines.append("end")
    return (".ofd", "\n".join(lines) + "\n")


def read_ofd_csv(d):
    rows = []
    with open(os.path.join(d, "rcwa_efficiency.csv"), encoding="utf-8") as fp:
        for r in csv.DictReader(fp):
            rows.append({k: float(v) for k, v in r.items()})
    return rows


def layers_pillar(fx, fy):
    sub = f"{EPS_SUB} {EPS_SUB} 0.5 0"
    grat = f"{EPS_PILLAR} 1 {fx} {THICK} 0 0 {fy}"
    return ["1 1 0.5 0", grat, sub]


def check(cond, msg, failures):
    print(("  ok  " if cond else "  NG  ") + msg)
    if not cond:
        failures.append(msg)


def maxdiff(a, b, keys):
    return max(abs(x[k] - y[k]) for x, y in zip(a, b) for k in keys)


def main():
    # Windows の CI はパイプ出力が cp1252 になり、日本語の print で落ちるため
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    ap = argparse.ArgumentParser()
    ap.add_argument("--orcwa", required=True)
    ap.add_argument("--orcwa-rcwa", default=None)
    args = ap.parse_args()
    orcwa = os.path.abspath(args.orcwa)
    failures = []
    keys = ("R_TE", "T_TE", "R_TM", "T_TM")
    two_d = f"{NH} {PERIOD} {NH} {PERIOD}"
    one_d = f"{NH} {PERIOD}"

    with tempfile.TemporaryDirectory() as work:
        print("[0] 層ごとに fill が違っても解ける (一様層の fill は結果に効かない)")
        # 以前は fill の違う層が混ざると層が捨てられて orcwa が落ちていた
        for rc, tag in ((one_d, "1D"), (two_d, "2D")):
            fy = "" if tag == "1D" else " 0 0 0.7"
            same = ["1 1 0.3 0", f"{EPS_PILLAR} 1 0.3 {THICK}{fy}", f"{EPS_SUB} {EPS_SUB} 0.3 0"]
            mixed = ["1 1 0.5 0", f"{EPS_PILLAR} 1 0.3 {THICK}{fy}", f"{EPS_SUB} {EPS_SUB} 0.8 0"]
            if tag == "2D":
                same = [same[0] + " 0 0 0.7", same[1], same[2] + " 0 0 0.7"]
                mixed = [mixed[0] + " 0 0 0.2", mixed[1], mixed[2]]
            a = read_ofd_csv(run(orcwa, work, f"same_{tag}", ofd("same", rc, same, 20, 30), ["-n", "2"]))
            b = read_ofd_csv(run(orcwa, work, f"mixed_{tag}", ofd("mixed", rc, mixed, 20, 30), ["-n", "2"]))
            d = maxdiff(a, b, keys)
            check(d < 1e-7, f"{tag}: max|diff| = {d:.2e} (< 1e-7)", failures)

        print("[1] fill_y = 1 の 2D 入力 = 1D 入力")
        for theta, phi in ((0, 0), (20, 30)):
            l1 = ["1 1 0.5 0", f"{EPS_PILLAR} 1 0.4 {THICK}", f"{EPS_SUB} {EPS_SUB} 0.5 0"]
            l2 = ["1 1 0.5 0", f"{EPS_PILLAR} 1 0.4 {THICK} 0 0 1", f"{EPS_SUB} {EPS_SUB} 0.5 0"]
            a = read_ofd_csv(run(orcwa, work, f"oned_{theta}", ofd("1d", one_d, l1, theta, phi), ["-n", "2"]))
            b = read_ofd_csv(run(orcwa, work, f"twod_{theta}", ofd("2d", two_d, l2, theta, phi), ["-n", "2"]))
            d = maxdiff(a, b, keys)
            check(d < TOL, f"theta={theta} phi={phi}: max|diff| = {d:.2e} (< {TOL:g})", failures)

        print("[2] 正方形柱 (C4 対称) の法線入射: TE = TM, R+T = 1")
        sq = read_ofd_csv(run(orcwa, work, "square", ofd("sq", two_d, layers_pillar(0.6, 0.6)), ["-n", "2"]))
        d = max(max(abs(r["R_TE"] - r["R_TM"]), abs(r["T_TE"] - r["T_TM"])) for r in sq)
        check(d < TOL, f"max|TE - TM| = {d:.2e} (< {TOL:g})", failures)
        e = max(max(abs(r["R_TE"] + r["T_TE"] - 1), abs(r["R_TM"] + r["T_TM"] - 1)) for r in sq)
        check(e < TOL, f"max|R+T-1| = {e:.2e} (< {TOL:g})", failures)
        nontrivial = max(abs(r["R_TE"] - sq[0]["R_TE"]) for r in sq)
        check(nontrivial > 1e-4, f"R が波長で変化している (幅 {nontrivial:.2e})", failures)

        print("[3] 90° 回転: y 方向の線格子と x 方向の線格子で TE/TM が入れ替わる")
        ly = read_ofd_csv(run(orcwa, work, "lines_y", ofd("ly", two_d, layers_pillar(0.4, 1)), ["-n", "2"]))
        lx = read_ofd_csv(run(orcwa, work, "lines_x", ofd("lx", two_d, layers_pillar(1, 0.4)), ["-n", "2"]))
        d = max(max(abs(a["R_TE"] - b["R_TM"]), abs(a["T_TE"] - b["T_TM"]),
                    abs(a["R_TM"] - b["R_TE"]), abs(a["T_TM"] - b["T_TE"])) for a, b in zip(ly, lx))
        check(d < TOL, f"max|TE(y 線) - TM(x 線)| ほか = {d:.2e} (< {TOL:g})", failures)
        split = max(abs(r["R_TE"] - r["R_TM"]) for r in ly)
        check(split > 1e-3, f"線格子は偏波依存がある (|R_TE-R_TM| 最大 {split:.2e})", failures)

        if args.orcwa_rcwa:
            print("[4] .orcwa 経路 (orcwa_rcwa) との一致: 0.3x0.2 μm の矩形柱")
            orcwa_rcwa = os.path.abspath(args.orcwa_rcwa)
            fx, fy = 0.6, 0.4
            ref = read_ofd_csv(run(orcwa, work, "rect", ofd("rect", two_d, layers_pillar(fx, fy)), ["-n", "2"]))
            for pol, col in ((1, "TM"), (2, "TE")):
                # .orcwa: 20x20 セルで柱の辺 (0.3, 0.2 μm) がセル境界 (0.025 μm 刻み) に乗る
                text = "\n".join([
                    "OpenRCWA 4 2", "title = rect pillar",
                    f"xmesh = 0 20 {PERIOD}", f"ymesh = 0 20 {PERIOD}",
                    f"zmesh = -4e-07 10 {-THICK} 10 0.0 10 4e-07",
                    f"material = 1 {EPS_PILLAR} 0 1 0",
                    f"geometry = 2 1 0 {PERIOD} 0 {PERIOD} -4e-07 {-THICK}",
                    f"geometry = 2 1 0 {fx * PERIOD} 0 {fy * PERIOD} {-THICK} 0",
                    f"planewave = 0 0 {pol}", "pbc = 1 1 0", f"rcwaorder = {NH} {NH}",
                    FREQ_LINE, "end"]) + "\n"
                d = run(orcwa_rcwa, work, f"orcwa_{col}", (".orcwa", text), ["-o", "out.csv"])
                # 書式: "# lambda[um], R, T, A, R+T"
                with open(os.path.join(d, "out.csv"), encoding="utf-8") as fp:
                    data = [[float(v) for v in line.split(",")]
                            for line in fp if line.strip() and not line.startswith("#")]
                worst, matched = 0.0, 0
                for lam_um, r, t, *_ in data:
                    o = min(ref, key=lambda x: abs(x["lambda[m]"] - lam_um * 1e-6))
                    if abs(o["lambda[m]"] / (lam_um * 1e-6) - 1) > 1e-6:
                        continue
                    matched += 1
                    worst = max(worst, abs(r - o["R_" + col]), abs(t - o["T_" + col]))
                check(matched == len(ref) and worst < 1e-6,
                      f"{col}: max|diff| = {worst:.2e} ({matched}/{len(ref)} 点が対応, < 1e-6)",
                      failures)

        print("[5] 入力エラーは終了コード != 0 とメッセージで止まる")
        bad = (
            ("4 5e-07 4", "2.25 1 0.5 2e-07", "invalid rcwa data"),
            (one_d, "2.25 1 0.5 2e-07 0 0 0.5", "fill_y < 1 needs a 2D grating"),
            (two_d, "2.25 1 0.5 2e-07 0 0 1.5", "invalid rcwalayer data"),
            (two_d, "2.25 1 0.5 2e-07 0 0 0.5 9", "too many values"),
            (f"{NH} {PERIOD} {NH} 0", "2.25 1 0.5 2e-07", "invalid rcwa data"),
        )
        for i, (rc, layer, msg) in enumerate(bad):
            text = ofd("bad", rc, ["1 1 0.5 0", layer, f"{EPS_SUB} {EPS_SUB} 0.5 0"])
            d = os.path.join(work, f"bad{i}")
            os.makedirs(d, exist_ok=True)
            with open(os.path.join(d, "b.ofd"), "w", encoding="utf-8") as fp:
                fp.write(text[1])
            r = subprocess.run([orcwa, "-n", "1", "b.ofd"], cwd=d, stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT, encoding="utf-8", errors="replace")
            check(r.returncode != 0 and msg in r.stdout,
                  f"rcwa = {rc} / rcwalayer = {layer} -> rc={r.returncode}, '{msg}'", failures)

    if failures:
        print(f"NG: {len(failures)} 件")
        return 1
    print("OK: .ofd の 2D 格子")
    return 0


if __name__ == "__main__":
    sys.exit(main())
