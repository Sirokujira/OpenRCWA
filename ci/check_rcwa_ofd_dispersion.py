#!/usr/bin/env python3
""".ofd の RCWA 分散材料 (rcwadispersion) の回帰チェック。

    rcwadispersion = <layer(1..)> <side(1|2)> <einf> <ae> <be> <ce>   [rad/s]
    eps(ω) = einf + Σ ae² / (ce² − ω² − i·be·ω)

orcwa を実際に走らせて次を確認する:

  1. 解析解: 空気 / Drude 金属薄膜 (30 nm) / ガラスの R, T を薄膜の特性行列法
     (Airy の式) と比べる。法線入射と θ=30° の TE / TM。
  2. .orcwa 経路との一致 (--orcwa-rcwa を渡したとき): 同じ構造の
     data/sample/rcwa_metal_drude.orcwa を orcwa_rcwa で解いた TM の R, T, A。
  3. 極の加算: 1 つの極を ae² の和が同じ 2 行に分けても結果は同じ。
  4. eps1 / eps2 の割り当て: 分散材料を side 1 (fill f) に置いた格子と、
     side 2 (fill 1−f) に置いた格子は平行移動しただけなので一致する。
  5. 衝突なし (be = 0) の Drude 金属は無損失なので R+T = 1。
  6. 入力エラー (層番号の範囲外、side 3、be < 0、値の数) は終了コード != 0。

    python3 ci/check_rcwa_ofd_dispersion.py --orcwa bin/orcwa \\
        [--orcwa-rcwa bin/orcwa_rcwa --sample data/sample/rcwa_metal_drude.orcwa]

すべて満たせば exit 0、外れれば内容を表示して exit 1。
"""
import argparse
import cmath
import csv
import math
import os
import subprocess
import sys
import tempfile

C0 = 2.99792458e8
# rcwa_metal_drude.orcwa と同じ Drude 金属: einf=1, ωp=1.36e16, γ=3e13 [rad/s]
WP, GAMMA = 1.36e16, 3.0e13
FILM = 3.0e-8                     # 膜厚 [m]
EPS_SUB = 2.25
F_MIN, F_MAX, NDIV = 4.0e14, 7.0e14, 6
FREQ_LINE = f"frequency1 = {F_MIN:.10e} {F_MAX:.10e} {NDIV}"


def drude(f, wp=WP, gamma=GAMMA):
    w = 2 * math.pi * f
    return 1.0 + wp * wp / complex(-w * w, -gamma * w)


def film_rt(eps_film, f, theta_deg, pol):
    """空気 / 膜 / 基板の R, T (特性行列法)。

    exp(-iωt) の誘電率 (損失 = 正の虚部) を、Macleod の exp(+iωt) 規約
    N = n − ik に直して計算する (|r|², T は規約によらない)。"""
    lam = C0 / f
    n0 = 1.0
    s = n0 * math.sin(math.radians(theta_deg))

    def eta(N):
        c = cmath.sqrt(1 - (s / N) ** 2)
        if (N * c).imag > 0:          # exp(+iωt) で減衰する枝 (Im(N cosθ) <= 0)
            c = -c
        return (N * c if pol == "TE" else N / c), N * c

    N1 = complex(cmath.sqrt(eps_film)).conjugate()
    Ns = complex(math.sqrt(EPS_SUB))
    e0, _ = eta(complex(n0))
    e1, nc1 = eta(N1)
    es, _ = eta(Ns)
    d = 2 * math.pi * nc1 * FILM / lam
    B = cmath.cos(d) + 1j * cmath.sin(d) / e1 * es
    C = 1j * e1 * cmath.sin(d) + cmath.cos(d) * es
    R = abs((e0 * B - C) / (e0 * B + C)) ** 2
    T = 4 * e0.real * es.real / abs(e0 * B + C) ** 2
    return R, T


def run(binary, workdir, name, ext, text, args=()):
    d = os.path.join(workdir, name)
    os.makedirs(d, exist_ok=True)
    path = os.path.join(d, name + ext)
    with open(path, "w", encoding="utf-8") as fp:
        fp.write(text)
    r = subprocess.run([binary, *args, os.path.basename(path)], cwd=d,
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                       encoding="utf-8", errors="replace")
    return d, r


def ofd(rcwa, layers, disp, theta=0.0):
    lines = ["OpenRCWA 4 2", "title = dispersion check", f"rcwa = {rcwa}"]
    lines += [f"rcwalayer = {l}" for l in layers]
    lines += [f"rcwadispersion = {p}" for p in disp]
    lines += [f"planewave = {theta} 0 1", FREQ_LINE, "end"]
    return "\n".join(lines) + "\n"


def solve_ofd(orcwa, work, name, text):
    d, r = run(orcwa, work, name, ".ofd", text, ["-n", "2"])
    if r.returncode != 0:
        print(r.stdout[-2000:])
        raise RuntimeError(f"{name}: orcwa exited {r.returncode}")
    with open(os.path.join(d, "rcwa_efficiency.csv"), encoding="utf-8") as fp:
        return [{k: float(v) for k, v in row.items()} for row in csv.DictReader(fp)], r.stdout


def check(cond, msg, failures):
    print(("  ok  " if cond else "  NG  ") + msg)
    if not cond:
        failures.append(msg)


FILM_LAYERS = ["1 1 0.5 0", f"1 1 0.5 {FILM}", f"{EPS_SUB} {EPS_SUB} 0.5 0"]
DRUDE_POLE = f"1 {WP:.6e} {GAMMA:.6e} 0"   # einf ae be ce


def main():
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    ap = argparse.ArgumentParser()
    ap.add_argument("--orcwa", required=True)
    ap.add_argument("--orcwa-rcwa", default=None)
    ap.add_argument("--sample", default=None, help="rcwa_metal_drude.orcwa のパス")
    args = ap.parse_args()
    orcwa = os.path.abspath(args.orcwa)
    failures = []
    keys = ("R_TE", "T_TE", "R_TM", "T_TM")

    with tempfile.TemporaryDirectory() as work:
        print("[1] Drude 金属薄膜の R, T を特性行列法 (解析解) と比べる")
        film = None
        for theta in (0.0, 30.0):
            rows, log = solve_ofd(orcwa, work, f"film{int(theta)}",
                                  ofd("2 5e-07", FILM_LAYERS,
                                      [f"2 1 {DRUDE_POLE}", f"2 2 {DRUDE_POLE}"], theta))
            if theta == 0.0:
                film = rows
                check("dispersive (Lorentz) materials : 2" in log,
                      "ログに分散材料の数が出る", failures)
            worst = 0.0
            for r in rows:
                eps = drude(r["frequency[Hz]"])
                for pol in ("TE", "TM"):
                    R, T = film_rt(eps, r["frequency[Hz]"], theta, pol)
                    worst = max(worst, abs(r["R_" + pol] - R), abs(r["T_" + pol] - T))
            check(worst < 1e-6, f"theta={theta:g}: max|RCWA - 解析解| = {worst:.2e} (< 1e-6)",
                  failures)
        a = [1 - r["R_TE"] - r["T_TE"] for r in film]
        check(min(a) > 1e-3, f"損失があるので吸収 A > 0 (最小 {min(a):.3e})", failures)
        rs = [r["R_TE"] for r in film]
        check(max(rs) - min(rs) > 1e-2, f"R が周波数で変わる (幅 {max(rs) - min(rs):.3e})",
              failures)

        if args.orcwa_rcwa and args.sample:
            print("[2] .orcwa 経路 (rcwa_metal_drude.orcwa) との一致")
            with open(args.sample, encoding="utf-8") as fp:
                text = fp.read()
            d, r = run(os.path.abspath(args.orcwa_rcwa), work, "orcwa_drude", ".orcwa",
                       text, ["-o", "out.csv"])
            check(r.returncode == 0, f"orcwa_rcwa rc={r.returncode}", failures)
            with open(os.path.join(d, "out.csv"), encoding="utf-8") as fp:
                data = [[float(v) for v in line.split(",")]
                        for line in fp if line.strip() and not line.startswith("#")]
            worst, matched = 0.0, 0
            for lam_um, R, T, A, *_ in data:
                o = min(film, key=lambda x: abs(x["lambda[m]"] - lam_um * 1e-6))
                if abs(o["lambda[m]"] / (lam_um * 1e-6) - 1) > 1e-6:
                    continue
                matched += 1
                worst = max(worst, abs(R - o["R_TM"]), abs(T - o["T_TM"]),
                            abs(A - (1 - o["R_TM"] - o["T_TM"])))
            check(matched == len(film) and worst < 1e-6,
                  f"TM: max|diff| = {worst:.2e} ({matched}/{len(film)} 点, < 1e-6)", failures)

        print("[3] 極の加算: 1 つの極を ae² の和が同じ 2 行に分けても同じ")
        ae1, ae2 = WP * math.sqrt(0.3), WP * math.sqrt(0.7)
        split = [f"2 1 1 {ae1:.10e} {GAMMA:.6e} 0", f"2 1 1 {ae2:.10e} {GAMMA:.6e} 0",
                 f"2 2 {DRUDE_POLE}"]
        rows, _ = solve_ofd(orcwa, work, "split", ofd("2 5e-07", FILM_LAYERS, split))
        d = max(abs(x[k] - y[k]) for x, y in zip(rows, film) for k in keys)
        check(d < 1e-8, f"max|diff| = {d:.2e} (< 1e-8)", failures)

        print("[4] eps1 / eps2 の割り当て: side 1 (fill f) と side 2 (fill 1-f) が一致")
        lor = "2.0 3.0e15 2.0e14 4.0e15"          # 損失のある Lorentz 極
        g1 = ["1 1 0.5 0", "1 2.25 0.3 1e-07", f"{EPS_SUB} {EPS_SUB} 0.5 0"]
        g2 = ["1 1 0.5 0", "2.25 1 0.7 1e-07", f"{EPS_SUB} {EPS_SUB} 0.5 0"]
        a1, _ = solve_ofd(orcwa, work, "side1", ofd("6 5e-07", g1, [f"2 1 {lor}"], 20))
        a2, _ = solve_ofd(orcwa, work, "side2", ofd("6 5e-07", g2, [f"2 2 {lor}"], 20))
        d = max(abs(x[k] - y[k]) for x, y in zip(a1, a2) for k in keys)
        check(d < 1e-6, f"max|diff| = {d:.2e} (< 1e-6)", failures)
        a = min(1 - r["R_TM"] - r["T_TM"] for r in a1)
        check(a > 1e-4, f"Lorentz 極の損失で A > 0 (最小 {a:.3e})", failures)

        print("[5] 衝突なし (be = 0) の Drude 金属は無損失: R+T = 1")
        lossless = [f"2 1 1 {WP:.6e} 0 0", f"2 2 1 {WP:.6e} 0 0"]
        rows, log = solve_ofd(orcwa, work, "lossless", ofd("2 5e-07", FILM_LAYERS, lossless))
        e = max(max(abs(r["R_TE"] + r["T_TE"] - 1), abs(r["R_TM"] + r["T_TM"] - 1)) for r in rows)
        check(e < 1e-6, f"max|R+T-1| = {e:.2e} (< 1e-6)", failures)
        check("energy not conserved" not in log, "エネルギー保存の警告が出ない", failures)

        print("[6] 入力エラーは終了コード != 0 で止まる")
        bad = (
            (["4 1 1 1e16 1e13 0"], "invalid rcwadispersion data"),
            (["2 3 1 1e16 1e13 0"], "invalid rcwadispersion data"),
            (["2 1 1 1e16 -1e13 0"], "negative be"),
            (["2 1 1 1e16 1e13"], "rcwadispersion"),
            (["2 1 1 1e16 1e13 -1"], "invalid rcwadispersion data"),
        )
        for i, (disp, msg) in enumerate(bad):
            _, r = run(orcwa, work, f"bad{i}", ".ofd", ofd("2 5e-07", FILM_LAYERS, disp), ["-n", "1"])
            check(r.returncode != 0 and msg in r.stdout,
                  f"rcwadispersion = {disp[0]} -> rc={r.returncode}, '{msg}'", failures)

    if failures:
        print(f"NG: {len(failures)} 件")
        return 1
    print("OK: .ofd の分散材料")
    return 0


if __name__ == "__main__":
    sys.exit(main())
