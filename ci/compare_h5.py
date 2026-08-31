#!/usr/bin/env python3
"""直列版と MPI 版の time_series_data.h5 を比較する (CI の回帰チェック用)。

並列 HDF5 では H5Gcreate / H5Dcreate が集団操作で、データ空間の次元が rank
ごとに食い違うとデータセットが 0 要素になったり H5Fclose でデッドロックしたり
する。ハングしないことだけでなく、中身が直列版と一致することまで確認する。

E/H/P/P_loss/Surface と時刻・格子などのメタデータは完全一致 (差 0) を要求する。
Eiter/Hiter は MPI の総和順序が変わるぶん丸め誤差が出るので緩い許容差を使う。

    python3 ci/compare_h5.py <serial.h5> <mpi.h5>

一致すれば exit 0、差があれば内容を表示して exit 1。
"""
import sys

import h5py
import numpy as np

# MPI の総和順序依存で丸め誤差が出るデータセット (名前 -> 絶対許容差)
TOLERANT = {"Eiter": 1e-12, "Hiter": 1e-10}

# 直列版のみが出力するデータセット (現在は無し)。
# 片方にしか無いデータセットが増えたら、除外する前に「そもそも両方で出すべき
# ものではないか」を先に考えること。GUI は同じファイルを読むので、ビルド構成で
# データセットが欠けると機能が黙って落ちる。
SERIAL_ONLY = set()


def dataset_names(f):
    names = []
    f.visit(names.append)
    return {n for n in names if isinstance(f[n], h5py.Dataset)}


def maxdiff(a, b):
    if a.size == 0:
        return 0.0
    if a.dtype.names:
        return max(
            float(np.max(np.abs(a[k].astype(float) - b[k].astype(float))))
            for k in a.dtype.names
        )
    if a.dtype.kind in "fiu":
        return float(np.max(np.abs(a.astype(float) - b.astype(float)))) if a.size else 0.0
    return 0.0 if np.array_equal(a, b) else float("inf")


def main(path_ref, path_mpi):
    ref = h5py.File(path_ref, "r")
    mpi = h5py.File(path_mpi, "r")

    kr, km = dataset_names(ref), dataset_names(mpi)
    problems = []

    missing = {k for k in kr - km if k.split("/")[-1] not in SERIAL_ONLY}
    for k in sorted(missing):
        problems.append(f"missing in MPI output: {k}")
    for k in sorted(km - kr):
        problems.append(f"unexpected in MPI output: {k}")

    for k in sorted(kr & km):
        a, b = ref[k], mpi[k]
        if a.shape != b.shape:
            problems.append(f"shape differs: {k} serial={a.shape} mpi={b.shape}")
            continue
        name = k.split("/")[-1]
        d = maxdiff(a[()], b[()])
        if d > TOLERANT.get(name, 0.0):
            problems.append(f"value differs: {k} maxdiff={d:.3e}")

    if problems:
        print(f"FAILED: {path_ref} vs {path_mpi}")
        for p in problems[:40]:
            print("  " + p)
        if len(problems) > 40:
            print(f"  ... and {len(problems) - 40} more")
        return 1

    print(f"OK: {path_mpi} matches {path_ref} "
          f"({len(kr & km)} datasets)")
    return 0


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    sys.exit(main(sys.argv[1], sys.argv[2]))
