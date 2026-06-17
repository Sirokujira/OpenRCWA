#include "rcwa/RCWADriver.h"
#include "rcwa/RCWASolver.h"
#include "rcwa/Layer.h"

#include <algorithm>
#include <set>
#include <iostream>

using namespace Eigen;

namespace {

// 全スラブで共有する x 方向ブレークポイント (μm) を作る。
// addLayer は全層が同一の境界座標を持つことを要求するため、各直方体の x エッジを
// 集めた共通格子を用いる。
std::vector<scalar> buildCommonXGrid(const RCWAProblem& prob)
{
    std::set<scalar> xs;
    xs.insert(prob.xmin);
    xs.insert(prob.xmax);
    for (const auto& b : prob.boxes) {
        scalar x0 = std::max(b.x0, prob.xmin);
        scalar x1 = std::min(b.x1, prob.xmax);
        if (x0 > prob.xmin && x0 < prob.xmax) xs.insert(x0);
        if (x1 > prob.xmin && x1 < prob.xmax) xs.insert(x1);
    }
    return std::vector<scalar>(xs.begin(), xs.end());
}

// 指定セル中央 (xc) と z スラブ中央 (zc) における誘電率を返す
scalex epsAt(const RCWAProblem& prob, scalar xc, scalar zc)
{
    scalex e = prob.backgroundEps;
    for (const auto& b : prob.boxes) {
        if (zc < b.z0 || zc > b.z1) continue;
        if (xc >= b.x0 && xc <= b.x1) {
            if (b.material >= 0 &&
                b.material < static_cast<int>(prob.materialEps.size())) {
                e = prob.materialEps[b.material];
            }
        }
    }
    return e;
}

// 共通 x 格子・所与の z スラブから 1D 用の Layer を生成する
Layer makeLayer(const RCWAProblem& prob,
                const std::vector<scalar>& xgrid,
                scalar zc)
{
    const int nCellX = static_cast<int>(xgrid.size()) - 1;

    VectorXs coordX(xgrid.size());
    for (size_t i = 0; i < xgrid.size(); ++i) coordX[i] = xgrid[i];

    // 1D 格子: y 方向は単一セル
    VectorXs coordY(2);
    coordY[0] = prob.ymin;
    coordY[1] = prob.ymax;

    // Layer は eps を (nCellY 行 x nCellX 列) で要求する。1D 格子なので y は 1 セル。
    MatrixXcs eps(1, nCellX);
    for (int i = 0; i < nCellX; ++i) {
        scalar xc = 0.5 * (xgrid[i] + xgrid[i + 1]);
        eps(0, i) = epsAt(prob, xc, zc);
    }

    return Layer(coordX, coordY, eps);
}

} // namespace

std::vector<RCWAResult> runRCWA(const RCWAProblem& prob, std::string& err)
{
    std::vector<RCWAResult> results;

    if (prob.nHy != 0) {
        err = "現在のドライバは 1D 格子 (nHy=0) のみ対応しています";
        return results;
    }

    // z 層境界 (昇順) -> 上端(最大 z)から下端へ並べ替えてスタックを作る
    std::vector<scalar> zedges = deriveZEdges(prob);
    if (zedges.size() < 2) {
        err = "z 方向の層を構成できません (zmesh / geometry を確認してください)";
        return results;
    }

    std::vector<scalar> xgrid = buildCommonXGrid(prob);

    // 各スラブの中心 z と厚さ (上端から下端の順)
    std::vector<scalar> slabCenter, slabThick;
    for (int i = static_cast<int>(zedges.size()) - 1; i > 0; --i) {
        scalar zhi = zedges[i];
        scalar zlo = zedges[i - 1];
        slabCenter.push_back(0.5 * (zhi + zlo));
        slabThick.push_back(zhi - zlo);
    }
    const int nSlab = static_cast<int>(slabCenter.size());

    // 入射偏波 (面内 E ベクトル方向)。pol=1 -> x 偏波, pol=2 -> y 偏波。
    const scalar px = (prob.pol == 2) ? 0.0 : 1.0;
    const scalar py = (prob.pol == 2) ? 1.0 : 0.0;

    for (scalar lambda : prob.lambdas) {
        try {
            RCWASolver solver(prob.nHx, prob.nHy);
            solver.disablePML();  // 周期境界 (回折格子)

            for (int s = 0; s < nSlab; ++s) {
                Layer layer = makeLayer(prob, xgrid, slabCenter[s]);
                solver.addLayer(layer);
            }

            std::vector<int>    layerStack(nSlab);
            std::vector<scalar> thickness(nSlab);
            for (int s = 0; s < nSlab; ++s) {
                layerStack[s] = s;
                thickness[s]  = slabThick[s];
            }

            solver.solve(lambda, layerStack, thickness);

            // 入射側 = スタック上端 (layerStack[0]), 透過側 = 下端
            const int refLayer = 0;
            const int trnLayer = nSlab - 1;
            const scalar k0 = 2.0 * Pi / lambda;

            VectorXcs cInc;
            solver.generateHorizontalPlaneWave(px, py, refLayer, cInc);

            VectorXs REF, TRN;
            solver.scatterPlaneWave(cInc, refLayer, trnLayer, k0, REF, TRN);

            RCWAResult r;
            r.lambda = lambda;
            r.R = REF.sum();   // 全回折次数の電力反射率の和
            r.T = TRN.sum();   // 全回折次数の電力透過率の和
            results.push_back(r);
        }
        catch (const std::exception& e) {
            err = std::string("RCWA 計算中に例外 (lambda=") +
                  std::to_string(lambda) + "): " + e.what();
            return results;
        }
    }

    return results;
}
