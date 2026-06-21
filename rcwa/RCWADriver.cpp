#include "rcwa/RCWADriver.h"
#include "rcwa/RCWASolver.h"
#include "rcwa/Layer.h"

#include <algorithm>
#include <set>
#include <cmath>
#include <iostream>

using namespace Eigen;

namespace {

// 光速 [μm/s] と真空誘電率 [F/m] — 複素誘電率の導電率項計算に使用。
constexpr scalar C0_UM  = 2.99792458e14;   // c [μm/s]
constexpr scalar EPS0   = 8.854187817e-12; // ε₀ [F/m]

// 全スラブ共通の x エッジを収集する。
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

// 全スラブ共通の y エッジを収集する。
std::vector<scalar> buildCommonYGrid(const RCWAProblem& prob)
{
    std::set<scalar> ys;
    ys.insert(prob.ymin);
    ys.insert(prob.ymax);
    for (const auto& b : prob.boxes) {
        scalar y0 = std::max(b.y0, prob.ymin);
        scalar y1 = std::min(b.y1, prob.ymax);
        if (y0 > prob.ymin && y0 < prob.ymax) ys.insert(y0);
        if (y1 > prob.ymin && y1 < prob.ymax) ys.insert(y1);
    }
    return std::vector<scalar>(ys.begin(), ys.end());
}

// 点 (xc, yc, zc) における複素誘電率。
// lambda [μm] は損失材料 (esgm≠0) の導電率項の計算に必要。
scalex epsAt(const RCWAProblem& prob,
             scalar xc, scalar yc, scalar zc,
             scalar lambda)
{
    scalex e = prob.backgroundEps;
    for (const auto& b : prob.boxes) {
        if (zc < b.z0 || zc > b.z1) continue;
        if (!b.containsXY(xc, yc)) continue;
        int m = b.material;
        if (m < 0 || m >= static_cast<int>(prob.materialEps.size())) continue;
        e = prob.materialEps[m];
        // 導電率による損失項: Δε_im = -σ/(ω·ε₀)
        // Lorentz 分散モデル: eps(omega) = einf + ae^2/(ce^2-omega^2-i*be*omega)
        if (!prob.materialDispersion.empty() &&
            m < static_cast<int>(prob.materialDispersion.size())) {
            const auto& d = prob.materialDispersion[m];
            if (d.ae != 0.0) {
                scalar omega = 2.0 * Pi * C0_UM / lambda;
                scalex denom = scalex(d.ce*d.ce - omega*omega, -d.be*omega);
                e = scalex(d.einf, 0.0) + scalex(d.ae*d.ae, 0.0) / denom;
                return e;  // dispersion overrides static eps; skip conductivity term
            }
        }
        // 導電率による損失項: Δε_im = -σ/(ω·ε₀)
        if (!prob.materialSigma.empty() && m < static_cast<int>(prob.materialSigma.size())) {
            scalar sigma = prob.materialSigma[m];
            if (sigma != 0.0) {
                scalar omega = 2.0 * Pi * C0_UM / lambda; // [rad/s]
                e += scalex(0.0, -sigma / (omega * EPS0));
            }
        }
    }
    return e;
}

// xgrid×ygrid セル格子と z スラブ中央 zc から Layer を生成する。
// eps 行列は (nCellY 行 × nCellX 列) で Layer が要求する形式。
Layer makeLayer(const RCWAProblem& prob,
                const std::vector<scalar>& xgrid,
                const std::vector<scalar>& ygrid,
                scalar zc, scalar lambda)
{
    const int nCellX = static_cast<int>(xgrid.size()) - 1;
    const int nCellY = static_cast<int>(ygrid.size()) - 1;

    VectorXs coordX(xgrid.size());
    for (size_t i = 0; i < xgrid.size(); ++i) coordX[i] = xgrid[i];

    VectorXs coordY(ygrid.size());
    for (size_t j = 0; j < ygrid.size(); ++j) coordY[j] = ygrid[j];

    MatrixXcs eps(nCellY, nCellX);
    for (int j = 0; j < nCellY; ++j) {
        scalar yc = 0.5 * (ygrid[j] + ygrid[j + 1]);
        for (int i = 0; i < nCellX; ++i) {
            scalar xc = 0.5 * (xgrid[i] + xgrid[i + 1]);
            eps(j, i) = epsAt(prob, xc, yc, zc, lambda);
        }
    }

    return Layer(coordX, coordY, eps);
}

} // namespace

std::vector<RCWAResult> runRCWA(const RCWAProblem& prob, std::string& err)
{
    std::vector<RCWAResult> results;

    // z 層境界 (昇順) -> 上端から下端の順にスタックを作る
    std::vector<scalar> zedges = deriveZEdges(prob);
    if (zedges.size() < 2) {
        err = "z 方向の層を構成できません (zmesh / geometry を確認してください)";
        return results;
    }

    std::vector<scalar> xgrid = buildCommonXGrid(prob);
    std::vector<scalar> ygrid = buildCommonYGrid(prob);

    // 各スラブの中心 z と厚さ (上端から下端の順)
    std::vector<scalar> slabCenter, slabThick;
    for (int i = static_cast<int>(zedges.size()) - 1; i > 0; --i) {
        scalar zhi = zedges[i];
        scalar zlo = zedges[i - 1];
        slabCenter.push_back(0.5 * (zhi + zlo));
        slabThick.push_back(zhi - zlo);
    }
    const int nSlab = static_cast<int>(slabCenter.size());

    // 入射角 [rad]
    const scalar thetaRad = prob.theta * Pi / 180.0;
    const scalar phiRad   = prob.phi   * Pi / 180.0;

    // 入射 E 場の横 (xy) 成分ベクトル。単位振幅 (|E_inc_3D|² = 1) に正規化する。
    //   pol=1 (TM / p 偏波): E は入射面内 → 横成分 = (cosθ·cosφ, cosθ·sinφ)
    //   pol=2 (TE / s 偏波): E ⊥ 入射面   → 横成分 = (−sinφ, cosφ)
    // 法線入射 (θ≈0) では TE/TM が縮退するため単純な x/y 方向を使う。
    scalar px, py;
    if (std::abs(std::sin(thetaRad)) < 1e-9) {
        px = (prob.pol == 2) ? 0.0 : 1.0;
        py = (prob.pol == 2) ? 1.0 : 0.0;
    } else if (prob.pol == 2) {
        // TE: E ⊥ 入射面 → (-sinφ, cosφ)、z 成分なし
        px = -std::sin(phiRad);
        py =  std::cos(phiRad);
    } else {
        // TM: E ∥ 入射面 → 横成分 (cosθ·cosφ, cosθ·sinφ), Ez = −sinθ
        px = std::cos(thetaRad) * std::cos(phiRad);
        py = std::cos(thetaRad) * std::sin(phiRad);
    }

    for (scalar lambda : prob.lambdas) {
        try {
            RCWASolver solver(prob.nHx, prob.nHy);
            solver.disablePML();  // 周期境界 (回折格子)

            for (int s = 0; s < nSlab; ++s) {
                Layer layer = makeLayer(prob, xgrid, ygrid, slabCenter[s], lambda);
                solver.addLayer(layer);
            }

            // 入射側スラブ (最上段) の誘電率から Bloch 位相を決定
            const int refLayer = 0;
            const int trnLayer = nSlab - 1;
            const scalar k0    = 2.0 * Pi / lambda;

            scalex eps_inc = epsAt(prob,
                                   0.5 * (prob.xmin + prob.xmax),
                                   0.5 * (prob.ymin + prob.ymax),
                                   slabCenter[refLayer], lambda);
            const scalar n_inc = std::sqrt(std::max(eps_inc.real(), scalar(1.0)));
            const scalar kx0 = k0 * n_inc * std::sin(thetaRad) * std::cos(phiRad);
            const scalar ky0 = k0 * n_inc * std::sin(thetaRad) * std::sin(phiRad);
            solver.setBlochWavevector(kx0, ky0);

            std::vector<int>    layerStack(nSlab);
            std::vector<scalar> thickness(nSlab);
            for (int s = 0; s < nSlab; ++s) {
                layerStack[s] = s;
                thickness[s]  = slabThick[s];
            }

            solver.solve(lambda, layerStack, thickness);

            VectorXcs cInc;
            solver.generateHorizontalPlaneWave(px, py, refLayer, cInc);

            VectorXs REF, TRN;
            solver.scatterPlaneWave(cInc, refLayer, trnLayer, k0, REF, TRN);

            RCWAResult r;
            r.lambda = lambda;
            r.R = REF.sum();
            r.T = TRN.sum();
            r.A = 1.0 - r.R - r.T;
            r.REF_orders.assign(REF.data(), REF.data() + REF.size());
            r.TRN_orders.assign(TRN.data(), TRN.data() + TRN.size());
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
