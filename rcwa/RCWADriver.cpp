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
constexpr scalar MU0    = 1.25663706212e-6; // μ₀ [H/m]

// 曲面形状を階段近似するための面内分割数。
// 直方体はエッジだけで厳密に表現できるので分割しない。
constexpr int N_CURVED_CELLS = 8;

// 形状 b の断面が x / y 方向に曲率を持つか (= その軸を細分する必要があるか)。
//   SPHERE / CYL_Z : XY 断面が円 → x, y とも細分
//   CYL_Y          : XZ 断面が円 → x を細分 (y は軸方向で平坦)
//   CYL_X          : YZ 断面が円 → y を細分 (x は軸方向で平坦)
bool curvedInX(const RCWABox& b)
{
    return b.shape == RCWA_SHAPE_SPHERE || b.shape == RCWA_SHAPE_CYL_Z ||
           b.shape == RCWA_SHAPE_CYL_Y;
}

bool curvedInY(const RCWABox& b)
{
    return b.shape == RCWA_SHAPE_SPHERE || b.shape == RCWA_SHAPE_CYL_Z ||
           b.shape == RCWA_SHAPE_CYL_X;
}

// 区間 [lo, hi] を n 等分する内部エッジを、[cmin, cmax] に収まる範囲で追加する。
void insertSubdivisions(std::set<scalar>& out, scalar lo, scalar hi,
                        scalar cmin, scalar cmax, int n)
{
    for (int i = 1; i < n; ++i) {
        scalar v = lo + (hi - lo) * i / static_cast<scalar>(n);
        if (v > cmin && v < cmax) out.insert(v);
    }
}

// 全スラブ共通の x エッジを収集する。
// 曲面形状はバウンディングボックスのエッジだけでは断面を表現できないため
// 内部も細分する。細分しないとセル中心が数点しか取れず、球・円柱・直方体が
// まったく同じ誘電率分布になってしまう (形状指定が無視されるのと同じ)。
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
        if (curvedInX(b) && x0 < x1)
            insertSubdivisions(xs, x0, x1, prob.xmin, prob.xmax, N_CURVED_CELLS);
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
        if (curvedInY(b) && y0 < y1)
            insertSubdivisions(ys, y0, y1, prob.ymin, prob.ymax, N_CURVED_CELLS);
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
        // z 依存断面 (球・x/y 軸円柱の z ごとの半径変化) を含めて判定する
        if (!b.containsXYZ(xc, yc, zc)) continue;
        int m = b.material;
        if (m < 0 || m >= static_cast<int>(prob.materialEps.size())) continue;
        e = prob.materialEps[m];

        const scalar omega = 2.0 * Pi * C0_UM / lambda;  // [rad/s]

        // 多極 Lorentz 分散: eps = einf + Σ_p ae_p²/(ce_p²-ω²-i·be_p·ω)
        // 分散が指定された材料では静的 eps ではなくこの式を用いる。
        if (m < static_cast<int>(prob.materialDispersion.size())) {
            const auto& d = prob.materialDispersion[m];
            if (d.active()) {
                e = scalex(d.einf, 0.0);
                for (const auto& p : d.poles)
                    e += scalex(p.ae * p.ae, 0.0)
                       / scalex(p.ce * p.ce - omega * omega, -p.be * omega);
            }
        }

        // 導電率による損失項: Δε_im = +σ/(ω·ε₀)
        // 本ソルバの時間規約 exp(-iωt) では損失媒質は正の虚部を持つ
        // (Lorentz 分散式 ce²-ω²-i·be·ω が正虚部を返すのと同一規約)。
        // 分散極とは独立に加算する (両方指定すれば両方効く)。
        if (m < static_cast<int>(prob.materialSigma.size())) {
            scalar sigma = prob.materialSigma[m];
            if (sigma != 0.0)
                e += scalex(0.0, sigma / (omega * EPS0));
        }
    }
    return e;
}

// 点 (xc, yc, zc) における複素比透磁率。epsAt と同じ形状判定を使う。
// 磁気導電率 msgm があれば損失として虚部に加える (exp(-iωt) 規約で正)。
scalex muAt(const RCWAProblem& prob,
            scalar xc, scalar yc, scalar zc,
            scalar lambda)
{
    scalex m = scalex(1.0, 0.0);   // 背景は非磁性
    for (const auto& b : prob.boxes) {
        if (!b.containsXYZ(xc, yc, zc)) continue;
        int mi = b.material;
        if (mi < 0 || mi >= static_cast<int>(prob.materialMu.size())) continue;
        m = prob.materialMu[mi];

        if (mi < static_cast<int>(prob.materialMagSigma.size())) {
            scalar msig = prob.materialMagSigma[mi];
            if (msig != 0.0) {
                scalar omega = 2.0 * Pi * C0_UM / lambda;  // [rad/s]
                m += scalex(0.0, msig / (omega * MU0));
            }
        }
    }
    return m;
}

// この問題に磁性材料が含まれるか (含まれなければ μ 行列を作らず高速経路)。
bool hasMagneticMaterial(const RCWAProblem& prob)
{
    for (const auto& m : prob.materialMu)
        if (std::abs(m - scalex(1.0, 0.0)) > 1e-12) return true;
    for (scalar s : prob.materialMagSigma)
        if (s != 0.0) return true;
    return false;
}

// xgrid×ygrid セル格子と z スラブ中央 zc から Layer を生成する。
// eps 行列は (nCellY 行 × nCellX 列) で Layer が要求する形式。
Layer makeLayer(const RCWAProblem& prob,
                const std::vector<scalar>& xgrid,
                const std::vector<scalar>& ygrid,
                scalar zc, scalar lambda, bool magnetic)
{
    const int nCellX = static_cast<int>(xgrid.size()) - 1;
    const int nCellY = static_cast<int>(ygrid.size()) - 1;

    VectorXs coordX(xgrid.size());
    for (size_t i = 0; i < xgrid.size(); ++i) coordX[i] = xgrid[i];

    VectorXs coordY(ygrid.size());
    for (size_t j = 0; j < ygrid.size(); ++j) coordY[j] = ygrid[j];

    MatrixXcs eps(nCellY, nCellX);
    MatrixXcs mu(nCellY, nCellX);
    for (int j = 0; j < nCellY; ++j) {
        scalar yc = 0.5 * (ygrid[j] + ygrid[j + 1]);
        for (int i = 0; i < nCellX; ++i) {
            scalar xc = 0.5 * (xgrid[i] + xgrid[i + 1]);
            eps(j, i) = epsAt(prob, xc, yc, zc, lambda);
            mu(j, i)  = magnetic ? muAt(prob, xc, yc, zc, lambda)
                                 : scalex(1.0, 0.0);
        }
    }

    // 非磁性なら μ を渡さない (Layer 側で高速経路が選ばれる)
    return magnetic ? Layer(coordX, coordY, eps, mu)
                    : Layer(coordX, coordY, eps);
}

} // namespace

std::vector<RCWAResult> runRCWA(const RCWAProblem& prob, std::string& err,
                                const RCWAFieldRequest* fieldReq)
{
    std::vector<RCWAResult> results;
    bool fieldSaved = false;

    // z 層境界 (昇順) -> 上端から下端の順にスタックを作る
    std::vector<scalar> zedges = deriveZEdges(prob);
    if (zedges.size() < 2) {
        err = "z 方向の層を構成できません (zmesh / geometry を確認してください)";
        return results;
    }

    // 場イメージは半無限層に挟まれた有限厚の内部層に対してのみ定義される。
    // 層数 = zedges.size()-1 なので、内部層を持つには 3 層以上が必要。
    if (fieldReq && fieldReq->wantsAnything() && zedges.size() - 1 < 3) {
        err = "断面出力には内部層が必要です "
              "(zmesh / geometry を 3 層以上に分割してください)";
        return results;
    }

    std::vector<scalar> xgrid = buildCommonXGrid(prob);
    std::vector<scalar> ygrid = buildCommonYGrid(prob);

    // 磁性材料が一切なければ μ 行列を作らず、Layer/P/Q の高速経路を通す。
    const bool magnetic = hasMagneticMaterial(prob);

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
    // 直交する単位偏波基底の横成分:
    //   TM (p 偏波): E は入射面内 → (cosθ·cosφ, cosθ·sinφ)   [Ez = −sinθ]
    //   TE (s 偏波): E ⊥ 入射面   → (−sinφ, cosφ)            [Ez = 0]
    // 法線入射 (θ≈0) では TE/TM が縮退するため単純な x/y 方向を使う。
    scalar tmX, tmY, teX, teY;
    if (std::abs(std::sin(thetaRad)) < 1e-9) {
        tmX = 1.0; tmY = 0.0;
        teX = 0.0; teY = 1.0;
    } else {
        tmX = std::cos(thetaRad) * std::cos(phiRad);
        tmY = std::cos(thetaRad) * std::sin(phiRad);
        teX = -std::sin(phiRad);
        teY =  std::cos(phiRad);
    }

    // 偏波種別ごとの複素展開係数 (a_TM, a_TE)。|a_TM|² + |a_TE|² = 1 に保つ。
    // 基底は直交するので、この規格化で |E_inc_3D| = 1 が維持される。
    const scalar invSqrt2 = 1.0 / std::sqrt(scalar(2.0));
    scalex aTM, aTE;
    switch (prob.pol) {
    case 2:  // TE (s)
        aTM = 0.0; aTE = 1.0;
        break;
    case 3: { // psi [deg] の直線偏波 (0=TM, 90=TE)
        const scalar psiRad = prob.psi * Pi / 180.0;
        aTM = std::cos(psiRad);
        aTE = std::sin(psiRad);
        break;
    }
    case 4:  // 右円偏波
        aTM = invSqrt2; aTE = scalex(0.0, invSqrt2);
        break;
    case 5:  // 左円偏波
        aTM = invSqrt2; aTE = scalex(0.0, -invSqrt2);
        break;
    case 1:
    default: // TM (p)
        aTM = 1.0; aTE = 0.0;
        break;
    }

    const scalex px = aTM * tmX + aTE * teX;
    const scalex py = aTM * tmY + aTE * teY;

    for (scalar lambda : prob.lambdas) {
        try {
            RCWASolver solver(prob.nHx, prob.nHy);
            solver.disablePML();  // 周期境界 (回折格子)

            for (int s = 0; s < nSlab; ++s) {
                Layer layer = makeLayer(prob, xgrid, ygrid, slabCenter[s],
                                        lambda, magnetic);
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
            // 入射媒質の屈折率。複素平方根の実部を用いることで eps<1 (プラズマ等)
            // にも対応する。実部を 1 でクランプしていた旧実装は eps<1 で入射角を
            // 誤り、Bloch 波数がずれていた。
            scalex mu_inc = magnetic
                ? muAt(prob, 0.5 * (prob.xmin + prob.xmax),
                       0.5 * (prob.ymin + prob.ymax),
                       slabCenter[refLayer], lambda)
                : scalex(1.0, 0.0);
            // 屈折率は n = sqrt(ε·μ)。磁性入射媒質でも角度が正しくなる。
            const scalar n_inc =
                std::max(std::sqrt(eps_inc * mu_inc).real(), scalar(1e-6));
            // 入射半無限媒質が損失を持つと R/T の規格化 (実 kz 前提) が崩れる。
            if (std::abs(eps_inc.imag()) > 1e-6 * std::abs(eps_inc.real())) {
                static bool warned = false;
                if (!warned) {
                    std::cerr << "*** 警告: 入射側媒質に損失があります"
                                 " (R/T は無損失入射を前提とした規格化です)\n";
                    warned = true;
                }
            }
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

            // 断面出力 (最初の波長のみ)
            if (fieldReq && !fieldSaved && fieldReq->wantsAnything()) {
                if (!fieldReq->path.empty()) {
                    std::vector<std::pair<int, scalex>> inputCoeffs;
                    for (int i = 0; i < cInc.size(); ++i)
                        if (std::abs(cInc(i)) > 1e-14)
                            inputCoeffs.emplace_back(i, cInc(i));
                    solver.saveFieldImage(fieldReq->path,
                                          fieldReq->slice, fieldReq->coord,
                                          fieldReq->component, fieldReq->opt,
                                          inputCoeffs, layerStack, thickness);
                }
                if (!fieldReq->devicePath.empty()) {
                    solver.saveDeviceImage(fieldReq->devicePath,
                                           fieldReq->slice, fieldReq->coord,
                                           layerStack, thickness);
                }
                fieldSaved = true;
            }

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
