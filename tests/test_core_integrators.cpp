// test_core_integrators.cpp — core/ の積分器 (RCWAIntegrator / DifferentialIntegrator /
// RedhefferIntegrator) の短時間回帰テスト。
//
// test_cpu は 8192 区間の RCWA を正解データとして作る研究用ベンチマーク
// (コメントによれば 1 日がかり) で、判定も閾値も持たない。ここでは小さい問題
// (7×5 高調波) で、正解データ無しに成り立つ関係を数秒で確かめる:
//
//   1. 解析解: z 方向に一様な区間の S 行列は、両端の固有基底で
//      Tuu = Tdd = diag(exp(i·γ·L/k0)), Rud = Rdu = 0 になる (区間数によらない)。
//      γ は k0² でスケールした PQ の固有値の平方根 (rcwa/WaveEquationSolver)、
//      伝搬長は L/k0 に正規化される (rcwa/RCWAScatterMatrix の normalizedLz_)。
//   2. 収束: x 方向の境界が z とともに動く区間で、RCWAIntegrator の区間数を
//      増やすと Tuu が収束する (誤差が区間数に対して単調に減る)。
//   3. 積分器の一致: DifferentialIntegrator と RedhefferIntegrator (VarRCWA) が
//      十分細かい RCWAIntegrator と一致する。
//
// 終了コード = 失敗数 (test_rcwa_input と同じ規約)。

#include <cmath>
#include <complex>
#include <iostream>
#include <memory>

#include <Eigen/Core>
#include "rcwa/defns.h"
#include "core/VarLayerSampler.h"
#include "core/RCWAIntegrator.h"
#include "core/DifferentialIntegrator.h"
#include "core/RedhefferIntegrator.h"

using namespace Eigen;

namespace {

int g_failures = 0;

#define CHECK(cond, msg)                                                   \
    do {                                                                   \
        if (cond) {                                                        \
            std::cout << "  ok  " << msg << "\n";                          \
        } else {                                                           \
            std::cout << "  NG  " << msg << "\n";                          \
            ++g_failures;                                                  \
        }                                                                  \
    } while (0)

const scalar Z0 = -0.1, Z1 = 1.1;      // 区間 [μm]
const int NHX = 3, NHY = 2;            // 高調波の半数 (7 × 5)
const scalar LAMBDA = 1.55;            // 波長 [μm]

// test_cpu と同じ導波路断面 (クラッド中のシリコンコア) を小さくしたもの。
// moving = true なら x[3] (コア右端) を z ∈ [0, 1] で 2.6 → 3.7 に動かす。
std::shared_ptr<VarLayerSampler> makeSampler(bool moving)
{
    const scalar cl = 1.445 * 1.445, co = 3.48 * 3.48;
    VectorXs x(8), y(7);
    x << .0, .5, 1.5, 2.6, 3.7, 3.9, 5., 5.5;
    y << .0, .5, 1.5, 1.58, 1.8, 2.8, 3.3;
    MatrixXcs eps(6, 7);
    eps.setConstant(cl);
    eps(3, 2) = co;
    auto s = std::make_shared<VarLayerSampler>(LAMBDA, 2 * NHX + 1, 2 * NHY + 1, x, y, eps);
    if (moving) {
        s->addXRule(3, [](scalar z) {
            if (z < 0.) return scalar(2.6);
            if (z > 1.) return scalar(3.7);
            return (1. - z) * 2.6 + z * 3.7;
        });
    }
    return s;
}

scalar maxAbs(const MatrixXcs& m) { return m.cwiseAbs().maxCoeff(); }

}  // namespace

int main()
{
    std::cout << "[1] 一様な区間: Tuu = diag(exp(i·γ·L/k0)), Rud = 0\n";
    {
        auto s = makeSampler(false);
        const scalar L = (Z1 - Z0) * LAMBDA / (2 * Pi);   // L / k0
        for (int N : {1, 8}) {
            RCWAIntegrator in(Z0, Z1, s);
            in.compute(N);
            const VectorXcs& g = in.gam0();
            MatrixXcs expect = MatrixXcs::Zero(g.size(), g.size());
            for (int i = 0; i < g.size(); ++i)
                expect(i, i) = std::exp(scalex(0, 1) * g(i) * L);
            const scalar eT = std::max(maxAbs(in.Tuu() - expect), maxAbs(in.Tdd() - expect));
            const scalar eR = std::max(maxAbs(in.Rud()), maxAbs(in.Rdu()));
            CHECK(eT < 1e-8, "RCWA N=" << N << ": max|Tuu - exp(iγL/k0)| = " << eT);
            CHECK(eR < 1e-8, "RCWA N=" << N << ": max|Rud| = " << eR);
        }
        // Differential は数値積分なので区間数とともに解析解へ収束する
        scalar prevD = 1e300;
        bool monoD = true;
        for (int N : {64, 128, 256}) {
            DifferentialIntegrator d(Z0, Z1, s);
            d.compute(N);
            const VectorXcs& g = d.gam0();
            MatrixXcs expect = MatrixXcs::Zero(g.size(), g.size());
            for (int i = 0; i < g.size(); ++i)
                expect(i, i) = std::exp(scalex(0, 1) * g(i) * L);
            const scalar eT = maxAbs(d.Tuu() - expect);
            std::cout << "      Differential N=" << N << ": max|Tuu - exp(iγL/k0)| = " << eT << "\n";
            if (eT >= prevD) monoD = false;
            prevD = eT;
        }
        CHECK(monoD, "Differential: 区間数を増やすと解析解との差が単調に減る");
        CHECK(prevD < 1e-5, "Differential N=256: 解析解との差 " << prevD << " (< 1e-5)");
    }

    std::cout << "[2] 境界が動く区間: RCWAIntegrator の区間数に対する収束\n";
    auto s = makeSampler(true);
    RCWAIntegrator ref(Z0, Z1, s);
    ref.compute(512);
    const MatrixXcs Tref = ref.Tuu();
    const scalar scale = maxAbs(Tref);
    scalar prev = 1e300, minRatio = 1e300;
    bool monotone = true;
    for (int N : {16, 32, 64, 128}) {
        RCWAIntegrator in(Z0, Z1, s);
        in.compute(N);
        const scalar e = maxAbs(in.Tuu() - Tref) / scale;
        std::cout << "      N=" << N << ": max|Tuu - Tuu(512)| / max|Tuu| = " << e << "\n";
        if (e >= prev) monotone = false;
        if (prev < 1e300) minRatio = std::min(minRatio, prev / e);
        prev = e;
    }
    CHECK(monotone, "区間数を増やすと誤差が単調に減る");
    // 階段近似の RCWA は区間幅の 2 乗で収束する (実測で区間数 2 倍ごとに約 1/4)
    CHECK(minRatio > 3.0, "2 次の収束: 区間数 2 倍あたりの誤差比の最小値 " << minRatio << " (> 3)");
    CHECK(prev < 1e-4, "N=128 の相対誤差 " << prev << " (< 1e-4)");

    std::cout << "[3] 積分器の一致 (基準: RCWAIntegrator N=512)\n";
    {
        DifferentialIntegrator d(Z0, Z1, s);
        d.compute(512);
        const scalar e = maxAbs(d.Tuu() - Tref) / scale;
        CHECK(e < 1e-5, "Differential N=512: 相対誤差 " << e << " (< 1e-5)");

        RedhefferIntegrator r(Z0, Z1, s);
        r.compute(0, 100000, 1e-3);
        const scalar er = maxAbs(r.Tuu() - Tref) / scale;
        CHECK(er < 1e-3, "Redheffer (pe=1e-3, " << r.N() << " 区間): 相対誤差 " << er << " (< 1e-3)");
    }

    if (g_failures == 0) std::cout << "All tests PASSED\n";
    else std::cout << g_failures << " test(s) FAILED\n";
    return g_failures;
}
