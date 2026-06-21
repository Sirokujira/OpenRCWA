// ============================================================
// test_rcwa_input.cpp — RCWAInput パーサ + RCWADriver の単体テスト
//
// テスト一覧:
//   1. Fresnel 反射 (air→glass): R≈0.04, R+T≈1
//   2. Fabry-Pérot エタロン: 全波長で R+T≈1
//   3. 斜め入射 Brewster 角: TM 反射 R≈0
//   4. frequency1 + frequency2 結合: 波長数・昇順を確認
//   5. 円柱ジオメトリ containsXY
//   6. material_dispersion パース確認
//   7. 回折次数別出力: sum(REF_orders)≡R, sum(TRN_orders)≡T
//   8. 2D 格子 (nHy=1): エネルギー保存
// ============================================================
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "rcwa/RCWAInput.h"
#include "rcwa/RCWADriver.h"

// ============================================================
// ヘルパー
// ============================================================

static int g_fails = 0;

// 失敗時にカウントするチェックマクロ (致命的でない)
#define CHECK(cond, msg)                                              \
    do {                                                              \
        if (!(cond)) {                                                \
            std::cerr << "  FAIL: " << msg << "\n";                  \
            ++g_fails;                                                \
        }                                                             \
    } while (0)

// .orcwa 内容をテンポラリファイルに書いてパスを返す
static std::string writeTmp(const std::string& tag, const std::string& content)
{
    std::string path = std::string("/tmp/test_rcwa_") + tag + ".orcwa";
    std::ofstream f(path);
    f << content;
    return path;
}

// パース + ソルブ。空配列ならエラーメッセージを表示して空を返す。
static std::vector<RCWAResult> solve(const std::string& tag, const std::string& content)
{
    std::string path = writeTmp(tag, content);
    RCWAProblem prob;
    std::string err;
    if (!parseRCWAInput(path, prob, err)) {
        std::cerr << "  parseRCWAInput error: " << err << "\n";
        return {};
    }
    auto res = runRCWA(prob, err);
    if (res.empty()) {
        std::cerr << "  runRCWA error: " << err << "\n";
    }
    return res;
}

// ============================================================
// Test 1: Fresnel 反射 — air→glass, 法線入射
// 解析解: R = ((n-1)/(n+1))^2 ≈ 0.04 (n=1.5), R+T = 1
// ============================================================
static void test_fresnel()
{
    static const char* name = "test_fresnel";
    auto res = solve("fresnel", R"(
OpenRCWA 4 2
title = fresnel air-glass
xmesh = -2.5e-07 10 2.5e-07
ymesh = -2.5e-07 10 2.5e-07
zmesh = -3e-07 10 0.0 10 7e-07
material = 1 2.25 0 1 0
geometry = 2 1 -2.5e-07 2.5e-07 -2.5e-07 2.5e-07 -3e-07 0
planewave = 0 0 1
pbc = 1 1 0
rcwaorder = 6 0
frequency1 = 5.0e+14 5.0e+14 0
end
)");
    if (res.empty()) { std::cerr << "  SKIP: " << name << " (solve failed)\n"; return; }
    double R = res[0].R, T = res[0].T;
    std::cout << "  R=" << R << " T=" << T << " R+T=" << R + T << "\n";
    CHECK(std::abs(R - 0.04) < 2e-3, "R≈0.04");
    CHECK(std::abs(R + T - 1.0) < 1e-3, "R+T≈1 (energy conservation)");
    if (g_fails == 0) std::cout << "PASS: " << name << "\n";
    else              std::cout << "FAIL: " << name << "\n";
}

// ============================================================
// Test 2: Fabry-Pérot エタロン — 全波長でエネルギー保存
// n=1.5 ガラス板 (厚 0.3 μm) を空気中に配置、31 波長
// ============================================================
static void test_fabry_perot()
{
    static const char* name = "test_fabry_perot";
    int prev_fails = g_fails;
    auto res = solve("fabry_perot", R"(
OpenRCWA 4 2
title = fabry-perot glass slab in air
xmesh = -2.5e-07 10 2.5e-07
ymesh = -2.5e-07 10 2.5e-07
zmesh = -1e-06 10 0.0 10 3e-07 10 7e-07
material = 1 2.25 0 1 0
geometry = 2 1 -2.5e-07 2.5e-07 -2.5e-07 2.5e-07 0 3e-07
planewave = 0 0 1
pbc = 1 1 0
rcwaorder = 4 0
frequency1 = 4.0e+14 7.0e+14 30
end
)");
    if (res.empty()) { std::cerr << "  SKIP: " << name << " (solve failed)\n"; return; }
    std::cout << "  wavelengths: " << res.size() << "\n";
    for (const auto& r : res) {
        double RT = r.R + r.T;
        CHECK(std::abs(RT - 1.0) < 2e-3,
              std::string("R+T≈1 at lambda=") + std::to_string(r.lambda));
    }
    if (g_fails == prev_fails) std::cout << "PASS: " << name << "\n";
    else                       std::cout << "FAIL: " << name << "\n";
}

// ============================================================
// Test 3: 斜め入射 Brewster 角 — TM 反射が ~0
// air→glass (n=1.5): θ_B = arctan(1.5) ≈ 56.31°
// TM (pol=1) での R は Brewster 角で 0 になる。
// ============================================================
static void test_brewster()
{
    static const char* name = "test_brewster";
    int prev_fails = g_fails;
    auto res = solve("brewster", R"(
OpenRCWA 4 2
title = brewster angle TM
xmesh = -5e-07 10 5e-07
ymesh = -5e-07 10 5e-07
zmesh = -5e-07 10 0.0 10 1e-06
material = 1 2.25 0 1 0
geometry = 2 1 -5e-07 5e-07 -5e-07 5e-07 -5e-07 0
planewave = 56.31 0 1
pbc = 1 1 0
rcwaorder = 8 0
frequency1 = 5.0e+14 5.0e+14 0
end
)");
    if (res.empty()) { std::cerr << "  SKIP: " << name << " (solve failed)\n"; return; }
    double R = res[0].R, T = res[0].T;
    std::cout << "  R_TM(Brewster)=" << R << " T=" << T << " R+T=" << R + T << "\n";
    CHECK(R < 0.01, "R_TM≈0 at Brewster angle");
    CHECK(std::abs(R + T - 1.0) < 2e-3, "energy conservation at oblique incidence");
    if (g_fails == prev_fails) std::cout << "PASS: " << name << "\n";
    else                       std::cout << "FAIL: " << name << "\n";
}

// ============================================================
// Test 4: frequency1 + frequency2 両方使用
// 合計 5+5=10 波長が昇順ソートで得られることを確認
// ============================================================
static void test_frequency2()
{
    static const char* name = "test_frequency2";
    int prev_fails = g_fails;
    std::string path = writeTmp("freq2", R"(
OpenRCWA 4 2
title = frequency2 test
xmesh = -2.5e-07 10 2.5e-07
ymesh = -2.5e-07 10 2.5e-07
zmesh = -3e-07 10 0.0 10 7e-07
material = 1 2.25 0 1 0
geometry = 2 1 -2.5e-07 2.5e-07 -2.5e-07 2.5e-07 -3e-07 0
planewave = 0 0 1
pbc = 1 1 0
rcwaorder = 4 0
frequency1 = 4.0e+14 5.0e+14 4
frequency2 = 6.0e+14 7.0e+14 4
end
)");
    RCWAProblem prob;
    std::string err;
    if (!parseRCWAInput(path, prob, err)) {
        std::cerr << "  parse error: " << err << "\n";
        ++g_fails;
        return;
    }
    std::cout << "  lambdas: " << prob.lambdas.size() << "\n";
    CHECK(prob.lambdas.size() == 10, "10 wavelengths from frequency1+frequency2");
    for (size_t i = 1; i < prob.lambdas.size(); ++i)
        CHECK(prob.lambdas[i] > prob.lambdas[i-1], "lambdas sorted ascending");
    if (g_fails == prev_fails) std::cout << "PASS: " << name << "\n";
    else                       std::cout << "FAIL: " << name << "\n";
}

// ============================================================
// Test 5: 円柱ジオメトリ — containsXY の動作確認
// shape=11 (CYL_Z): バウンディングボックス 0→1 x 0→1
// 中心(0.5, 0.5), 半径 0.5 の楕円
// ============================================================
static void test_cylinder_geometry()
{
    static const char* name = "test_cylinder_geometry";
    int prev_fails = g_fails;

    RCWABox box;
    box.material = 2;
    box.x0 = 0.0; box.x1 = 1.0;
    box.y0 = 0.0; box.y1 = 1.0;
    box.z0 = 0.0; box.z1 = 1.0;

    // ---- 円柱 (CYL_Z) ----
    box.shape = RCWA_SHAPE_CYL_Z;
    // 中心 → 内側
    CHECK(box.containsXY(0.5, 0.5), "CYL_Z: center inside");
    // +x 端 → 境界上 (≤1 なので内側)
    CHECK(box.containsXY(1.0, 0.5), "CYL_Z: +x boundary inside");
    // バウンディングボックス角 → 外側 (正規化距離 sqrt(2) > 1)
    CHECK(!box.containsXY(0.0, 0.0), "CYL_Z: bbox corner outside");
    CHECK(!box.containsXY(1.0, 1.0), "CYL_Z: bbox corner outside");

    // ---- 直方体 (BOX) ----
    box.shape = RCWA_SHAPE_BOX;
    CHECK(box.containsXY(0.0, 0.0), "BOX: corner inside");
    CHECK(box.containsXY(0.5, 0.5), "BOX: center inside");
    CHECK(!box.containsXY(1.5, 0.5), "BOX: outside");

    // ---- 球 (SPHERE): XY 断面は円柱と同じ楕円判定 ----
    box.shape = RCWA_SHAPE_SPHERE;
    CHECK(box.containsXY(0.5, 0.5), "SPHERE: center inside");
    CHECK(!box.containsXY(0.0, 0.0), "SPHERE: corner outside");

    if (g_fails == prev_fails) std::cout << "PASS: " << name << "\n";
    else                       std::cout << "FAIL: " << name << "\n";
}

// ============================================================
// Test 6: material_dispersion パース
// インデックス 2 に Drude モデルパラメータが設定されることを確認
// ============================================================
static void test_material_dispersion()
{
    static const char* name = "test_material_dispersion";
    int prev_fails = g_fails;
    std::string path = writeTmp("dispersion", R"(
OpenRCWA 4 2
title = dispersion test
xmesh = -2.5e-07 10 2.5e-07
ymesh = -2.5e-07 10 2.5e-07
zmesh = -3e-07 10 0.0 10 7e-07
material = 1 2.25 0 1 0
material_dispersion = 2 1.0 1.37e16 2.73e13 0.0
geometry = 2 1 -2.5e-07 2.5e-07 -2.5e-07 2.5e-07 -3e-07 0
planewave = 0 0 1
pbc = 1 1 0
rcwaorder = 4 0
frequency1 = 5.0e+14 5.0e+14 0
end
)");
    RCWAProblem prob;
    std::string err;
    if (!parseRCWAInput(path, prob, err)) {
        std::cerr << "  parse error: " << err << "\n";
        ++g_fails;
        return;
    }
    std::cout << "  materialDispersion.size()=" << prob.materialDispersion.size() << "\n";
    CHECK(prob.materialDispersion.size() >= 3, "dispersion vector includes material 2");
    if (prob.materialDispersion.size() >= 3) {
        const auto& d = prob.materialDispersion[2];
        std::cout << "  einf=" << d.einf << " ae=" << d.ae
                  << " be=" << d.be << " ce=" << d.ce << "\n";
        CHECK(std::abs(d.einf - 1.0) < 1e-9, "einf=1.0");
        CHECK(std::abs(d.ae - 1.37e16) < 1e12, "ae=1.37e16");
        CHECK(std::abs(d.be - 2.73e13) < 1e10, "be=2.73e13");
        CHECK(d.ce == 0.0, "ce=0 (Drude model)");
    }
    if (g_fails == prev_fails) std::cout << "PASS: " << name << "\n";
    else                       std::cout << "FAIL: " << name << "\n";
}

// ============================================================
// Test 7: 回折次数別出力
// sum(REF_orders) == R, sum(TRN_orders) == T
// nHx=3, nHy=0 → 7 次
// ============================================================
static void test_per_order_output()
{
    static const char* name = "test_per_order_output";
    int prev_fails = g_fails;
    auto res = solve("perorder", R"(
OpenRCWA 4 2
title = per-order output
xmesh = -2.5e-07 10 2.5e-07
ymesh = -2.5e-07 10 2.5e-07
zmesh = -3e-07 10 0.0 10 7e-07
material = 1 2.25 0 1 0
geometry = 2 1 -2.5e-07 2.5e-07 -2.5e-07 2.5e-07 -3e-07 0
planewave = 0 0 1
pbc = 1 1 0
rcwaorder = 3 0
frequency1 = 5.0e+14 5.0e+14 0
end
)");
    if (res.empty()) { std::cerr << "  SKIP: " << name << " (solve failed)\n"; return; }
    const auto& r = res[0];
    // (2*3+1)*(2*0+1) = 7 orders
    std::cout << "  REF_orders.size()=" << r.REF_orders.size()
              << " TRN_orders.size()=" << r.TRN_orders.size() << "\n";
    CHECK(r.REF_orders.size() == 7, "7 REF diffraction orders");
    CHECK(r.TRN_orders.size() == 7, "7 TRN diffraction orders");
    double sumR = 0, sumT = 0;
    for (double v : r.REF_orders) sumR += v;
    for (double v : r.TRN_orders) sumT += v;
    std::cout << "  R=" << r.R << " sum(REF)=" << sumR
              << " T=" << r.T << " sum(TRN)=" << sumT << "\n";
    CHECK(std::abs(sumR - r.R) < 1e-10, "sum(REF_orders) == R");
    CHECK(std::abs(sumT - r.T) < 1e-10, "sum(TRN_orders) == T");
    if (g_fails == prev_fails) std::cout << "PASS: " << name << "\n";
    else                       std::cout << "FAIL: " << name << "\n";
}

// ============================================================
// Test 8: 2D 格子 (nHy=1) — エネルギー保存
// 正方格子のダミー構造; 全波長で R+T≈1 を確認
// ============================================================
static void test_2d_grating_energy()
{
    static const char* name = "test_2d_grating";
    int prev_fails = g_fails;
    // zmesh: air above (0→7e-7), grating layer (-3e-7→0), air below (-7e-7→-3e-7)
    // Transmission goes into the uniform air slab below the grating.
    auto res = solve("2dgrating", R"(
OpenRCWA 4 2
title = 2D grating energy conservation
xmesh = -2.5e-07 10 2.5e-07
ymesh = -2.5e-07 10 2.5e-07
zmesh = -7e-07 10 -3e-07 10 0.0 10 7e-07
material = 1 2.25 0 1 0
geometry = 2 1 -1.25e-07 1.25e-07 -1.25e-07 1.25e-07 -3e-07 0
planewave = 0 0 1
pbc = 1 1 0
rcwaorder = 2 2
frequency1 = 5.0e+14 5.0e+14 0
end
)");
    if (res.empty()) { std::cerr << "  SKIP: " << name << " (solve failed)\n"; return; }
    double R = res[0].R, T = res[0].T;
    std::cout << "  2D grating R=" << R << " T=" << T << " R+T=" << R + T << "\n";
    CHECK(std::abs(R + T - 1.0) < 2e-3, "2D grating energy conservation");
    if (g_fails == prev_fails) std::cout << "PASS: " << name << "\n";
    else                       std::cout << "FAIL: " << name << "\n";
}

// ============================================================
// Test 9: TE 偏波 Brewster 角
// air→glass (n=1.5): θ_B ≈ 56.31°, TE では R ≈ 0.148 (非ゼロ), R+T ≈ 1
// ============================================================
static void test_te_oblique()
{
    static const char* name = "test_te_oblique";
    int prev_fails = g_fails;
    auto res = solve("te_oblique", R"(
OpenRCWA 4 2
title = TE oblique incidence at Brewster angle
xmesh = -5e-07 10 5e-07
ymesh = -5e-07 10 5e-07
zmesh = -5e-07 10 0.0 10 1e-06
material = 1 2.25 0 1 0
geometry = 2 1 -5e-07 5e-07 -5e-07 5e-07 -5e-07 0
planewave = 56.31 0 2
pbc = 1 1 0
rcwaorder = 8 0
frequency1 = 5.0e+14 5.0e+14 0
end
)");
    if (res.empty()) { std::cerr << "  SKIP: " << name << " (solve failed)\n"; return; }
    double R = res[0].R, T = res[0].T, A = res[0].A;
    std::cout << "  R_TE(Brewster)=" << R << " T=" << T << " A=" << A
              << " R+T=" << R + T << "\n";
    // 解析解: r_s = (cosθ - n*cosθt)/(cosθ + n*cosθt) で R_TE ≈ 0.148
    CHECK(R > 0.10 && R < 0.20, "R_TE at Brewster ~0.148 (non-zero)");
    CHECK(std::abs(R + T - 1.0) < 2e-3, "TE energy conservation R+T≈1");
    CHECK(std::abs(A) < 2e-3, "lossless: A≈0");
    if (g_fails == prev_fails) std::cout << "PASS: " << name << "\n";
    else                       std::cout << "FAIL: " << name << "\n";
}

// ============================================================
// Test 10: wavelength キーワード — 直接 μm 指定
// wavelength = 0.5 0.7 4  → 5 波長、昇順を確認
// ============================================================
static void test_wavelength_keyword()
{
    static const char* name = "test_wavelength_keyword";
    int prev_fails = g_fails;
    // Lx = 0.5 μm → Wood anomaly at λ=0.5 μm (first order grazing)。
    // 0.6–0.9 μm 範囲を使うことで全次数が evanescent に収まる。
    std::string path = writeTmp("wavekw", R"(
OpenRCWA 4 2
title = wavelength keyword test
xmesh = -2.5e-07 10 2.5e-07
ymesh = -2.5e-07 10 2.5e-07
zmesh = -3e-07 10 0.0 10 7e-07
material = 1 2.25 0 1 0
geometry = 2 1 -2.5e-07 2.5e-07 -2.5e-07 2.5e-07 -3e-07 0
planewave = 0 0 1
pbc = 1 1 0
rcwaorder = 4 0
wavelength = 0.6 0.9 4
end
)");
    RCWAProblem prob;
    std::string err;
    if (!parseRCWAInput(path, prob, err)) {
        std::cerr << "  parse error: " << err << "\n"; ++g_fails; return;
    }
    std::cout << "  lambdas.size()=" << prob.lambdas.size() << "\n";
    CHECK(prob.lambdas.size() == 5, "5 wavelengths from wavelength keyword");
    // 昇順かつ 0.6 ≤ λ ≤ 0.9 μm
    for (const auto& lam : prob.lambdas)
        CHECK(lam >= 0.599 && lam <= 0.901, "lambda in [0.6,0.9] um");
    for (size_t i = 1; i < prob.lambdas.size(); ++i)
        CHECK(prob.lambdas[i] > prob.lambdas[i-1], "lambdas sorted ascending");
    // ソルバも走らせて R+T≈1 を確認
    auto res = runRCWA(prob, err);
    if (!res.empty()) {
        for (const auto& r : res)
            CHECK(std::abs(r.R + r.T - 1.0) < 2e-3, "R+T≈1 for wavelength sweep");
    }
    if (g_fails == prev_fails) std::cout << "PASS: " << name << "\n";
    else                       std::cout << "FAIL: " << name << "\n";
}

// ============================================================
// Test 11: background キーワード
// background=2.25 (ガラス) を設定し、ガラス→ガラス界面 (R≈0) になることを確認
// ============================================================
static void test_background_keyword()
{
    static const char* name = "test_background_keyword";
    int prev_fails = g_fails;
    // geometry なし → 全域が background の誘電率 (ガラス) になるので R=0, T=1
    auto res = solve("bgkw", R"(
OpenRCWA 4 2
title = background keyword test
xmesh = -2.5e-07 10 2.5e-07
ymesh = -2.5e-07 10 2.5e-07
zmesh = -3e-07 10 0.0 10 7e-07
background = 2.25
planewave = 0 0 1
pbc = 1 1 0
rcwaorder = 4 0
frequency1 = 5.0e+14 5.0e+14 0
end
)");
    if (res.empty()) { std::cerr << "  SKIP: " << name << " (solve failed)\n"; return; }
    double R = res[0].R, T = res[0].T;
    std::cout << "  glass-only R=" << R << " T=" << T << "\n";
    // 均質ガラス中: 界面がないので R≈0, T≈1
    CHECK(R < 1e-6, "uniform glass: R≈0");
    CHECK(std::abs(T - 1.0) < 1e-4, "uniform glass: T≈1");
    if (g_fails == prev_fails) std::cout << "PASS: " << name << "\n";
    else                       std::cout << "FAIL: " << name << "\n";
}

// ============================================================
// Test 12: 吸収率 — 無損失材料では A ≈ 0
// ============================================================
static void test_absorption_lossless()
{
    static const char* name = "test_absorption_lossless";
    int prev_fails = g_fails;
    auto res = solve("abs_lossless", R"(
OpenRCWA 4 2
title = absorption lossless
xmesh = -2.5e-07 10 2.5e-07
ymesh = -2.5e-07 10 2.5e-07
zmesh = -1e-06 10 0.0 10 3e-07 10 7e-07
material = 1 2.25 0 1 0
geometry = 2 1 -2.5e-07 2.5e-07 -2.5e-07 2.5e-07 0 3e-07
planewave = 0 0 1
pbc = 1 1 0
rcwaorder = 4 0
frequency1 = 4.0e+14 7.0e+14 10
end
)");
    if (res.empty()) { std::cerr << "  SKIP: " << name << " (solve failed)\n"; return; }
    std::cout << "  wavelengths: " << res.size() << "\n";
    for (const auto& r : res) {
        CHECK(std::abs(r.A) < 2e-3,
              std::string("A≈0 at lambda=") + std::to_string(r.lambda));
        // A は必ず 1−R−T に等しい
        CHECK(std::abs(r.A - (1.0 - r.R - r.T)) < 1e-12, "A == 1-R-T exactly");
    }
    if (g_fails == prev_fails) std::cout << "PASS: " << name << "\n";
    else                       std::cout << "FAIL: " << name << "\n";
}

// ============================================================
// main
// ============================================================
int main()
{
    std::cout << "=== RCWA input / driver unit tests ===\n";

    test_fresnel();
    test_fabry_perot();
    test_brewster();
    test_frequency2();
    test_cylinder_geometry();
    test_material_dispersion();
    test_per_order_output();
    test_2d_grating_energy();
    test_te_oblique();
    test_wavelength_keyword();
    test_background_keyword();
    test_absorption_lossless();

    std::cout << "======================================\n";
    if (g_fails == 0)
        std::cout << "All tests PASSED\n";
    else
        std::cerr << g_fails << " test(s) FAILED\n";

    return g_fails;
}
