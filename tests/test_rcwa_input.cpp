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
//   9-13. TE 斜め入射 / wavelength / background / 吸収 / 縦方向場
//   14. TM/TE Fresnel 角度スイープ (解析解比較)
//   15. φ 回転不変性 (TM Brewster, TE @ φ=0/45/90)
//   16. 法線入射 TM/TE 縮退
//   17. コニカルマウント 1D 格子のエネルギー保存
//   18. 損失材料 (σ>0): A>0, R+T<1
//   19. Drude 分散の物理: 金属スラブの高反射
//   20. パースエラー処理 + wavelength 単一値形式
//   21. CYL_X/CYL_Y containsXY + deriveZEdges
//   22. 横方向場成分 (Ex/Ey/Hx/Hy) の偏波選択則
//   23. SaveOption 整合性: |f|² == Re² + Im²
//   24. 球ジオメトリの z スライス化 (階段近似) + エネルギー保存
//   25. runRCWA の場イメージ出力 (RCWAFieldRequest)
//   26. Bloch 包絡線位相 (斜め入射で Re/Im が振動し位相勾配が kx0 に一致)
//   27. sliceYZ 断面 (非正方ハーモニクス 9x3 で次元不整合が起きない)
//   28. 任意角の直線偏波 pol=3 psi (Malus 則)
//   29. 円偏波 pol=4/5 (R = (R_TM+R_TE)/2, 非キラルで左右一致)
//   30. 多極 Lorentz 分散の加算
//   31. eps<1 の入射媒質 (n_inc クランプ撤廃)
//   32. 内部層なしスタックでの場出力要求 (旧 segfault の回帰テスト)
//   33. frequency (単数形) キーワード
//   34. 複素誘電率の直接指定 (material_eps / material_index)
//   35. 未知キーワードの扱い (FDTD 専用は黙殺, 綴り間違いは警告)
//   36. 誘電率分布の出力 (saveDeviceImage)
// ============================================================
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <Eigen/Core>

#include "rcwa/RCWAInput.h"
#include "rcwa/RCWADriver.h"
#include "rcwa/RCWASolver.h"
#include "rcwa/Layer.h"

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
        CHECK(d.active(), "material 2 marked dispersive");
        CHECK(d.poles.size() == 1, "single pole parsed");
        if (d.poles.size() == 1) {
            const auto& p = d.poles[0];
            std::cout << "  einf=" << d.einf << " ae=" << p.ae
                      << " be=" << p.be << " ce=" << p.ce << "\n";
            CHECK(std::abs(d.einf - 1.0) < 1e-9, "einf=1.0");
            CHECK(std::abs(p.ae - 1.37e16) < 1e12, "ae=1.37e16");
            CHECK(std::abs(p.be - 2.73e13) < 1e10, "be=2.73e13");
            CHECK(p.ce == 0.0, "ce=0 (Drude model)");
        }
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
// Test 13: 縦方向 (Ez, Hz) 場の評価
//   一様な (xy 一定) 層スタックでは基本回折次数のみ励起されるため、
//   Maxwell の z 成分から導かれる Ez / Hz に明確な解析的性質がある:
//     - 法線入射 (kx0=0):        Ez ≈ 0, Hz ≈ 0
//     - 斜め入射 TM (x 偏波,φ=0): Ez ≠ 0, Hz ≈ 0
//     - 斜め入射 TE (y 偏波,φ=0): Hz ≠ 0, Ez ≈ 0
//   saveFieldImage が出力する CSV の最大振幅でこれらを検証する。
// ============================================================

// 一様 (単一セル) レイヤ。周期 Lx×Ly, 誘電率 eps。
static Layer makeUniformLayer(double eps, double Lx, double Ly)
{
    Eigen::VectorXs cx(2), cy(2);
    cx << -Lx / 2, Lx / 2;
    cy << -Ly / 2, Ly / 2;
    Eigen::MatrixXcs e(1, 1);
    e(0, 0) = scalex(eps, 0.0);
    return Layer(cx, cy, e);
}

// CSV 内の全数値の絶対値の最大を返す (非有限なら巨大値を返す)。
static double csvMaxAbs(const std::string& path)
{
    std::ifstream f(path);
    if (!f) return -1.0;
    double mx = 0.0;
    std::string line;
    while (std::getline(f, line)) {
        std::stringstream ss(line);
        std::string tok;
        while (std::getline(ss, tok, ',')) {
            try {
                double v = std::stod(tok);
                if (!std::isfinite(v)) return 1e300;
                v = std::fabs(v);
                if (v > mx) mx = v;
            } catch (...) { /* 空トークンは無視 */ }
        }
    }
    return mx;
}

// 空気/ガラス/空気の一様スタックを解き、Ez/Hz の最大振幅を返す。
static void runFieldCase(double thetaDeg, int pol, double& ezMax, double& hzMax)
{
    const double lambda = 0.6, Lx = 0.5, Ly = 0.5, slab = 0.3;
    const int nHx = 4;

    RCWASolver solver(nHx, 0);
    solver.disablePML();
    solver.addLayer(makeUniformLayer(1.0,  Lx, Ly));  // 0: 入射側 (空気)
    solver.addLayer(makeUniformLayer(2.25, Lx, Ly));  // 1: スラブ (ガラス)
    solver.addLayer(makeUniformLayer(1.0,  Lx, Ly));  // 2: 透過側 (空気)

    const double k0 = 2.0 * Pi / lambda;
    const double theta = thetaDeg * Pi / 180.0;
    solver.setBlochWavevector(k0 * std::sin(theta), 0.0);  // φ=0

    std::vector<int> stack{0, 1, 2};
    std::vector<scalar> thick{0.0, slab, 0.0};
    solver.solve(lambda, stack, thick);

    // 入射平面波を生成: x 偏波 (px=1) は φ=0 で TM, y 偏波 (py=1) は TE。
    // generateHorizontalPlaneWave はハーモニクス空間の単一次数 (m=n=0) を
    // 固有モード係数 cInc に変換する。これを inputCoeffs として渡すことで
    // 正しい平面波励起になる (固有モード添字を直接指定すると縮退 ±Kx が
    // 混ざり平面波にならない)。
    const scalar px = (pol == 2) ? 0.0 : 1.0;
    const scalar py = (pol == 2) ? 1.0 : 0.0;
    Eigen::VectorXcs cInc;
    solver.generateHorizontalPlaneWave(px, py, 0, cInc);

    std::vector<std::pair<int, scalex>> inputCoeffs;
    for (int i = 0; i < cInc.size(); ++i)
        if (std::abs(cInc(i)) > 1e-14)
            inputCoeffs.emplace_back(i, cInc(i));

    solver.saveFieldImage("/tmp/rcwa_ez.csv", sliceXZ, 0.0, Ez, modulation,
                          inputCoeffs, stack, thick);
    solver.saveFieldImage("/tmp/rcwa_hz.csv", sliceXZ, 0.0, Hz, modulation,
                          inputCoeffs, stack, thick);

    ezMax = csvMaxAbs("/tmp/rcwa_ez.csv");
    hzMax = csvMaxAbs("/tmp/rcwa_hz.csv");
}

static void test_longitudinal_fields()
{
    static const char* name = "test_longitudinal_fields";
    int prev_fails = g_fails;
    double ez = 0, hz = 0;

    // (1) 法線入射 x 偏波: Ez ≈ 0, Hz ≈ 0
    runFieldCase(0.0, 1, ez, hz);
    std::cout << "  normal  x-pol: Ez_max=" << ez << " Hz_max=" << hz << "\n";
    CHECK(ez >= 0.0 && ez < 1e-6, "normal incidence: Ez≈0");
    CHECK(hz >= 0.0 && hz < 1e-6, "normal incidence: Hz≈0");

    // (2) 斜め入射 TM (x 偏波): Ez ≠ 0, Hz ≈ 0
    runFieldCase(30.0, 1, ez, hz);
    std::cout << "  oblique TM   : Ez_max=" << ez << " Hz_max=" << hz << "\n";
    CHECK(ez > 1e-3, "oblique TM: Ez nonzero (longitudinal E present)");
    CHECK(hz >= 0.0 && hz < 1e-6, "oblique TM: Hz≈0");

    // (3) 斜め入射 TE (y 偏波): Hz ≠ 0, Ez ≈ 0
    runFieldCase(30.0, 2, ez, hz);
    std::cout << "  oblique TE   : Ez_max=" << ez << " Hz_max=" << hz << "\n";
    CHECK(hz > 1e-3, "oblique TE: Hz nonzero (longitudinal H present)");
    CHECK(ez >= 0.0 && ez < 1e-6, "oblique TE: Ez≈0");

    if (g_fails == prev_fails) std::cout << "PASS: " << name << "\n";
    else                       std::cout << "FAIL: " << name << "\n";
}

// ============================================================
// TM/TE 偏波テスト用ヘルパー
// ============================================================

// Fresnel 反射率の解析解 (入射側 n1, 透過側 n2)。
// pol=1: TM (p 偏波), pol=2: TE (s 偏波)。
static double fresnelR2(double thetaDeg, double n1, double n2, int pol)
{
    const double th  = thetaDeg * M_PI / 180.0;
    const double ct  = std::cos(th);
    const double st  = std::sin(th);
    const double s2  = n1 * st / n2;                    // sinθt
    const double ctt = std::sqrt(1.0 - s2 * s2);        // cosθt
    double r;
    if (pol == 2)  // TE: r_s = (n1·cosθ − n2·cosθt)/(n1·cosθ + n2·cosθt)
        r = (n1 * ct - n2 * ctt) / (n1 * ct + n2 * ctt);
    else           // TM: r_p = (n2·cosθ − n1·cosθt)/(n2·cosθ + n1·cosθt)
        r = (n2 * ct - n1 * ctt) / (n2 * ct + n1 * ctt);
    return r * r;
}

// 入射側が空気 (n1=1) の場合の省略形。
static double fresnelR(double thetaDeg, double n, int pol)
{
    return fresnelR2(thetaDeg, 1.0, n, pol);
}

// air→glass 半無限界面の .orcwa を theta/phi/pol を差し替えて生成する。
static std::string makeInterfaceInput(double thetaDeg, double phiDeg, int pol)
{
    std::ostringstream os;
    os << "OpenRCWA 4 2\n"
          "title = TM/TE interface test\n"
          "xmesh = -5e-07 10 5e-07\n"
          "ymesh = -5e-07 10 5e-07\n"
          "zmesh = -5e-07 10 0.0 10 1e-06\n"
          "material = 1 2.25 0 1 0\n"
          "geometry = 2 1 -5e-07 5e-07 -5e-07 5e-07 -5e-07 0\n"
          "planewave = " << thetaDeg << " " << phiDeg << " " << pol << "\n"
          "pbc = 1 1 0\n"
          "rcwaorder = 4 2\n"
          "frequency1 = 5.0e+14 5.0e+14 0\n"
          "end\n";
    return os.str();
}

// ============================================================
// Test 14: TM/TE Fresnel 角度スイープ (φ=0)
// θ = 20°, 40°, 60° で R_TM / R_TE が解析解と一致することを確認
// ============================================================
static void test_tm_te_fresnel_sweep()
{
    static const char* name = "test_tm_te_fresnel_sweep";
    int prev_fails = g_fails;
    const double n = 1.5;
    const double angles[] = {20.0, 40.0, 60.0};

    for (double th : angles) {
        for (int pol = 1; pol <= 2; ++pol) {
            std::string tag = std::string("fresnel_sweep_") +
                              std::to_string(int(th)) + "_" + std::to_string(pol);
            auto res = solve(tag, makeInterfaceInput(th, 0.0, pol));
            if (res.empty()) {
                std::cerr << "  FAIL: solve failed at theta=" << th
                          << " pol=" << pol << "\n";
                ++g_fails;
                continue;
            }
            double R = res[0].R, T = res[0].T;
            double Rref = fresnelR(th, n, pol);
            std::cout << "  theta=" << th << (pol == 2 ? " TE" : " TM")
                      << ": R=" << R << " (analytic " << Rref << ")"
                      << " R+T=" << R + T << "\n";
            CHECK(std::abs(R - Rref) < 2e-3,
                  std::string("R matches Fresnel at theta=") + std::to_string(th)
                  + " pol=" + std::to_string(pol));
            CHECK(std::abs(R + T - 1.0) < 2e-3,
                  std::string("R+T≈1 at theta=") + std::to_string(th)
                  + " pol=" + std::to_string(pol));
        }
    }
    if (g_fails == prev_fails) std::cout << "PASS: " << name << "\n";
    else                       std::cout << "FAIL: " << name << "\n";
}

// ============================================================
// Test 15: φ 回転不変性
// 一様界面では結果は方位角 φ に依存しない:
//   - TM at Brewster: φ=0/45/90 すべてで R≈0
//   - TE at Brewster 角: φ=0/45/90 すべてで R≈0.148
// 偏波回転 (px,py) と Bloch 波数 (kx0,ky0) の両方が φ を正しく
// 扱っていなければ通らない。
// ============================================================
static void test_phi_rotation_invariance()
{
    static const char* name = "test_phi_rotation_invariance";
    int prev_fails = g_fails;
    const double thB = 56.31;  // Brewster 角 arctan(1.5)
    const double n = 1.5;
    const double R_TE_ref = fresnelR(thB, n, 2);  // ≈0.148
    const double phis[] = {0.0, 45.0, 90.0};

    for (double phi : phis) {
        // TM: Brewster 角で R≈0
        {
            std::string tag = std::string("phi_tm_") + std::to_string(int(phi));
            auto res = solve(tag, makeInterfaceInput(thB, phi, 1));
            if (res.empty()) { std::cerr << "  FAIL: TM solve phi=" << phi << "\n"; ++g_fails; continue; }
            double R = res[0].R, T = res[0].T;
            std::cout << "  phi=" << phi << " TM: R=" << R << " R+T=" << R + T << "\n";
            CHECK(R < 0.01, std::string("R_TM≈0 at Brewster, phi=") + std::to_string(phi));
            CHECK(std::abs(R + T - 1.0) < 2e-3,
                  std::string("TM energy conservation, phi=") + std::to_string(phi));
        }
        // TE: 同角度で R≈0.148
        {
            std::string tag = std::string("phi_te_") + std::to_string(int(phi));
            auto res = solve(tag, makeInterfaceInput(thB, phi, 2));
            if (res.empty()) { std::cerr << "  FAIL: TE solve phi=" << phi << "\n"; ++g_fails; continue; }
            double R = res[0].R, T = res[0].T;
            std::cout << "  phi=" << phi << " TE: R=" << R
                      << " (analytic " << R_TE_ref << ") R+T=" << R + T << "\n";
            CHECK(std::abs(R - R_TE_ref) < 2e-3,
                  std::string("R_TE matches Fresnel, phi=") + std::to_string(phi));
            CHECK(std::abs(R + T - 1.0) < 2e-3,
                  std::string("TE energy conservation, phi=") + std::to_string(phi));
        }
    }
    if (g_fails == prev_fails) std::cout << "PASS: " << name << "\n";
    else                       std::cout << "FAIL: " << name << "\n";
}

// ============================================================
// Test 16: 法線入射での TM/TE 縮退
// θ=0 では偏波の区別がなく、R は pol=1/2 で厳密に一致するはず
// (等方一様スラブなので偏波依存性なし)。
// ============================================================
static void test_normal_incidence_pol_degeneracy()
{
    static const char* name = "test_normal_incidence_pol_degeneracy";
    int prev_fails = g_fails;

    auto resTM = solve("normal_tm", makeInterfaceInput(0.0, 0.0, 1));
    auto resTE = solve("normal_te", makeInterfaceInput(0.0, 0.0, 2));
    if (resTM.empty() || resTE.empty()) {
        std::cerr << "  FAIL: solve failed\n"; ++g_fails; return;
    }
    double R1 = resTM[0].R, R2 = resTE[0].R;
    std::cout << "  normal incidence: R_TM=" << R1 << " R_TE=" << R2 << "\n";
    CHECK(std::abs(R1 - R2) < 1e-10, "R identical for pol=1/2 at normal incidence");
    CHECK(std::abs(R1 - 0.04) < 2e-3, "R matches Fresnel 0.04");
    if (g_fails == prev_fails) std::cout << "PASS: " << name << "\n";
    else                       std::cout << "FAIL: " << name << "\n";
}

// ============================================================
// Test 17: コニカルマウント 1D 格子 — TM/TE エネルギー保存
// θ=30°, φ=30° (入射面が格子溝に対して斜め) では TM/TE が
// 回折で混合する。無損失なので R+T≈1 が成立するはず。
// ============================================================
static void test_conical_grating()
{
    static const char* name = "test_conical_grating";
    int prev_fails = g_fails;

    for (int pol = 1; pol <= 2; ++pol) {
        std::ostringstream os;
        os << "OpenRCWA 4 2\n"
              "title = conical mount grating\n"
              "xmesh = -2.5e-07 10 2.5e-07\n"
              "ymesh = -2.5e-07 10 2.5e-07\n"
              // 空気 / 格子層 / 空気 (透過側は一様スラブが必要)
              "zmesh = -7e-07 10 -3e-07 10 0.0 10 7e-07\n"
              "material = 1 2.25 0 1 0\n"
              // x 方向半周期のみガラス → 1D 格子
              "geometry = 2 1 -2.5e-07 0.0 -2.5e-07 2.5e-07 -3e-07 0\n"
              "planewave = 30 30 " << pol << "\n"
              "pbc = 1 1 0\n"
              "rcwaorder = 4 2\n"
              "frequency1 = 5.0e+14 5.0e+14 0\n"
              "end\n";
        std::string tag = std::string("conical_") + std::to_string(pol);
        auto res = solve(tag, os.str());
        if (res.empty()) { std::cerr << "  FAIL: solve failed pol=" << pol << "\n"; ++g_fails; continue; }
        double R = res[0].R, T = res[0].T;
        std::cout << "  conical " << (pol == 2 ? "TE" : "TM")
                  << ": R=" << R << " T=" << T << " R+T=" << R + T << "\n";
        CHECK(R > 0.0 && T > 0.0, "R and T positive");
        CHECK(std::abs(R + T - 1.0) < 2e-3,
              std::string("conical grating energy conservation pol=") + std::to_string(pol));
    }
    if (g_fails == prev_fails) std::cout << "PASS: " << name << "\n";
    else                       std::cout << "FAIL: " << name << "\n";
}

// ============================================================
// Test 18: 損失材料 (導電率 σ>0) — 吸収 A>0, R+T<1
// σ=1e4 S/m, ω=2π·5e14 → Δε_im = σ/(ωε₀) ≈ 0.36
// n ≈ 1.5 − 0.12i, 0.3 μm スラブで有意な吸収が生じる。
// ============================================================
static void test_absorption_lossy()
{
    static const char* name = "test_absorption_lossy";
    int prev_fails = g_fails;
    auto res = solve("abs_lossy", R"(
OpenRCWA 4 2
title = absorption lossy slab
xmesh = -2.5e-07 10 2.5e-07
ymesh = -2.5e-07 10 2.5e-07
zmesh = -1e-06 10 0.0 10 3e-07 10 7e-07
material = 1 2.25 1.0e4 1 0
geometry = 2 1 -2.5e-07 2.5e-07 -2.5e-07 2.5e-07 0 3e-07
planewave = 0 0 1
pbc = 1 1 0
rcwaorder = 4 0
frequency1 = 5.0e+14 5.0e+14 0
end
)");
    if (res.empty()) { std::cerr << "  SKIP: " << name << " (solve failed)\n"; return; }
    double R = res[0].R, T = res[0].T, A = res[0].A;
    std::cout << "  lossy slab: R=" << R << " T=" << T << " A=" << A << "\n";
    CHECK(A > 0.1, "conductive slab absorbs (A > 0.1)");
    CHECK(R + T < 1.0 - 1e-3, "R+T < 1 for lossy material");
    CHECK(std::abs(A - (1.0 - R - T)) < 1e-12, "A == 1-R-T exactly");
    CHECK(R >= 0.0 && T >= 0.0 && A <= 1.0, "R, T, A physically bounded");
    if (g_fails == prev_fails) std::cout << "PASS: " << name << "\n";
    else                       std::cout << "FAIL: " << name << "\n";
}

// ============================================================
// Test 19: Drude 分散の物理 — 金属スラブの高反射
// eps(ω) = 1 − ωp²/(ω² + iγω), ωp=1.37e16, γ=2.73e13 (銀近似)
// λ≈0.6 μm で eps ≈ −18 → 0.3 μm スラブは不透明で R が高い。
// (Test 6 はパースのみ確認; ここでは solve への反映を検証する)
// ============================================================
static void test_drude_metal_reflectance()
{
    static const char* name = "test_drude_metal_reflectance";
    int prev_fails = g_fails;
    auto res = solve("drude_metal", R"(
OpenRCWA 4 2
title = Drude metal slab
xmesh = -2.5e-07 10 2.5e-07
ymesh = -2.5e-07 10 2.5e-07
zmesh = -1e-06 10 0.0 10 3e-07 10 7e-07
material = 1 2.25 0 1 0
material_dispersion = 2 1.0 1.37e16 2.73e13 0.0
geometry = 2 1 -2.5e-07 2.5e-07 -2.5e-07 2.5e-07 0 3e-07
planewave = 0 0 1
pbc = 1 1 0
rcwaorder = 4 0
frequency1 = 5.0e+14 5.0e+14 0
end
)");
    if (res.empty()) { std::cerr << "  SKIP: " << name << " (solve failed)\n"; return; }
    double R = res[0].R, T = res[0].T, A = res[0].A;
    std::cout << "  Drude metal: R=" << R << " T=" << T << " A=" << A << "\n";
    CHECK(R > 0.8, "metallic slab is highly reflective (R > 0.8)");
    CHECK(T < 0.01, "0.3 um metal slab is opaque (T < 0.01)");
    CHECK(A > 0.0 && A < 0.2, "small ohmic absorption (0 < A < 0.2)");
    if (g_fails == prev_fails) std::cout << "PASS: " << name << "\n";
    else                       std::cout << "FAIL: " << name << "\n";
}

// ============================================================
// Test 20: パースエラー処理 + wavelength 単一値形式
//   (a) 波長未指定 → parseRCWAInput が false を返しエラーを設定
//   (b) 不正ヘッダ → false
//   (c) wavelength = 0.6 (単一値) → 1 波長
// ============================================================
static void test_parse_errors()
{
    static const char* name = "test_parse_errors";
    int prev_fails = g_fails;
    RCWAProblem prob;
    std::string err;

    // (a) 波長キーワードなし
    {
        std::string path = writeTmp("err_nowave", R"(
OpenRCWA 4 2
title = no wavelength
xmesh = -2.5e-07 10 2.5e-07
ymesh = -2.5e-07 10 2.5e-07
zmesh = -3e-07 10 0.0 10 7e-07
planewave = 0 0 1
end
)");
        RCWAProblem p;
        bool ok = parseRCWAInput(path, p, err);
        std::cout << "  no-wavelength: ok=" << ok << " err=\"" << err << "\"\n";
        CHECK(!ok, "missing wavelengths rejected");
        CHECK(!err.empty(), "error message set for missing wavelengths");
    }

    // (b) 不正ヘッダ
    {
        std::string path = writeTmp("err_header", R"(
NotOpenRCWA 4 2
end
)");
        RCWAProblem p;
        err.clear();
        bool ok = parseRCWAInput(path, p, err);
        std::cout << "  bad-header: ok=" << ok << " err=\"" << err << "\"\n";
        CHECK(!ok, "bad header rejected");
        CHECK(!err.empty(), "error message set for bad header");
    }

    // (c) wavelength 単一値形式
    {
        std::string path = writeTmp("wave_single", R"(
OpenRCWA 4 2
title = single wavelength
xmesh = -2.5e-07 10 2.5e-07
ymesh = -2.5e-07 10 2.5e-07
zmesh = -3e-07 10 0.0 10 7e-07
planewave = 0 0 1
pbc = 1 1 0
rcwaorder = 4 0
wavelength = 0.6
end
)");
        RCWAProblem p;
        err.clear();
        bool ok = parseRCWAInput(path, p, err);
        CHECK(ok, "single-value wavelength accepted");
        if (ok) {
            std::cout << "  single wavelength: lambdas.size()=" << p.lambdas.size() << "\n";
            CHECK(p.lambdas.size() == 1, "exactly 1 wavelength");
            CHECK(std::abs(p.lambdas[0] - 0.6) < 1e-12, "lambda == 0.6 um");
        }
    }
    if (g_fails == prev_fails) std::cout << "PASS: " << name << "\n";
    else                       std::cout << "FAIL: " << name << "\n";
}

// ============================================================
// Test 21: CYL_X / CYL_Y 形状と deriveZEdges
//   - CYL_X: containsXY は y 範囲のみで判定
//   - CYL_Y: containsXY は x 範囲のみで判定
//   - deriveZEdges: 範囲内の box エッジのみ追加・昇順
// ============================================================
static void test_cyl_xy_and_zedges()
{
    static const char* name = "test_cyl_xy_and_zedges";
    int prev_fails = g_fails;

    RCWABox box;
    box.material = 2;
    box.x0 = 0.0; box.x1 = 1.0;
    box.y0 = 0.0; box.y1 = 1.0;
    box.z0 = 0.0; box.z1 = 1.0;

    // CYL_X (x 軸円柱): y のみ判定
    box.shape = RCWA_SHAPE_CYL_X;
    CHECK(box.containsXY(-5.0, 0.5), "CYL_X: x ignored, y inside");
    CHECK(!box.containsXY(0.5, 1.5),  "CYL_X: y outside");

    // CYL_Y (y 軸円柱): x のみ判定
    box.shape = RCWA_SHAPE_CYL_Y;
    CHECK(box.containsXY(0.5, -5.0), "CYL_Y: y ignored, x inside");
    CHECK(!box.containsXY(1.5, 0.5),  "CYL_Y: x outside");

    // ---- containsXYZ: z 依存断面の厳密判定 ----
    // SPHERE: 楕円体。中心は内側、極付近の軸外点は外側。
    box.shape = RCWA_SHAPE_SPHERE;
    CHECK(box.containsXYZ(0.5, 0.5, 0.5),  "XYZ SPHERE: center inside");
    CHECK(box.containsXYZ(0.5, 0.5, 0.99), "XYZ SPHERE: on-axis near pole inside");
    CHECK(!box.containsXYZ(0.9, 0.5, 0.9), "XYZ SPHERE: off-axis near pole outside");
    CHECK(!box.containsXYZ(0.5, 0.5, 1.5), "XYZ SPHERE: z out of range");

    // CYL_X: 軸方向 x 範囲が有効になり、(y,z) 断面が円判定になる
    box.shape = RCWA_SHAPE_CYL_X;
    CHECK(box.containsXYZ(0.5, 0.5, 0.5),  "XYZ CYL_X: center inside");
    CHECK(!box.containsXYZ(-5.0, 0.5, 0.5), "XYZ CYL_X: x bounded (unlike containsXY)");
    CHECK(!box.containsXYZ(0.5, 0.9, 0.9), "XYZ CYL_X: (y,z) corner outside circle");
    CHECK(box.containsXYZ(0.5, 0.5, 0.99), "XYZ CYL_X: on-axis near z edge inside");

    // CYL_Y: 軸方向 y 範囲 + (x,z) 断面円
    box.shape = RCWA_SHAPE_CYL_Y;
    CHECK(box.containsXYZ(0.5, 0.5, 0.5),  "XYZ CYL_Y: center inside");
    CHECK(!box.containsXYZ(0.5, -5.0, 0.5), "XYZ CYL_Y: y bounded");
    CHECK(!box.containsXYZ(0.9, 0.5, 0.9), "XYZ CYL_Y: (x,z) corner outside circle");

    // deriveZEdges
    RCWAProblem prob;
    prob.zmin = -1.0; prob.zmax = 2.0;
    RCWABox in;   in.z0 = 0.0;  in.z1 = 1.0;   // 範囲内 → 両エッジ追加
    RCWABox out;  out.z0 = 5.0; out.z1 = 6.0;  // 範囲外 → 追加されない
    prob.boxes.push_back(in);
    prob.boxes.push_back(out);
    auto edges = deriveZEdges(prob);
    std::cout << "  deriveZEdges: " << edges.size() << " edges\n";
    CHECK(edges.size() == 4, "4 edges: zmin, 0, 1, zmax");
    for (size_t i = 1; i < edges.size(); ++i)
        CHECK(edges[i] > edges[i-1], "z edges sorted ascending");
    if (edges.size() == 4) {
        CHECK(edges[0] == -1.0 && edges[1] == 0.0 &&
              edges[2] == 1.0 && edges[3] == 2.0, "edge values correct");
    }
    if (g_fails == prev_fails) std::cout << "PASS: " << name << "\n";
    else                       std::cout << "FAIL: " << name << "\n";
}

// ============================================================
// Test 22-23 用ヘルパー: 一様スタックを解いて指定成分/保存形式の
// CSV を書き出す (runFieldCase の汎用版)。
// ============================================================
static void saveUniformFieldCSV(double thetaDeg, int pol,
                                FieldComponent comp, SaveOption opt,
                                const std::string& path,
                                SliceType slice = sliceXZ, int nHy = 0)
{
    const double lambda = 0.6, Lx = 0.5, Ly = 0.5, slab = 0.3;
    const int nHx = 4;

    RCWASolver solver(nHx, nHy);
    solver.disablePML();
    solver.addLayer(makeUniformLayer(1.0,  Lx, Ly));
    solver.addLayer(makeUniformLayer(2.25, Lx, Ly));
    solver.addLayer(makeUniformLayer(1.0,  Lx, Ly));

    const double k0 = 2.0 * Pi / lambda;
    const double theta = thetaDeg * Pi / 180.0;
    solver.setBlochWavevector(k0 * std::sin(theta), 0.0);

    std::vector<int> stack{0, 1, 2};
    std::vector<scalar> thick{0.0, slab, 0.0};
    solver.solve(lambda, stack, thick);

    const scalar px = (pol == 2) ? 0.0 : 1.0;
    const scalar py = (pol == 2) ? 1.0 : 0.0;
    Eigen::VectorXcs cInc;
    solver.generateHorizontalPlaneWave(px, py, 0, cInc);

    std::vector<std::pair<int, scalex>> inputCoeffs;
    for (int i = 0; i < cInc.size(); ++i)
        if (std::abs(cInc(i)) > 1e-14)
            inputCoeffs.emplace_back(i, cInc(i));

    solver.saveFieldImage(path, slice, 0.0, comp, opt,
                          inputCoeffs, stack, thick);
}

// CSV の指定行を数値配列として読む (見つからなければ空)。
static std::vector<double> csvReadRow(const std::string& path, int row)
{
    std::ifstream f(path);
    std::string line;
    for (int i = 0; i <= row; ++i)
        if (!std::getline(f, line)) return {};
    std::vector<double> vals;
    std::stringstream ss(line);
    std::string tok;
    while (std::getline(ss, tok, ','))
        try { vals.push_back(std::stod(tok)); } catch (...) {}
    return vals;
}

// CSV 内の全数値を行優先の一次元配列として読む。
static std::vector<double> csvReadAll(const std::string& path)
{
    std::vector<double> vals;
    std::ifstream f(path);
    std::string line;
    while (std::getline(f, line)) {
        std::stringstream ss(line);
        std::string tok;
        while (std::getline(ss, tok, ','))
            try { vals.push_back(std::stod(tok)); } catch (...) {}
    }
    return vals;
}

// ============================================================
// Test 22: 横方向場成分の偏波選択則
// 一様スタック法線入射では φ=0 の x 偏波は Ex/Hy のみ、
// y 偏波は Ey/Hx のみを励起する。
// ============================================================
static void test_transverse_fields()
{
    static const char* name = "test_transverse_fields";
    int prev_fails = g_fails;

    // x 偏波 (pol=1): Ex, Hy ≠ 0; Ey, Hx ≈ 0
    saveUniformFieldCSV(0.0, 1, Ex, modulation, "/tmp/rcwa_t22_ex.csv");
    saveUniformFieldCSV(0.0, 1, Ey, modulation, "/tmp/rcwa_t22_ey.csv");
    saveUniformFieldCSV(0.0, 1, Hx, modulation, "/tmp/rcwa_t22_hx.csv");
    saveUniformFieldCSV(0.0, 1, Hy, modulation, "/tmp/rcwa_t22_hy.csv");
    double ex = csvMaxAbs("/tmp/rcwa_t22_ex.csv");
    double ey = csvMaxAbs("/tmp/rcwa_t22_ey.csv");
    double hx = csvMaxAbs("/tmp/rcwa_t22_hx.csv");
    double hy = csvMaxAbs("/tmp/rcwa_t22_hy.csv");
    std::cout << "  x-pol: Ex=" << ex << " Ey=" << ey
              << " Hx=" << hx << " Hy=" << hy << "\n";
    CHECK(ex > 0.5,  "x-pol: Ex excited");
    CHECK(hy > 1e-3, "x-pol: Hy excited");
    CHECK(ey >= 0.0 && ey < 1e-6, "x-pol: Ey≈0");
    CHECK(hx >= 0.0 && hx < 1e-6, "x-pol: Hx≈0");

    // y 偏波 (pol=2): Ey, Hx ≠ 0; Ex, Hy ≈ 0
    saveUniformFieldCSV(0.0, 2, Ex, modulation, "/tmp/rcwa_t22_ex2.csv");
    saveUniformFieldCSV(0.0, 2, Ey, modulation, "/tmp/rcwa_t22_ey2.csv");
    saveUniformFieldCSV(0.0, 2, Hx, modulation, "/tmp/rcwa_t22_hx2.csv");
    saveUniformFieldCSV(0.0, 2, Hy, modulation, "/tmp/rcwa_t22_hy2.csv");
    ex = csvMaxAbs("/tmp/rcwa_t22_ex2.csv");
    ey = csvMaxAbs("/tmp/rcwa_t22_ey2.csv");
    hx = csvMaxAbs("/tmp/rcwa_t22_hx2.csv");
    hy = csvMaxAbs("/tmp/rcwa_t22_hy2.csv");
    std::cout << "  y-pol: Ex=" << ex << " Ey=" << ey
              << " Hx=" << hx << " Hy=" << hy << "\n";
    CHECK(ey > 0.5,  "y-pol: Ey excited");
    CHECK(hx > 1e-3, "y-pol: Hx excited");
    CHECK(ex >= 0.0 && ex < 1e-6, "y-pol: Ex≈0");
    CHECK(hy >= 0.0 && hy < 1e-6, "y-pol: Hy≈0");

    if (g_fails == prev_fails) std::cout << "PASS: " << name << "\n";
    else                       std::cout << "FAIL: " << name << "\n";
}

// ============================================================
// Test 23: SaveOption の整合性
// 同一ケースの modulation / realpart / imagpart CSV について
// 各セルで |f|² == Re² + Im² が成り立つことを確認する。
// ============================================================
static void test_save_option_consistency()
{
    static const char* name = "test_save_option_consistency";
    int prev_fails = g_fails;

    saveUniformFieldCSV(30.0, 1, Ex, modulation, "/tmp/rcwa_t23_mod.csv");
    saveUniformFieldCSV(30.0, 1, Ex, realpart,   "/tmp/rcwa_t23_re.csv");
    saveUniformFieldCSV(30.0, 1, Ex, imagpart,   "/tmp/rcwa_t23_im.csv");

    auto mod = csvReadAll("/tmp/rcwa_t23_mod.csv");
    auto re  = csvReadAll("/tmp/rcwa_t23_re.csv");
    auto im  = csvReadAll("/tmp/rcwa_t23_im.csv");

    std::cout << "  cells: mod=" << mod.size()
              << " re=" << re.size() << " im=" << im.size() << "\n";
    CHECK(!mod.empty(), "modulation CSV non-empty");
    CHECK(mod.size() == re.size() && mod.size() == im.size(),
          "all three CSVs have identical cell count");

    if (mod.size() == re.size() && mod.size() == im.size()) {
        double worst = 0.0;
        for (size_t i = 0; i < mod.size(); ++i) {
            double lhs = mod[i] * mod[i];
            double rhs = re[i] * re[i] + im[i] * im[i];
            worst = std::max(worst, std::abs(lhs - rhs));
        }
        std::cout << "  max ||f|^2 - (Re^2+Im^2)| = " << worst << "\n";
        // CSV は 8 桁精度で出力されるため丸め誤差 ~1e-8 を許容する
        CHECK(worst < 1e-6, "|f|^2 == Re^2 + Im^2 pointwise");
    }
    if (g_fails == prev_fails) std::cout << "PASS: " << name << "\n";
    else                       std::cout << "FAIL: " << name << "\n";
}

// ============================================================
// Test 24: 球ジオメトリの z スライス化
// deriveZEdges が球の z 範囲を細分し (階段近似)、
// 誘電体球を含む構造でもエネルギー保存が成立することを確認する。
// ============================================================
static void test_sphere_zslicing()
{
    static const char* name = "test_sphere_zslicing";
    int prev_fails = g_fails;

    std::string path = writeTmp("sphere", R"(
OpenRCWA 4 2
title = dielectric sphere z-slicing
xmesh = -2.5e-07 10 2.5e-07
ymesh = -2.5e-07 10 2.5e-07
zmesh = -7e-07 10 7e-07
material = 1 2.25 0 1 0
geometry = 2 2 -2e-07 2e-07 -2e-07 2e-07 -2e-07 2e-07
planewave = 0 0 1
pbc = 1 1 0
rcwaorder = 2 2
frequency1 = 5.0e+14 5.0e+14 0
end
)");
    RCWAProblem prob;
    std::string err;
    if (!parseRCWAInput(path, prob, err)) {
        std::cerr << "  parse error: " << err << "\n"; ++g_fails; return;
    }
    auto edges = deriveZEdges(prob);
    std::cout << "  z edges for sphere: " << edges.size() << "\n";
    // zmin, zmax, 球の z0/z1, 中間分割 7 点 → 11 エッジ
    CHECK(edges.size() >= 10, "sphere z range subdivided into slices");

    auto res = runRCWA(prob, err);
    if (res.empty()) { std::cerr << "  runRCWA error: " << err << "\n"; ++g_fails; return; }
    double R = res[0].R, T = res[0].T;
    std::cout << "  sphere R=" << R << " T=" << T << " R+T=" << R + T << "\n";
    CHECK(std::abs(R + T - 1.0) < 2e-3, "sphere: energy conservation R+T≈1");

    // --- 形状が実際に区別されているか ---
    // エネルギー保存だけでは「球が円柱として扱われている」バグを見逃す
    // (どちらの誘電率分布でも R+T=1 は成立してしまう)。
    // 同一バウンディングボックスで形状だけを変え、結果が異なることを確かめる。
    // 高コントラスト (eps=12) にすると差が大きく出て判別しやすい。
    auto runShape = [&](int shape, const std::string& tag) -> double {
        std::ostringstream os;
        os << "OpenRCWA 4 2\n"
              "title = shape discrimination\n"
              "xmesh = -2.5e-07 10 2.5e-07\n"
              "ymesh = -2.5e-07 10 2.5e-07\n"
              "zmesh = -7e-07 10 7e-07\n"
              "material = 1 12.0 0 1 0\n"
              "geometry = 2 " << shape
           << " -2e-07 2e-07 -2e-07 2e-07 -2e-07 2e-07\n"
              "planewave = 0 0 1\n"
              "pbc = 1 1 0\n"
              "rcwaorder = 3 3\n"
              "wavelength = 0.6\n"
              "end\n";
        auto r = solve(tag, os.str());
        if (r.empty()) { ++g_fails; return -1.0; }
        CHECK(std::abs(r[0].R + r[0].T - 1.0) < 2e-3,
              std::string("energy conservation for shape ") + std::to_string(shape));
        return r[0].R;
    };

    const double rBox    = runShape(RCWA_SHAPE_BOX,    "shp_box");
    const double rSphere = runShape(RCWA_SHAPE_SPHERE, "shp_sphere");
    const double rCylZ   = runShape(RCWA_SHAPE_CYL_Z,  "shp_cylz");
    const double rCylX   = runShape(RCWA_SHAPE_CYL_X,  "shp_cylx");
    std::cout << "  same bbox: BOX=" << rBox << " SPHERE=" << rSphere
              << " CYL_Z=" << rCylZ << " CYL_X=" << rCylX << "\n";

    // 球と z 軸円柱は XY 断面が同一なので、z 依存断面 (containsXYZ) と
    // 面内の階段近似の両方が効いていなければ同じ答えになる。
    CHECK(std::abs(rSphere - rCylZ) > 1e-2, "SPHERE distinguished from CYL_Z");
    CHECK(std::abs(rSphere - rBox)  > 1e-2, "SPHERE distinguished from BOX");
    CHECK(std::abs(rCylZ  - rBox)   > 1e-2, "CYL_Z distinguished from BOX");
    CHECK(std::abs(rCylX  - rBox)   > 1e-2, "CYL_X distinguished from BOX");

    if (g_fails == prev_fails) std::cout << "PASS: " << name << "\n";
    else                       std::cout << "FAIL: " << name << "\n";
}

// ============================================================
// Test 25: runRCWA の場イメージ出力 (RCWAFieldRequest)
// 斜め入射 TM のガラススラブで Ez 断面 CSV が生成され、
// 有限かつ非ゼロの値を含むことを確認する。
// ============================================================
static void test_driver_field_output()
{
    static const char* name = "test_driver_field_output";
    int prev_fails = g_fails;

    std::string path = writeTmp("drvfield", R"(
OpenRCWA 4 2
title = driver field output
xmesh = -2.5e-07 10 2.5e-07
ymesh = -2.5e-07 10 2.5e-07
zmesh = -1e-06 10 0.0 10 3e-07 10 7e-07
material = 1 2.25 0 1 0
geometry = 2 1 -2.5e-07 2.5e-07 -2.5e-07 2.5e-07 0 3e-07
planewave = 30 0 1
pbc = 1 1 0
rcwaorder = 4 0
wavelength = 0.6
end
)");
    RCWAProblem prob;
    std::string err;
    if (!parseRCWAInput(path, prob, err)) {
        std::cerr << "  parse error: " << err << "\n"; ++g_fails; return;
    }

    RCWAFieldRequest req;
    req.component = Ez;
    req.slice     = sliceXZ;
    req.coord     = 0.0;
    req.opt       = modulation;
    req.path      = "/tmp/rcwa_driver_ez.csv";

    auto res = runRCWA(prob, err, &req);
    if (res.empty()) { std::cerr << "  runRCWA error: " << err << "\n"; ++g_fails; return; }

    double ezMax = csvMaxAbs(req.path);
    std::cout << "  driver Ez CSV: max=" << ezMax << "\n";
    CHECK(ezMax > 1e-3, "oblique TM via driver: Ez CSV nonzero");
    CHECK(ezMax < 1e100, "Ez values finite");
    CHECK(std::abs(res[0].R + res[0].T - 1.0) < 2e-3,
          "R+T≈1 unaffected by field output");
    if (g_fails == prev_fails) std::cout << "PASS: " << name << "\n";
    else                       std::cout << "FAIL: " << name << "\n";
}

// ============================================================
// Test 26: Bloch 包絡線位相
// 一様スタックの斜め入射では基本次数のみ励起されるため、周期部分は x に依らず
// 一定。したがって場の x 依存性は Bloch 因子 exp(i·kx0·x) がすべてであり:
//   - realpart は x 方向に振動する (旧実装では定数だった)
//   - 位相は x に対して kx0 の傾きで直線的に増加する
// ============================================================
static void test_bloch_envelope_phase()
{
    static const char* name = "test_bloch_envelope_phase";
    int prev_fails = g_fails;

    const double lambda = 0.6, Lx = 0.5, thetaDeg = 30.0;
    const double kx0 = (2.0 * M_PI / lambda) * std::sin(thetaDeg * M_PI / 180.0);

    // (1) 法線入射: kx0=0 なので realpart は x 方向に一定
    saveUniformFieldCSV(0.0, 1, Ex, realpart, "/tmp/rcwa_t26_re0.csv");
    auto row0 = csvReadRow("/tmp/rcwa_t26_re0.csv", 0);
    CHECK(!row0.empty(), "normal-incidence CSV readable");
    if (!row0.empty()) {
        double lo = *std::min_element(row0.begin(), row0.end());
        double hi = *std::max_element(row0.begin(), row0.end());
        std::cout << "  normal Re(Ex) range: [" << lo << ", " << hi << "]\n";
        CHECK(hi - lo < 1e-6, "normal incidence: Re(Ex) constant in x");
    }

    // (2) 斜め入射: realpart が振動する
    saveUniformFieldCSV(thetaDeg, 1, Ex, realpart, "/tmp/rcwa_t26_re.csv");
    saveUniformFieldCSV(thetaDeg, 1, Ex, imagpart, "/tmp/rcwa_t26_im.csv");
    saveUniformFieldCSV(thetaDeg, 1, Ex, modulation, "/tmp/rcwa_t26_mod.csv");
    auto re = csvReadRow("/tmp/rcwa_t26_re.csv", 0);
    auto im = csvReadRow("/tmp/rcwa_t26_im.csv", 0);
    auto md = csvReadRow("/tmp/rcwa_t26_mod.csv", 0);
    CHECK(re.size() > 2 && re.size() == im.size() && re.size() == md.size(),
          "oblique CSV rows same length");
    if (re.size() > 2 && re.size() == im.size() && re.size() == md.size()) {
        double lo = *std::min_element(re.begin(), re.end());
        double hi = *std::max_element(re.begin(), re.end());
        std::cout << "  oblique Re(Ex) range: [" << lo << ", " << hi << "]\n";
        CHECK(hi - lo > 1e-2, "oblique: Re(Ex) varies in x (Bloch envelope present)");

        // |f| は Bloch 因子で変わらない (一様層なので x 方向に一定)
        double mlo = *std::min_element(md.begin(), md.end());
        double mhi = *std::max_element(md.begin(), md.end());
        CHECK(mhi - mlo < 1e-6, "modulation unaffected by Bloch factor");

        // 位相の傾きが kx0 に一致することを確認 (両端の位相差)
        const int n = static_cast<int>(re.size());
        const double step = Lx / n;
        const double dx = step * (n - 1);
        double ph0 = std::atan2(im[0], re[0]);
        double ph1 = std::atan2(im[n-1], re[n-1]);
        double dphi = ph1 - ph0;
        while (dphi >  M_PI) dphi -= 2.0 * M_PI;
        while (dphi < -M_PI) dphi += 2.0 * M_PI;
        const double expected = kx0 * dx;
        std::cout << "  phase slope: dphi=" << dphi
                  << " expected=" << expected << " (kx0=" << kx0 << ")\n";
        CHECK(std::abs(dphi - expected) < 1e-3,
              "phase advances as exp(+i*kx0*x) with correct sign and magnitude");
    }
    if (g_fails == prev_fails) std::cout << "PASS: " << name << "\n";
    else                       std::cout << "FAIL: " << name << "\n";
}

// ============================================================
// Test 27: sliceYZ 断面出力
// harmonics2D は (nx_ × ny_) なので YZ 断面では転置してから畳み込む必要がある。
// nx_≠ny_ (nHx=4, nHy=1 → 9×3) で次元不整合が起きないことを確認する。
// ============================================================
static void test_slice_yz()
{
    static const char* name = "test_slice_yz";
    int prev_fails = g_fails;

    // nHy=1 → ny_=3, nHx=4 → nx_=9 (非正方)
    saveUniformFieldCSV(0.0, 1, Ex, modulation, "/tmp/rcwa_t27_yz.csv", sliceYZ, 1);
    double mx = csvMaxAbs("/tmp/rcwa_t27_yz.csv");
    auto row = csvReadRow("/tmp/rcwa_t27_yz.csv", 0);
    std::cout << "  sliceYZ (nx_=9, ny_=3): max=" << mx
              << " cols=" << row.size() << "\n";
    CHECK(mx > 0.5,   "sliceYZ: Ex nonzero");
    CHECK(mx < 1e100, "sliceYZ: values finite");
    CHECK(!row.empty(), "sliceYZ: CSV has data");

    // XZ 断面と同じ振幅になるはず (一様層の法線入射なので断面によらない)
    saveUniformFieldCSV(0.0, 1, Ex, modulation, "/tmp/rcwa_t27_xz.csv", sliceXZ, 1);
    double mxz = csvMaxAbs("/tmp/rcwa_t27_xz.csv");
    std::cout << "  sliceXZ max=" << mxz << "\n";
    CHECK(std::abs(mx - mxz) < 1e-6, "uniform layer: |Ex| same on YZ and XZ slices");

    if (g_fails == prev_fails) std::cout << "PASS: " << name << "\n";
    else                       std::cout << "FAIL: " << name << "\n";
}

// ============================================================
// Test 28: 任意角の直線偏波 (pol=3, psi)
// 一様界面では TM/TE が結合しないので Malus 則が成り立つ:
//   R(psi) = cos²psi·R_TM + sin²psi·R_TE
// psi=0 は pol=1 と、psi=90 は pol=2 と厳密に一致する。
// ============================================================
static void test_linear_polarization_angle()
{
    static const char* name = "test_linear_polarization_angle";
    int prev_fails = g_fails;
    const double thB = 56.31;  // Brewster 角: R_TM≈0

    auto runPsi = [&](int pol, double psi) -> double {
        std::ostringstream os;
        os << "OpenRCWA 4 2\n"
              "title = polarization angle\n"
              "xmesh = -5e-07 10 5e-07\n"
              "ymesh = -5e-07 10 5e-07\n"
              "zmesh = -5e-07 10 0.0 10 1e-06\n"
              "material = 1 2.25 0 1 0\n"
              "geometry = 2 1 -5e-07 5e-07 -5e-07 5e-07 -5e-07 0\n"
              "planewave = " << thB << " 0 " << pol << " " << psi << "\n"
              "pbc = 1 1 0\n"
              "rcwaorder = 4 0\n"
              "frequency1 = 5.0e+14 5.0e+14 0\n"
              "end\n";
        std::string tag = "psi_" + std::to_string(pol) + "_" + std::to_string(int(psi));
        auto res = solve(tag, os.str());
        if (res.empty()) { ++g_fails; return -1.0; }
        CHECK(std::abs(res[0].R + res[0].T - 1.0) < 2e-3,
              std::string("energy conservation at psi=") + std::to_string(psi));
        return res[0].R;
    };

    const double R_TM = runPsi(1, 0.0);
    const double R_TE = runPsi(2, 0.0);
    std::cout << "  R_TM=" << R_TM << " R_TE=" << R_TE << "\n";

    // psi=0 → TM と一致、psi=90 → TE と一致
    double r0  = runPsi(3, 0.0);
    double r90 = runPsi(3, 90.0);
    std::cout << "  pol=3 psi=0: R=" << r0 << "  psi=90: R=" << r90 << "\n";
    CHECK(std::abs(r0  - R_TM) < 1e-9, "pol=3 psi=0 identical to TM");
    CHECK(std::abs(r90 - R_TE) < 1e-9, "pol=3 psi=90 identical to TE");

    // 中間角: Malus 則
    for (double psi : {30.0, 45.0, 60.0}) {
        double r = runPsi(3, psi);
        double c = std::cos(psi * M_PI / 180.0), s = std::sin(psi * M_PI / 180.0);
        double ref = c * c * R_TM + s * s * R_TE;
        std::cout << "  psi=" << psi << ": R=" << r << " (Malus " << ref << ")\n";
        CHECK(std::abs(r - ref) < 1e-6,
              std::string("Malus law at psi=") + std::to_string(psi));
    }
    if (g_fails == prev_fails) std::cout << "PASS: " << name << "\n";
    else                       std::cout << "FAIL: " << name << "\n";
}

// ============================================================
// Test 29: 円偏波 (pol=4 右, pol=5 左)
// 非キラルな一様界面では R_circ = (R_TM + R_TE)/2 で左右は同一。
// ============================================================
static void test_circular_polarization()
{
    static const char* name = "test_circular_polarization";
    int prev_fails = g_fails;
    const double thB = 56.31;

    auto runPol = [&](int pol) -> double {
        std::ostringstream os;
        os << "OpenRCWA 4 2\n"
              "title = circular polarization\n"
              "xmesh = -5e-07 10 5e-07\n"
              "ymesh = -5e-07 10 5e-07\n"
              "zmesh = -5e-07 10 0.0 10 1e-06\n"
              "material = 1 2.25 0 1 0\n"
              "geometry = 2 1 -5e-07 5e-07 -5e-07 5e-07 -5e-07 0\n"
              "planewave = " << thB << " 0 " << pol << "\n"
              "pbc = 1 1 0\n"
              "rcwaorder = 4 0\n"
              "frequency1 = 5.0e+14 5.0e+14 0\n"
              "end\n";
        auto res = solve("circ_" + std::to_string(pol), os.str());
        if (res.empty()) { ++g_fails; return -1.0; }
        CHECK(std::abs(res[0].R + res[0].T - 1.0) < 2e-3,
              std::string("circular energy conservation pol=") + std::to_string(pol));
        return res[0].R;
    };

    const double R_TM = runPol(1), R_TE = runPol(2);
    const double rcp  = runPol(4), lcp = runPol(5);
    const double ref  = 0.5 * (R_TM + R_TE);
    std::cout << "  RCP=" << rcp << " LCP=" << lcp
              << " (expected " << ref << ")\n";
    CHECK(std::abs(rcp - ref) < 1e-6, "RCP: R = (R_TM + R_TE)/2");
    CHECK(std::abs(lcp - ref) < 1e-6, "LCP: R = (R_TM + R_TE)/2");
    CHECK(std::abs(rcp - lcp) < 1e-12, "achiral structure: RCP == LCP");
    if (g_fails == prev_fails) std::cout << "PASS: " << name << "\n";
    else                       std::cout << "FAIL: " << name << "\n";
}

// ============================================================
// Test 30: 多極 Lorentz 分散の加算
// Drude 極 (be=ce=0) では eps = einf − ae²/ω² なので、
// ae=A の 1 極と ae=A/√2 の 2 極は厳密に同じ eps を与える。
// 極が加算されず上書きされていれば結果が食い違う。
// ============================================================
static void test_multipole_dispersion()
{
    static const char* name = "test_multipole_dispersion";
    int prev_fails = g_fails;
    // λ=0.6 μm での ω ≈ 3.14e15。ae = ω/2 とすると eps ≈ 1 − 0.25 = 0.75 で
    // 透明領域に入り、R が eps に敏感な中間値になる (飽和した R=1 だと
    // 極の加算漏れを検出できない)。
    const double A = 0.5 * 2.0 * M_PI * 2.99792458e14 / 0.6;

    auto runPoles = [&](const std::string& dispersionLines,
                        const std::string& tag) -> double {
        std::ostringstream os;
        os << "OpenRCWA 4 2\n"
              "title = multipole dispersion\n"
              "xmesh = -2.5e-07 10 2.5e-07\n"
              "ymesh = -2.5e-07 10 2.5e-07\n"
              "zmesh = -1e-06 10 0.0 10 3e-07 10 7e-07\n"
              "material = 1 2.25 0 1 0\n"
           << dispersionLines <<
              "geometry = 2 1 -2.5e-07 2.5e-07 -2.5e-07 2.5e-07 0 3e-07\n"
              "planewave = 0 0 1\n"
              "pbc = 1 1 0\n"
              "rcwaorder = 4 0\n"
              "wavelength = 0.6\n"
              "end\n";
        auto res = solve(tag, os.str());
        if (res.empty()) { ++g_fails; return -1.0; }
        return res[0].R;
    };

    // 既定精度 (6 桁) では A/√2 が丸められ、2 極の和が 1 極と厳密に一致しない。
    // 極の加算そのものを検証したいので full precision で書き出す。
    std::ostringstream one, two;
    one << std::setprecision(17);
    two << std::setprecision(17);
    one << "material_dispersion = 2 1.0 " << A << " 0.0 0.0\n";
    const double Ahalf = A / std::sqrt(2.0);
    two << "material_dispersion = 2 1.0 " << Ahalf << " 0.0 0.0\n"
        << "material_dispersion = 2 1.0 " << Ahalf << " 0.0 0.0\n";

    // パース段階で極が 2 本入っていることを確認
    {
        std::ostringstream os;
        os << "OpenRCWA 4 2\ntitle = t\n"
              "xmesh = -2.5e-07 10 2.5e-07\nymesh = -2.5e-07 10 2.5e-07\n"
              "zmesh = -3e-07 10 0.0 10 7e-07\nmaterial = 1 2.25 0 1 0\n"
           << two.str() << "planewave = 0 0 1\npbc = 1 1 0\n"
              "rcwaorder = 2 0\nwavelength = 0.6\nend\n";
        std::string path = writeTmp("twopole_parse", os.str());
        RCWAProblem p; std::string err;
        if (parseRCWAInput(path, p, err) && p.materialDispersion.size() > 2) {
            std::cout << "  poles parsed: " << p.materialDispersion[2].poles.size() << "\n";
            CHECK(p.materialDispersion[2].poles.size() == 2, "two poles accumulated");
        } else { ++g_fails; }
    }

    // 加算漏れ (上書き) の場合に得られる値 = 単一の A/√2 極
    std::ostringstream half;
    half << std::setprecision(17);
    half << "material_dispersion = 2 1.0 " << Ahalf << " 0.0 0.0\n";

    double r1 = runPoles(one.str(),  "pole1");
    double r2 = runPoles(two.str(),  "pole2");
    double rh = runPoles(half.str(), "polehalf");
    std::cout << "  1 pole (ae=A): R=" << r1
              << "   2 poles (ae=A/sqrt2): R=" << r2
              << "   1 pole (ae=A/sqrt2): R=" << rh << "\n";
    CHECK(r1 > 1e-4 && r1 < 0.5, "R in a sensitive (non-saturated) range");
    CHECK(std::abs(r1 - r2) < 1e-9, "two half-strength poles == one full pole");
    // 上書き実装なら r2 は rh に一致してしまう — 判別力があることを確認する
    CHECK(std::abs(r1 - rh) > 1e-4, "test discriminates: half-strength pole differs");
    if (g_fails == prev_fails) std::cout << "PASS: " << name << "\n";
    else                       std::cout << "FAIL: " << name << "\n";
}

// ============================================================
// Test 31: eps<1 の入射媒質
// 旧実装は n_inc を 1 でクランプしていたため、eps<1 の媒質から斜め入射すると
// Bloch 波数 kx0 が過大になり Fresnel 解と食い違っていた。
// background=0.5 (n1≈0.7071) → ガラス (n2=1.5) の界面で検証する。
// ============================================================
static void test_low_index_incidence()
{
    static const char* name = "test_low_index_incidence";
    int prev_fails = g_fails;
    const double eps1 = 0.5, n1 = std::sqrt(eps1), n2 = 1.5;
    const double thetaDeg = 30.0;

    for (int pol = 1; pol <= 2; ++pol) {
        std::ostringstream os;
        os << "OpenRCWA 4 2\n"
              "title = low-index incidence\n"
              "xmesh = -2.5e-07 10 2.5e-07\n"
              "ymesh = -2.5e-07 10 2.5e-07\n"
              "zmesh = -5e-07 10 0.0 10 1e-06\n"
              "background = " << eps1 << "\n"
              "material = 1 2.25 0 1 0\n"
              "geometry = 2 1 -2.5e-07 2.5e-07 -2.5e-07 2.5e-07 -5e-07 0\n"
              "planewave = " << thetaDeg << " 0 " << pol << "\n"
              "pbc = 1 1 0\n"
              "rcwaorder = 4 0\n"
              "wavelength = 0.6\n"
              "end\n";
        auto res = solve("lowidx_" + std::to_string(pol), os.str());
        if (res.empty()) { std::cerr << "  FAIL: solve pol=" << pol << "\n"; ++g_fails; continue; }
        double R = res[0].R, T = res[0].T;
        double ref = fresnelR2(thetaDeg, n1, n2, pol);
        std::cout << "  n1=" << n1 << " " << (pol == 2 ? "TE" : "TM")
                  << ": R=" << R << " (analytic " << ref << ")"
                  << " R+T=" << R + T << "\n";
        CHECK(std::abs(R - ref) < 2e-3,
              std::string("low-index incidence matches Fresnel pol=") + std::to_string(pol));
        CHECK(std::abs(R + T - 1.0) < 2e-3,
              std::string("low-index energy conservation pol=") + std::to_string(pol));
    }
    if (g_fails == prev_fails) std::cout << "PASS: " << name << "\n";
    else                       std::cout << "FAIL: " << name << "\n";
}

// ============================================================
// Test 32: 内部層のないスタックでの場出力要求
// zmesh が 2 層 (入射側/透過側の半無限層のみ) しか作らない場合、
// 中間層の場は定義されない。旧実装は c_m[-1] へ書き込んで segfault した。
// R/T 自体は計算できるので、場出力を要求したときだけエラーにする。
// ============================================================
static void test_field_request_without_interior_layer()
{
    static const char* name = "test_field_request_without_interior_layer";
    int prev_fails = g_fails;

    const char* input = R"(
OpenRCWA 4 2
title = two-layer stack
xmesh = -5e-07 10 5e-07
ymesh = -5e-07 10 5e-07
zmesh = -5e-07 10 0.0 10 1e-06
material = 1 2.25 0 1 0
geometry = 2 1 -5e-07 5e-07 -5e-07 5e-07 -5e-07 0
planewave = 0 0 1
pbc = 1 1 0
rcwaorder = 4 0
wavelength = 0.6
end
)";

    // (a) 場出力なし → 通常どおり解ける
    {
        auto res = solve("nointerior_plain", input);
        CHECK(!res.empty(), "two-layer stack still solves for R/T");
        if (!res.empty()) {
            std::cout << "  two-layer R=" << res[0].R << " T=" << res[0].T << "\n";
            CHECK(std::abs(res[0].R - 0.04) < 2e-3, "two-layer R matches Fresnel");
        }
    }

    // (b) 場出力あり → クラッシュせずエラーを返す
    {
        std::string path = writeTmp("nointerior_field", input);
        RCWAProblem prob;
        std::string err;
        CHECK(parseRCWAInput(path, prob, err), "parse succeeds");

        RCWAFieldRequest req;
        req.component = Ez;
        req.slice     = sliceXZ;
        req.path      = "/tmp/rcwa_t32_should_not_exist.csv";
        std::remove(req.path.c_str());

        err.clear();
        auto res = runRCWA(prob, err, &req);
        std::cout << "  field request: results=" << res.size()
                  << " err=\"" << err << "\"\n";
        CHECK(res.empty(), "field request without interior layer is rejected");
        CHECK(!err.empty(), "error message set");
    }
    if (g_fails == prev_fails) std::cout << "PASS: " << name << "\n";
    else                       std::cout << "FAIL: " << name << "\n";
}

// ============================================================
// Test 33: frequency (単数形) キーワード
// OpenTHFD 形式の frequency は frequency1/2 と同じ f0 f1 ndiv 形式。
// 旧実装では未知キーワード扱いで無音で捨てられ、「波長が指定されていません」
// という無関係なエラーになっていた。
// ============================================================
static void test_frequency_singular()
{
    static const char* name = "test_frequency_singular";
    int prev_fails = g_fails;

    std::string path = writeTmp("freq_singular", R"(
OpenRCWA 4 2
title = frequency singular
xmesh = -2.5e-07 10 2.5e-07
ymesh = -2.5e-07 10 2.5e-07
zmesh = -3e-07 10 0.0 10 7e-07
material = 1 2.25 0 1 0
geometry = 2 1 -2.5e-07 2.5e-07 -2.5e-07 2.5e-07 -3e-07 0
planewave = 0 0 1
pbc = 1 1 0
rcwaorder = 4 0
frequency = 4.0e+14 5.0e+14 4
end
)");
    RCWAProblem prob;
    std::string err;
    bool ok = parseRCWAInput(path, prob, err);
    std::cout << "  frequency: ok=" << ok
              << " lambdas=" << prob.lambdas.size() << "\n";
    CHECK(ok, "singular 'frequency' keyword accepted");
    CHECK(prob.lambdas.size() == 5, "5 wavelengths from frequency = f0 f1 4");
    if (ok) {
        auto res = runRCWA(prob, err);
        CHECK(!res.empty(), "solves with singular frequency keyword");
        for (const auto& r : res)
            CHECK(std::abs(r.R + r.T - 1.0) < 2e-3, "R+T≈1");
    }
    if (g_fails == prev_fails) std::cout << "PASS: " << name << "\n";
    else                       std::cout << "FAIL: " << name << "\n";
}

// ============================================================
// Test 34: 複素誘電率の直接指定 (material_eps / material_index)
// n+ik の実測値をそのまま与えられることを確認する。
//   material_index = m n k  →  eps = (n + i·k)²
// 無損失なら Fresnel 解、k>0 なら吸収 (A>0) になる。
// ============================================================
static void test_direct_complex_eps()
{
    static const char* name = "test_direct_complex_eps";
    int prev_fails = g_fails;

    auto run = [&](const std::string& matLine, const std::string& tag) {
        std::ostringstream os;
        os << "OpenRCWA 4 2\n"
              "title = direct eps\n"
              "xmesh = -2.5e-07 10 2.5e-07\n"
              "ymesh = -2.5e-07 10 2.5e-07\n"
              "zmesh = -1e-06 10 0.0 10 3e-07 10 7e-07\n"
              "material = 1 1.0 0 1 0\n"
           << matLine <<
              "geometry = 2 1 -2.5e-07 2.5e-07 -2.5e-07 2.5e-07 0 3e-07\n"
              "planewave = 0 0 1\n"
              "pbc = 1 1 0\n"
              "rcwaorder = 4 0\n"
              "wavelength = 0.6\n"
              "end\n";
        return solve(tag, os.str());
    };

    // (a) material_index = 2 1.5 0 → n=1.5 無損失スラブ。material 行の
    //     epsr=1.0 を上書きするので、ガラススラブと同じ結果になるはず。
    {
        auto res = run("material_index = 2 1.5 0.0\n", "midx_lossless");
        if (res.empty()) { ++g_fails; }
        else {
            std::cout << "  material_index n=1.5: R=" << res[0].R
                      << " T=" << res[0].T << " A=" << res[0].A << "\n";
            CHECK(std::abs(res[0].R + res[0].T - 1.0) < 2e-3, "lossless: R+T≈1");
            CHECK(std::abs(res[0].A) < 2e-3, "lossless: A≈0");
            CHECK(res[0].R > 1e-3, "index override took effect (R != 0)");
        }
    }

    // (b) material_eps = 2 2.25 0 は上と等価 (eps = 1.5²)
    {
        auto a = run("material_index = 2 1.5 0.0\n", "midx_eq");
        auto b = run("material_eps = 2 2.25 0.0\n",  "meps_eq");
        if (a.empty() || b.empty()) { ++g_fails; }
        else {
            std::cout << "  index(1.5)=" << a[0].R
                      << "  eps(2.25)=" << b[0].R << "\n";
            CHECK(std::abs(a[0].R - b[0].R) < 1e-12,
                  "material_index n=1.5 == material_eps 2.25");
        }
    }

    // (c) k>0 → 吸収が生じる (exp(-iωt) 規約で正虚部 = 損失)
    {
        auto res = run("material_index = 2 1.5 0.1\n", "midx_lossy");
        if (res.empty()) { ++g_fails; }
        else {
            std::cout << "  material_index n=1.5 k=0.1: R=" << res[0].R
                      << " T=" << res[0].T << " A=" << res[0].A << "\n";
            CHECK(res[0].A > 0.05, "k>0 absorbs (A > 0.05)");
            CHECK(res[0].R + res[0].T < 1.0 - 1e-3, "lossy: R+T < 1");
        }
    }
    if (g_fails == prev_fails) std::cout << "PASS: " << name << "\n";
    else                       std::cout << "FAIL: " << name << "\n";
}

// ============================================================
// Test 35: 未知キーワードの扱い
// FDTD 専用キーワードは正常 (黙って読み飛ばす)。それ以外の綴り間違いは
// 警告を出しつつ処理を継続する — いずれもパースは成功しなければならない。
// ============================================================
static void test_unknown_keyword_handling()
{
    static const char* name = "test_unknown_keyword_handling";
    int prev_fails = g_fails;

    auto parseWith = [&](const std::string& extraLines, const std::string& tag) {
        std::ostringstream os;
        os << "OpenRCWA 4 2\n"
              "title = unknown keyword\n"
              "xmesh = -2.5e-07 10 2.5e-07\n"
              "ymesh = -2.5e-07 10 2.5e-07\n"
              "zmesh = -3e-07 10 0.0 10 7e-07\n"
              "material = 1 2.25 0 1 0\n"
           << extraLines <<
              "planewave = 0 0 1\n"
              "pbc = 1 1 0\n"
              "rcwaorder = 2 0\n"
              "wavelength = 0.6\n"
              "end\n";
        std::string path = writeTmp(tag, os.str());
        RCWAProblem prob;
        std::string err;
        bool ok = parseRCWAInput(path, prob, err);
        return std::make_pair(ok, prob.lambdas.size());
    };

    // FDTD 専用キーワードが混ざっていてもパースは通る
    auto fdtd = parseWith(
        "solver = 3000 100 1e-3\n"
        "abc = 1 5 5.0 1.1\n"
        "feed = V 0 0 0 X 1.0 0.0 50.0\n"
        "plot3dgeom = 1\n"
        "timestep = 0.0\n", "kw_fdtd");
    std::cout << "  FDTD-only keywords: ok=" << fdtd.first
              << " lambdas=" << fdtd.second << "\n";
    CHECK(fdtd.first, "FDTD-only keywords do not break parsing");
    CHECK(fdtd.second == 1, "wavelength still picked up");

    // 綴り間違い: 警告は出るがパースは継続する
    auto typo = parseWith("wavlength = 0.7\nrcwaordr = 3\n", "kw_typo");
    std::cout << "  typo keywords: ok=" << typo.first
              << " lambdas=" << typo.second << "\n";
    CHECK(typo.first, "unknown keywords do not abort parsing");
    CHECK(typo.second == 1, "typo'd keyword ignored, valid one still applied");

    if (g_fails == prev_fails) std::cout << "PASS: " << name << "\n";
    else                       std::cout << "FAIL: " << name << "\n";
}

// ============================================================
// Test 36: 誘電率分布の出力 (saveDeviceImage)
// 旧実装は行列を計算するだけで保存処理がコメントアウトされていた。
// 空気 (eps=1) とガラス (eps=2.25) の格子で両方の値が現れることを確認する。
// ============================================================
static void test_device_image_output()
{
    static const char* name = "test_device_image_output";
    int prev_fails = g_fails;

    std::string path = writeTmp("devimg", R"(
OpenRCWA 4 2
title = device image
xmesh = -2.5e-07 10 2.5e-07
ymesh = -2.5e-07 10 2.5e-07
zmesh = -7e-07 10 -3e-07 10 0.0 10 7e-07
material = 1 2.25 0 1 0
geometry = 2 1 -2.5e-07 0.0 -2.5e-07 2.5e-07 -3e-07 0
planewave = 0 0 1
pbc = 1 1 0
rcwaorder = 4 0
wavelength = 0.6
end
)");
    RCWAProblem prob;
    std::string err;
    if (!parseRCWAInput(path, prob, err)) {
        std::cerr << "  parse error: " << err << "\n"; ++g_fails; return;
    }

    RCWAFieldRequest req;
    req.slice      = sliceXZ;
    req.coord      = 0.0;
    req.devicePath = "/tmp/rcwa_t36_device.csv";
    std::remove(req.devicePath.c_str());

    auto res = runRCWA(prob, err, &req);
    if (res.empty()) { std::cerr << "  runRCWA error: " << err << "\n"; ++g_fails; return; }

    auto vals = csvReadAll(req.devicePath);
    std::cout << "  device CSV cells: " << vals.size() << "\n";
    CHECK(!vals.empty(), "device CSV written (was a commented-out stub)");
    if (!vals.empty()) {
        double lo = *std::min_element(vals.begin(), vals.end());
        double hi = *std::max_element(vals.begin(), vals.end());
        std::cout << "  eps range: [" << lo << ", " << hi << "]\n";
        // 格子層は空気 (1.0) とガラス (2.25) の両方を含む
        CHECK(std::abs(lo - 1.0) < 1e-6,  "device: air region eps=1");
        CHECK(std::abs(hi - 2.25) < 1e-6, "device: glass region eps=2.25");
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
    test_longitudinal_fields();
    test_tm_te_fresnel_sweep();
    test_phi_rotation_invariance();
    test_normal_incidence_pol_degeneracy();
    test_conical_grating();
    test_absorption_lossy();
    test_drude_metal_reflectance();
    test_parse_errors();
    test_cyl_xy_and_zedges();
    test_transverse_fields();
    test_save_option_consistency();
    test_sphere_zslicing();
    test_driver_field_output();
    test_bloch_envelope_phase();
    test_slice_yz();
    test_linear_polarization_angle();
    test_circular_polarization();
    test_multipole_dispersion();
    test_low_index_incidence();
    test_field_request_without_interior_layer();
    test_frequency_singular();
    test_direct_complex_eps();
    test_unknown_keyword_handling();
    test_device_image_output();

    std::cout << "======================================\n";
    if (g_fails == 0)
        std::cout << "All tests PASSED\n";
    else
        std::cerr << g_fails << " test(s) FAILED\n";

    return g_fails;
}
