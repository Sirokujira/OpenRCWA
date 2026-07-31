// ============================================================
// orcwa_rcwa : .orcwa 入力ファイルを RCWA ソルバで解く実行ファイル。
//
//   Usage : orcwa_rcwa [-o <out.csv>] <datafile.orcwa>
//
// 入力ファイル (FDTD ソルバと共通形式) を層状周期構造として解釈し、
// 波長ごとの反射率 R・透過率 T を標準出力および CSV に書き出す。
// ============================================================
#include "rcwa/RCWAInput.h"
#include "rcwa/RCWADriver.h"

#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>

static void usage()
{
    std::cout << "Usage : orcwa_rcwa [-o <out.csv>] [-v] [field options] <datafile.orcwa>\n"
              << "  -o <out.csv>  : output CSV file (default: <input>_rcwa.csv)\n"
              << "  -v            : verbose — also write per-order efficiencies\n"
              << "  断面出力オプション (最初の波長について CSV 出力):\n"
              << "  -field <c>    : Ex|Ey|Ez|Hx|Hy|Hz\n"
              << "  -device       : 誘電率分布 (実部) も出力する\n"
              << "  -slice <p>    : xy|xz|yz (default: xz)\n"
              << "  -scoord <v>   : slice coordinate [um] (default: 0)\n"
              << "  -sopt <o>     : mod|re|im (default: mod)\n"
              << "  -fieldout <f> : field CSV file (default: <input>_field.csv)\n"
              << "  -deviceout <f>: device CSV file (default: <input>_device.csv)\n";
}

// -field 引数 (Ex..Hz) を FieldComponent へ変換。不正なら false。
static bool parseFieldComponent(const std::string& s, FieldComponent& c)
{
    if      (s == "Ex") c = Ex;
    else if (s == "Ey") c = Ey;
    else if (s == "Ez") c = Ez;
    else if (s == "Hx") c = Hx;
    else if (s == "Hy") c = Hy;
    else if (s == "Hz") c = Hz;
    else return false;
    return true;
}

int main(int argc, char* argv[])
{
    std::string infile;
    std::string outfile;
    bool verbose = false;

    RCWAFieldRequest fieldReq;
    bool wantField = false, wantDevice = false;
    std::string fieldOut, deviceOut;

    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if (a == "-o" && i + 1 < argc) {
            outfile = argv[++i];
        } else if (a == "-v" || a == "--verbose") {
            verbose = true;
        } else if (a == "-field" && i + 1 < argc) {
            if (!parseFieldComponent(argv[++i], fieldReq.component)) {
                std::cerr << "*** 不正な -field 引数: " << argv[i] << "\n";
                return 1;
            }
            wantField = true;
        } else if (a == "-slice" && i + 1 < argc) {
            std::string p = argv[++i];
            if      (p == "xy") fieldReq.slice = sliceXY;
            else if (p == "xz") fieldReq.slice = sliceXZ;
            else if (p == "yz") fieldReq.slice = sliceYZ;
            else { std::cerr << "*** 不正な -slice 引数: " << p << "\n"; return 1; }
        } else if (a == "-scoord" && i + 1 < argc) {
            fieldReq.coord = std::stod(argv[++i]);
        } else if (a == "-sopt" && i + 1 < argc) {
            std::string o = argv[++i];
            if      (o == "mod") fieldReq.opt = modulation;
            else if (o == "re")  fieldReq.opt = realpart;
            else if (o == "im")  fieldReq.opt = imagpart;
            else { std::cerr << "*** 不正な -sopt 引数: " << o << "\n"; return 1; }
        } else if (a == "-fieldout" && i + 1 < argc) {
            fieldOut = argv[++i];
        } else if (a == "-device") {
            wantDevice = true;
        } else if (a == "-deviceout" && i + 1 < argc) {
            deviceOut = argv[++i];
            wantDevice = true;
        } else if (a == "--help" || a == "-h") {
            usage();
            return 0;
        } else {
            infile = a;
        }
    }

    if (infile.empty()) {
        usage();
        return 1;
    }
    if (outfile.empty()) {
        // 入力名から拡張子を除いて _rcwa.csv を付ける
        size_t dot = infile.find_last_of('.');
        outfile = (dot == std::string::npos ? infile : infile.substr(0, dot))
                + "_rcwa.csv";
    }

    RCWAProblem prob;
    std::string err;
    if (!parseRCWAInput(infile, prob, err)) {
        std::cerr << "*** 入力エラー: " << err << "\n";
        return 1;
    }

    std::cout << "<<< OpenRCWA (RCWA solver) >>>\n";
    std::cout << "Title    = " << prob.title << "\n";
    std::cout << "Period   = " << prob.Lx() << " x " << prob.Ly() << " [um]\n";
    std::cout << "Harmonics= " << (2 * prob.nHx + 1) << " x "
              << (2 * prob.nHy + 1) << "\n";
    std::cout << "Incidence= theta " << prob.theta << ", phi " << prob.phi
              << ", pol " << prob.pol << "\n";
    std::cout << "Layers   = " << (deriveZEdges(prob).size() - 1) << "\n";
    std::cout << "Wavelens = " << prob.lambdas.size()
              << " (" << prob.lambdas.back() << " - " << prob.lambdas.front()
              << " um)\n";
    std::cout << "=== solving ===\n";

    // 拡張子を除いた入力名 (既定の出力名に使う)
    auto stem = [&]() {
        size_t dot = infile.find_last_of('.');
        return dot == std::string::npos ? infile : infile.substr(0, dot);
    };
    if (wantField) {
        if (fieldOut.empty()) fieldOut = stem() + "_field.csv";
        fieldReq.path = fieldOut;
    }
    if (wantDevice) {
        if (deviceOut.empty()) deviceOut = stem() + "_device.csv";
        fieldReq.devicePath = deviceOut;
    }

    std::vector<RCWAResult> results =
        runRCWA(prob, err, fieldReq.wantsAnything() ? &fieldReq : nullptr);
    if (results.empty()) {
        std::cerr << "*** 計算エラー: " << err << "\n";
        return 1;
    }

    const int nHx = prob.nHx;
    const int nHy = prob.nHy;
    const int nOrd = (2*nHx+1) * (2*nHy+1);

    std::ofstream fout(outfile);
    fout << std::setprecision(8);
    if (verbose) {
        // per-order header: lambda, R, T, A, R+T, R_order_0, ..., T_order_N
        fout << "# lambda[um], R, T, A, R+T";
        for (int i = 0; i < nOrd; ++i) fout << ", R_order_" << i;
        for (int i = 0; i < nOrd; ++i) fout << ", T_order_" << i;
        fout << "\n";
    } else {
        fout << "# lambda[um], R, T, A, R+T\n";
    }

    std::cout << std::setprecision(6);
    for (const auto& r : results) {
        fout << r.lambda << ", " << r.R << ", " << r.T << ", " << r.A
             << ", " << (r.R + r.T);
        if (verbose) {
            for (scalar v : r.REF_orders) fout << ", " << v;
            for (scalar v : r.TRN_orders) fout << ", " << v;
        }
        fout << "\n";
        std::cout << "  lambda = " << r.lambda << " um   R = " << r.R
                  << "   T = " << r.T << "   A = " << r.A
                  << "   R+T = " << (r.R + r.T) << "\n";
    }
    fout.close();

    if (wantField)
        std::cout << "=== field  -> " << fieldReq.path << " ===\n";
    if (wantDevice)
        std::cout << "=== device -> " << fieldReq.devicePath << " ===\n";
    std::cout << "=== output -> " << outfile << " ===\n";
    std::cout << "=== normal end ===\n";
    return 0;
}
