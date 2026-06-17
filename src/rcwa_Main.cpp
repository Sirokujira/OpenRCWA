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
    std::cout << "Usage : orcwa_rcwa [-o <out.csv>] <datafile.orcwa>\n";
}

int main(int argc, char* argv[])
{
    std::string infile;
    std::string outfile;

    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if (a == "-o" && i + 1 < argc) {
            outfile = argv[++i];
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

    std::vector<RCWAResult> results = runRCWA(prob, err);
    if (results.empty()) {
        std::cerr << "*** 計算エラー: " << err << "\n";
        return 1;
    }

    std::ofstream fout(outfile);
    fout << "# lambda[um], R, T, R+T\n";
    fout << std::setprecision(8);
    std::cout << std::setprecision(6);
    for (const auto& r : results) {
        fout << r.lambda << ", " << r.R << ", " << r.T << ", "
             << (r.R + r.T) << "\n";
        std::cout << "  lambda = " << r.lambda << " um   R = " << r.R
                  << "   T = " << r.T << "   R+T = " << (r.R + r.T) << "\n";
    }
    fout.close();

    std::cout << "=== output -> " << outfile << " ===\n";
    std::cout << "=== normal end ===\n";
    return 0;
}
