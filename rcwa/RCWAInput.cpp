#include "rcwa/RCWAInput.h"

#include <fstream>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <set>

// ============================================================
// RCWABox::containsXY
// ============================================================
bool RCWABox::containsXY(scalar xc, scalar yc) const
{
    switch (shape) {
    case RCWA_SHAPE_BOX:
        return xc >= x0 && xc <= x1 && yc >= y0 && yc <= y1;

    case RCWA_SHAPE_SPHERE:
    case RCWA_SHAPE_CYL_Z: {
        // 楕円 / 球 (XY 断面): バウンディングボックス中心・半径で正規化
        scalar cx = 0.5*(x0+x1), cy = 0.5*(y0+y1);
        scalar rx = 0.5*(x1-x0), ry = 0.5*(y1-y0);
        if (rx <= 0 || ry <= 0) return false;
        scalar dx = (xc-cx)/rx, dy = (yc-cy)/ry;
        return dx*dx + dy*dy <= 1.0;
    }

    case RCWA_SHAPE_CYL_X:
        // x 軸方向の円柱: YZ 断面が楕円
        return yc >= y0 && yc <= y1;  // x は z スラブ判定済みなので y のみ

    case RCWA_SHAPE_CYL_Y:
        // y 軸方向の円柱: XZ 断面が楕円
        return xc >= x0 && xc <= x1;

    default:
        return false;
    }
}

namespace {

constexpr scalar C0      = 2.99792458e8;  // 光速 [m/s]
constexpr scalar M_TO_UM = 1.0e6;         // メートル -> マイクロン

// 1 行を空白で分割する
std::vector<std::string> tokenize(const std::string& line)
{
    std::vector<std::string> tokens;
    std::istringstream iss(line);
    std::string t;
    while (iss >> t) tokens.push_back(t);
    return tokens;
}

scalar toScalar(const std::string& s) { return std::stod(s); }
int    toInt(const std::string& s)    { return std::stoi(s); }

} // namespace

bool parseRCWAInput(const std::string& path, RCWAProblem& prob, std::string& err)
{
    std::ifstream fin(path);
    if (!fin) {
        err = "ファイルを開けません: " + path;
        return false;
    }

    // 材料インデックス 0 = 空気, 1 = PEC (プレースホルダ)
    prob.materialEps.clear();
    prob.materialEps.push_back(scalex(1.0, 0.0));      // 0: 空気
    prob.materialEps.push_back(scalex(1.0, -1.0e8));   // 1: PEC (近似)
    prob.materialSigma.clear();
    prob.materialSigma.push_back(0.0);  // 0: 空気
    prob.materialSigma.push_back(0.0);  // 1: PEC

    bool haveFreq = false;
    int  nline = 0;

    std::string line;
    while (std::getline(fin, line)) {
        // コメント / 空行をスキップ
        if (line.empty()) continue;

        std::vector<std::string> tok = tokenize(line);
        if (tok.empty()) continue;

        // ヘッダ行 (OpenRCWA <major> <minor>)
        if (nline == 0) {
            if (tok[0] != "OpenRCWA" && tok[0] != "OpenTHFD") {
                err = "OpenRCWA 形式ではありません";
                return false;
            }
            nline++;
            continue;
        }
        nline++;

        const std::string& key = tok[0];

        // "key = ..." 形式以外 (end など) はキーワード単独で処理
        if (key == "end") break;

        // title は "=" 以降をそのまま使う
        if (key == "title") {
            auto pos = line.find('=');
            if (pos != std::string::npos && pos + 2 <= line.size())
                prob.title = line.substr(pos + 1);
            // 先頭空白を除去
            size_t b = prob.title.find_first_not_of(" \t");
            if (b != std::string::npos) prob.title = prob.title.substr(b);
            continue;
        }

        // それ以外は "key = v1 v2 ..." を想定。tok[1] が "=" でなければスキップ
        if (tok.size() < 3 || tok[1] != "=") continue;
        // 値トークン (tok[2] 以降)
        const int nv = static_cast<int>(tok.size()) - 2;
        auto V = [&](int i) -> const std::string& { return tok[2 + i]; };

        try {
            if (key == "xmesh") {
                // xmesh = x0 d0 x1 d1 ... xn -> 両端のみ使用
                prob.xmin = toScalar(V(0)) * M_TO_UM;
                prob.xmax = toScalar(V(nv - 1)) * M_TO_UM;
            }
            else if (key == "ymesh") {
                prob.ymin = toScalar(V(0)) * M_TO_UM;
                prob.ymax = toScalar(V(nv - 1)) * M_TO_UM;
            }
            else if (key == "zmesh") {
                prob.zmin = toScalar(V(0)) * M_TO_UM;
                prob.zmax = toScalar(V(nv - 1)) * M_TO_UM;
            }
            else if (key == "material") {
                // material = type epsr esgm amur msgm
                // type==1: 通常材料。epsr を実部、esgm を導電率として保存する。
                if (nv >= 2) {
                    scalar epsr = toScalar(V(1));
                    scalar esgm = (nv >= 3) ? toScalar(V(2)) : 0.0;
                    prob.materialEps.push_back(scalex(epsr, 0.0));
                    prob.materialSigma.push_back(esgm);
                }
            }
            else if (key == "geometry") {
                // geometry = m shape x0 x1 y0 y1 z0 z1
                // shape=1: 直方体, shape=2: 球, shape=11: z 軸円柱, etc.
                if (nv >= 8) {
                    int m     = toInt(V(0));
                    int shape = toInt(V(1));
                    RCWABox box;
                    box.material = m;
                    box.shape = static_cast<RCWAShape>(shape);
                    box.x0 = toScalar(V(2)) * M_TO_UM;
                    box.x1 = toScalar(V(3)) * M_TO_UM;
                    box.y0 = toScalar(V(4)) * M_TO_UM;
                    box.y1 = toScalar(V(5)) * M_TO_UM;
                    box.z0 = toScalar(V(6)) * M_TO_UM;
                    box.z1 = toScalar(V(7)) * M_TO_UM;
                    // 既知の形状のみ追加 (未対応形状は無視)
                    if (shape == RCWA_SHAPE_BOX    || shape == RCWA_SHAPE_SPHERE  ||
                        shape == RCWA_SHAPE_CYL_Z  || shape == RCWA_SHAPE_CYL_X  ||
                        shape == RCWA_SHAPE_CYL_Y) {
                        prob.boxes.push_back(box);
                    }
                }
            }
            else if (key == "planewave") {
                if (nv >= 3) {
                    prob.theta = toScalar(V(0));
                    prob.phi   = toScalar(V(1));
                    prob.pol   = toInt(V(2));
                }
            }
            else if (key == "pbc") {
                if (nv >= 2) {
                    prob.periodicX = (toInt(V(0)) != 0);
                    prob.periodicY = (toInt(V(1)) != 0);
                }
            }
            else if (key == "rcwaorder") {
                // rcwaorder = nHx [nHy]  (RCWA 固有: 高調波の片側次数)
                if (nv >= 1) prob.nHx = toInt(V(0));
                if (nv >= 2) prob.nHy = toInt(V(1));
            }
            else if (key == "frequency1" || key == "frequency2") {
                // frequency1 = f0 f1 ndiv  [Hz] -> 波長 [μm] に変換
                if (!haveFreq && nv >= 3) {
                    scalar f0   = toScalar(V(0));
                    scalar f1   = toScalar(V(1));
                    int    ndiv = toInt(V(2));
                    if (ndiv < 0) ndiv = 0;
                    prob.lambdas.clear();
                    for (int i = 0; i <= ndiv; ++i) {
                        scalar f = (ndiv == 0) ? f0
                                 : f0 + (f1 - f0) * i / static_cast<scalar>(ndiv);
                        if (f > 0)
                            prob.lambdas.push_back((C0 / f) * M_TO_UM);
                    }
                    haveFreq = true;
                }
            }
            // その他のキーワード (solver, point, plot* 等) は RCWA では無視
        }
        catch (const std::exception& e) {
            err = std::string("解析エラー (") + key + "): " + e.what();
            return false;
        }
    }

    if (prob.lambdas.empty()) {
        err = "周波数 (frequency1) が指定されていません";
        return false;
    }
    if (prob.Lx() <= 0) {
        err = "x 方向の周期が不正です (xmesh を確認してください)";
        return false;
    }

    return true;
}

std::vector<scalar> deriveZEdges(const RCWAProblem& prob)
{
    std::set<scalar> edges;
    edges.insert(prob.zmin);
    edges.insert(prob.zmax);
    for (const auto& b : prob.boxes) {
        if (b.z0 >= prob.zmin && b.z0 <= prob.zmax) edges.insert(b.z0);
        if (b.z1 >= prob.zmin && b.z1 <= prob.zmax) edges.insert(b.z1);
    }
    return std::vector<scalar>(edges.begin(), edges.end());
}

void sampleCrossSection1D(
    const RCWAProblem& prob,
    scalar zCenter,
    std::vector<scalar>& coordX,
    std::vector<scalex>& eps)
{
    // この z スラブと交差する直方体の x エッジを集めてブレークポイントを作る
    std::set<scalar> xs;
    xs.insert(prob.xmin);
    xs.insert(prob.xmax);
    for (const auto& b : prob.boxes) {
        if (zCenter < b.z0 || zCenter > b.z1) continue;
        scalar bx0 = std::max(b.x0, prob.xmin);
        scalar bx1 = std::min(b.x1, prob.xmax);
        if (bx0 < bx1) {
            xs.insert(bx0);
            xs.insert(bx1);
        }
    }

    coordX.assign(xs.begin(), xs.end());

    // 各セル中央で物体の有無を判定し、誘電率を決める (後勝ち)
    scalar yc = 0.5 * (prob.ymin + prob.ymax);
    eps.clear();
    for (size_t i = 0; i + 1 < coordX.size(); ++i) {
        scalar xc = 0.5 * (coordX[i] + coordX[i + 1]);
        scalex e  = prob.backgroundEps;
        for (const auto& b : prob.boxes) {
            if (zCenter < b.z0 || zCenter > b.z1) continue;
            if (b.containsXY(xc, yc)) {
                if (b.material >= 0 &&
                    b.material < static_cast<int>(prob.materialEps.size())) {
                    e = prob.materialEps[b.material];
                }
            }
        }
        eps.push_back(e);
    }
}
