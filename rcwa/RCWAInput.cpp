#include "rcwa/RCWAInput.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <cstring>
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

// ============================================================
// RCWABox::containsXYZ — z 依存断面を含む厳密判定
// ============================================================
bool RCWABox::containsXYZ(scalar xc, scalar yc, scalar zc) const
{
    if (zc < z0 || zc > z1) return false;

    const scalar cx = 0.5*(x0+x1), cy = 0.5*(y0+y1), cz = 0.5*(z0+z1);
    const scalar rx = 0.5*(x1-x0), ry = 0.5*(y1-y0), rz = 0.5*(z1-z0);

    switch (shape) {
    case RCWA_SHAPE_BOX:
        return xc >= x0 && xc <= x1 && yc >= y0 && yc <= y1;

    case RCWA_SHAPE_SPHERE: {
        // 楕円体: 正規化距離の 2 乗和 ≤ 1
        if (rx <= 0 || ry <= 0 || rz <= 0) return false;
        scalar dx = (xc-cx)/rx, dy = (yc-cy)/ry, dz = (zc-cz)/rz;
        return dx*dx + dy*dy + dz*dz <= 1.0;
    }

    case RCWA_SHAPE_CYL_Z: {
        // z 軸円柱: XY 断面楕円 (z 依存なし)
        if (rx <= 0 || ry <= 0) return false;
        scalar dx = (xc-cx)/rx, dy = (yc-cy)/ry;
        return dx*dx + dy*dy <= 1.0;
    }

    case RCWA_SHAPE_CYL_X: {
        // x 軸円柱: 軸方向 x 範囲 + YZ 断面楕円
        if (ry <= 0 || rz <= 0) return false;
        if (xc < x0 || xc > x1) return false;
        scalar dy = (yc-cy)/ry, dz = (zc-cz)/rz;
        return dy*dy + dz*dz <= 1.0;
    }

    case RCWA_SHAPE_CYL_Y: {
        // y 軸円柱: 軸方向 y 範囲 + XZ 断面楕円
        if (rx <= 0 || rz <= 0) return false;
        if (yc < y0 || yc > y1) return false;
        scalar dx = (xc-cx)/rx, dz = (zc-cz)/rz;
        return dx*dx + dz*dz <= 1.0;
    }

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

// .orcwa は FDTD ソルバと共有する形式であり、時間領域固有のキーワードが
// 含まれていても正常。RCWA では意味を持たないので黙って読み飛ばす。
// (ここに無いキーワードは綴り間違いの可能性があるため警告する)
bool isFDTDOnlyKeyword(const std::string& k)
{
    // 出力・作図系はすべて接頭辞で判定する (plot*, far1d*, far2d*, near*)。
    // 個別に列挙すると FDTD 側の追加に追従できず誤警告の原因になる。
    static const char* kPrefixes[] = {
        "plot", "far1d", "far2d", "near1d", "near2d", "near3d",
    };
    for (const char* p : kPrefixes)
        if (k.compare(0, std::strlen(p), p) == 0) return true;

    // 時間領域固有の設定 (励振・境界・時間刻み・回路素子など)。
    static const std::set<std::string> kFDTDOnly = {
        "abc", "feed", "rfeed", "point", "load", "inductor", "source",
        "solver", "timestep", "pulsewidth", "freqdiv", "matchingloss",
        "name",
    };
    return kFDTDOnly.count(k) != 0;
}

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
    prob.materialEps.push_back(scalex(1.0, 1.0e8));    // 1: PEC (近似; exp(-iωt) 規約で損失は正虚部)
    prob.materialMu.clear();
    prob.materialMu.push_back(scalex(1.0, 0.0));  // 0: 空気
    prob.materialMu.push_back(scalex(1.0, 0.0));  // 1: PEC
    prob.materialSigma.clear();
    prob.materialSigma.push_back(0.0);  // 0: 空気
    prob.materialSigma.push_back(0.0);  // 1: PEC
    prob.materialMagSigma.clear();
    prob.materialMagSigma.push_back(0.0);  // 0: 空気
    prob.materialMagSigma.push_back(0.0);  // 1: PEC
    prob.materialDispersion.clear();
    prob.materialDispersion.emplace_back();  // 0: 空気 (非分散)
    prob.materialDispersion.emplace_back();  // 1: PEC  (非分散)

    bool haveFreq1 = false;
    bool haveFreq2 = false;
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
                    // amur = 比透磁率, msgm = 磁気導電率 [Ω/m]。
                    // 磁気導電率は損失として μ の虚部に入る (exp(-iωt) 規約で正)。
                    scalar amur = (nv >= 4) ? toScalar(V(3)) : 1.0;
                    scalar msgm = (nv >= 5) ? toScalar(V(4)) : 0.0;
                    prob.materialEps.push_back(scalex(epsr, 0.0));
                    prob.materialMu.push_back(scalex(amur, 0.0));
                    prob.materialSigma.push_back(esgm);
                    prob.materialMagSigma.push_back(msgm);
                    prob.materialDispersion.emplace_back();  // 非分散として初期化
                }
            }
            else if (key == "material_eps" || key == "material_index") {
                // RCWA 拡張: 材料 m の複素誘電率を直接指定する。
                //   material_eps   = m epsr [epsi]   誘電率そのもの
                //   material_index = m n    [k]      屈折率 (eps = (n + i·k)^2)
                // 時間規約 exp(-iωt) では損失は正の虚部 (k>0 が吸収)。
                if (nv >= 2) {
                    int m = toInt(V(0));
                    if (m >= 0 && m < static_cast<int>(prob.materialEps.size())) {
                        scalar a = toScalar(V(1));
                        scalar b = (nv >= 3) ? toScalar(V(2)) : 0.0;
                        prob.materialEps[m] = (key == "material_index")
                            ? scalex(a, b) * scalex(a, b)
                            : scalex(a, b);
                    } else {
                        std::cerr << "*** 警告: " << key << " の材料番号 " << m
                                  << " は未定義です (無視します)\n";
                    }
                }
            }
            else if (key == "material_mu") {
                // RCWA 拡張: 材料 m の複素比透磁率を直接指定する。
                //   material_mu = m mur [mui]
                if (nv >= 2) {
                    int m = toInt(V(0));
                    if (m >= 0 && m < static_cast<int>(prob.materialMu.size())) {
                        scalar a = toScalar(V(1));
                        scalar b = (nv >= 3) ? toScalar(V(2)) : 0.0;
                        prob.materialMu[m] = scalex(a, b);
                    } else {
                        std::cerr << "*** 警告: material_mu の材料番号 " << m
                                  << " は未定義です (無視します)\n";
                    }
                }
            }
            else if (key == "material_dispersion") {
                // material_dispersion = m einf ae be ce  [SI 単位: rad/s]
                // 同じ m に複数行を書くと極が加算される (多極 Lorentz モデル):
                //   eps(omega) = einf + Σ_p ae_p^2/(ce_p^2-omega^2-i*be_p*omega)
                // einf は最後に指定された値を採用する。
                if (nv >= 5) {
                    int m = toInt(V(0));
                    if (m >= 0 && m < static_cast<int>(prob.materialDispersion.size())) {
                        auto& d = prob.materialDispersion[m];
                        d.einf = toScalar(V(1));
                        RCWAProblem::LorentzPole p;
                        p.ae = toScalar(V(2));
                        p.be = toScalar(V(3));
                        p.ce = toScalar(V(4));
                        if (p.ae != 0.0) d.poles.push_back(p);
                    } else {
                        std::cerr << "*** 警告: material_dispersion の材料番号 "
                                  << m << " は未定義です (無視します)\n";
                    }
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
                    // 既知の形状のみ追加 (未対応形状は警告して無視)
                    if (shape == RCWA_SHAPE_BOX    || shape == RCWA_SHAPE_SPHERE  ||
                        shape == RCWA_SHAPE_CYL_Z  || shape == RCWA_SHAPE_CYL_X  ||
                        shape == RCWA_SHAPE_CYL_Y) {
                        prob.boxes.push_back(box);
                    } else {
                        std::cerr << "*** 警告: 未対応の geometry 形状 shape="
                                  << shape << " を無視します\n";
                    }
                }
            }
            else if (key == "planewave") {
                // planewave = theta phi pol [psi]
                if (nv >= 3) {
                    prob.theta = toScalar(V(0));
                    prob.phi   = toScalar(V(1));
                    prob.pol   = toInt(V(2));
                    if (nv >= 4) prob.psi = toScalar(V(3));
                    if (prob.pol < 1 || prob.pol > 5) {
                        std::cerr << "*** 警告: 未対応の偏波指定 pol=" << prob.pol
                                  << " (1=TM, 2=TE, 3=直線[psi], 4=右円, 5=左円)。"
                                     "TM として扱います\n";
                        prob.pol = 1;
                    }
                }
            }
            else if (key == "pbc") {
                if (nv >= 2) {
                    prob.periodicX = (toInt(V(0)) != 0);
                    prob.periodicY = (toInt(V(1)) != 0);
                    // RCWA は本質的に周期境界。非周期指定は無視されるため警告する。
                    if (!prob.periodicX || !prob.periodicY) {
                        std::cerr << "*** 警告: RCWA ソルバは周期境界のみ対応です"
                                     " (pbc=0 は無視され、周期境界として扱われます)\n";
                    }
                }
            }
            else if (key == "rcwaorder") {
                // rcwaorder = nHx [nHy]  (RCWA 固有: 高調波の片側次数)
                if (nv >= 1) prob.nHx = toInt(V(0));
                if (nv >= 2) prob.nHy = toInt(V(1));
            }
            else if (key == "frequency1" || key == "frequency2" ||
                     key == "frequency") {
                // frequency1/2 = f0 f1 ndiv  [Hz] -> 波長 [μm] に変換
                // 両方指定された場合は波長を結合し昇順ソート・重複除去する。
                // 単数形 "frequency" は OpenTHFD 形式で frequency1/2 を同時に
                // 設定するキーワード。RCWA では frequency1 と同義に扱う。
                bool isFirst = (key != "frequency2");
                bool& already = isFirst ? haveFreq1 : haveFreq2;
                if (!already && nv >= 3) {
                    scalar f0   = toScalar(V(0));
                    scalar f1   = toScalar(V(1));
                    int    ndiv = toInt(V(2));
                    if (ndiv < 0) ndiv = 0;
                    for (int i = 0; i <= ndiv; ++i) {
                        scalar f = (ndiv == 0) ? f0
                                 : f0 + (f1 - f0) * i / static_cast<scalar>(ndiv);
                        if (f > 0)
                            prob.lambdas.push_back((C0 / f) * M_TO_UM);
                    }
                    already = true;
                }
            }
            else if (key == "wavelength" || key == "lambda") {
                // wavelength = lam0 lam1 ndiv  [μm] (直接波長指定)
                // ndiv=0 のとき lam0 の 1 点のみ。
                if (nv >= 3) {
                    scalar lam0 = toScalar(V(0));
                    scalar lam1 = toScalar(V(1));
                    int    ndiv = toInt(V(2));
                    if (ndiv < 0) ndiv = 0;
                    for (int i = 0; i <= ndiv; ++i) {
                        scalar lam = (ndiv == 0) ? lam0
                                   : lam0 + (lam1 - lam0) * i / static_cast<scalar>(ndiv);
                        if (lam > 0)
                            prob.lambdas.push_back(lam);
                    }
                } else if (nv >= 1) {
                    // wavelength = lam  (単一波長)
                    scalar lam = toScalar(V(0));
                    if (lam > 0) prob.lambdas.push_back(lam);
                }
            }
            else if (key == "background") {
                // background = epsr [epsi]  — 物体に覆われない領域の誘電率
                if (nv >= 1) {
                    scalar epsr = toScalar(V(0));
                    scalar epsi = (nv >= 2) ? toScalar(V(1)) : 0.0;
                    prob.backgroundEps = scalex(epsr, epsi);
                }
            }
            else if (isFDTDOnlyKeyword(key)) {
                // 時間領域 (FDTD) 専用のキーワード。.orcwa は両ソルバで共有する
                // 形式なので、これらが書かれていても正常。黙って読み飛ばす。
            }
            else {
                // 上記のいずれでもない = 綴り間違いか RCWA 未対応の機能。
                // 黙って捨てると原因不明のエラーになるため必ず知らせる。
                std::cerr << "*** 警告: 未知のキーワード \"" << key
                          << "\" を無視します\n";
            }
        }
        catch (const std::exception& e) {
            err = std::string("解析エラー (") + key + "): " + e.what();
            return false;
        }
    }

    if (prob.lambdas.empty()) {
        err = "波長が指定されていません (frequency1 / frequency2 または wavelength キーワードを使用してください)";
        return false;
    }
    // frequency1 と frequency2 を結合した場合に重複を除去し昇順に整列する
    std::sort(prob.lambdas.begin(), prob.lambdas.end());
    prob.lambdas.erase(std::unique(prob.lambdas.begin(), prob.lambdas.end()),
                       prob.lambdas.end());
    if (prob.Lx() <= 0) {
        err = "x 方向の周期が不正です (xmesh を確認してください)";
        return false;
    }

    return true;
}

std::vector<scalar> deriveZEdges(const RCWAProblem& prob)
{
    // 断面が z に依存する曲面形状 (球, x/y 軸円柱) は、その z 範囲を
    // 細分して階段近似する。BOX / CYL_Z は z 依存がないので端点のみ。
    constexpr int N_CURVED_SLICES = 8;

    std::set<scalar> edges;
    edges.insert(prob.zmin);
    edges.insert(prob.zmax);
    for (const auto& b : prob.boxes) {
        if (b.z0 >= prob.zmin && b.z0 <= prob.zmax) edges.insert(b.z0);
        if (b.z1 >= prob.zmin && b.z1 <= prob.zmax) edges.insert(b.z1);

        if (b.shape == RCWA_SHAPE_SPHERE || b.shape == RCWA_SHAPE_CYL_X ||
            b.shape == RCWA_SHAPE_CYL_Y) {
            for (int i = 1; i < N_CURVED_SLICES; ++i) {
                scalar z = b.z0 + (b.z1 - b.z0) * i
                         / static_cast<scalar>(N_CURVED_SLICES);
                if (z >= prob.zmin && z <= prob.zmax) edges.insert(z);
            }
        }
    }
    return std::vector<scalar>(edges.begin(), edges.end());
}
