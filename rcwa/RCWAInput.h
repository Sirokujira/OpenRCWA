#pragma once
#include "rcwa/defns.h"
#include <string>
#include <vector>

// ============================================================
// RCWAInput : OpenRCWA テキスト入力ファイル (.orcwa) を RCWA ソルバ向けの
// 問題定義に変換するパーサ。FDTD ソルバ (sol/input_data.c) と同じファイル
// 形式を共有しつつ、RCWA エンジン (rcwa/RCWASolver) を駆動するための層状
// 周期構造として解釈する。
//
// 座標・長さはパース時にメートル -> マイクロンへ正規化する (RCWA エンジンは
// スケール不変だが、慣用に合わせて μm を用いる)。波長も μm で保持する。
// ============================================================

// 直方体 (geometry shape=1) の領域。座標は μm。
struct RCWABox
{
    int    material;                 // 材料インデックス (materialEps の添字)
    scalar x0, x1, y0, y1, z0, z1;   // 直方体の範囲 [μm]
};

struct RCWAProblem
{
    std::string title;

    // 計算領域 / 周期 (μm)
    scalar xmin = 0, xmax = 0;
    scalar ymin = 0, ymax = 0;
    scalar zmin = 0, zmax = 0;

    scalar Lx() const { return xmax - xmin; }  // x 方向周期
    scalar Ly() const { return ymax - ymin; }  // y 方向周期

    // 高調波の片側次数 (全次数は 2*nH+1)
    int nHx = 10;
    int nHy = 0;

    // 材料インデックス -> 複素誘電率。0 = 空気, 1 = PEC, 2 以降がユーザ材料。
    std::vector<scalex> materialEps;

    // 物体形状 (直方体のみ対応)
    std::vector<RCWABox> boxes;

    // 平面波入射
    scalar theta = 0;   // 入射角 θ [deg]
    scalar phi   = 0;   // 入射角 φ [deg]
    int    pol   = 1;   // 偏波 (1: V, 2: H)

    // 背景 (物体に覆われない領域) の誘電率
    scalex backgroundEps = scalex(1.0, 0.0);

    // 計算する波長のリスト [μm] (frequency1 から生成)
    std::vector<scalar> lambdas;

    // 周期境界 (RCWA は本質的に周期境界。pbc が無い場合の既定も周期とする)
    bool periodicX = true;
    bool periodicY = true;
};

// .orcwa ファイルを解析して prob に格納する。成功時 true、失敗時 false を返し
// err にエラーメッセージを設定する。
bool parseRCWAInput(const std::string& path, RCWAProblem& prob, std::string& err);

// boxes と zmin/zmax から、z 方向の層境界 (昇順・重複なし) を導出する。
// 戻り値の隣り合う 2 値が 1 つの層 (z スラブ) を成す。
std::vector<scalar> deriveZEdges(const RCWAProblem& prob);

// 指定した z スラブ中央における誘電率の x 断面を、ブレークポイント座標 coordX
// (μm) とセルごとの誘電率 eps (サイズ coordX.size()-1) として返す (1D 格子向け)。
void sampleCrossSection1D(
    const RCWAProblem& prob,
    scalar zCenter,
    std::vector<scalar>& coordX,
    std::vector<scalex>& eps);
