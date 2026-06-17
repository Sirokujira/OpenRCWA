#pragma once
#include "rcwa/defns.h"
#include "rcwa/RCWAInput.h"
#include <string>
#include <vector>

// ============================================================
// RCWADriver : RCWAProblem (パース済みの .orcwa) を RCWASolver で解き、
// 波長ごとの反射率 R・透過率 T を返す。
//
// 構造は「上端(最大 z)= 入射側半無限媒質, 中間 = 有限層, 下端(最小 z)=
// 透過側半無限媒質」として層状に分解される。上端・下端のスラブは均質である
// 前提で半無限媒質として扱われ、その厚さは無視される (標準的な RCWA 設定)。
// ============================================================

struct RCWAResult
{
    scalar lambda;  // 波長 [μm]
    scalar R;       // 反射率 (0 次および全回折次数の和)
    scalar T;       // 透過率
};

// prob を解き、波長ごとの R/T を返す。失敗時は空配列を返し err を設定する。
std::vector<RCWAResult> runRCWA(const RCWAProblem& prob, std::string& err);
