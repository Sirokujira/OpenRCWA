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
    scalar R;       // 反射率 (全回折次数の和)
    scalar T;       // 透過率 (全回折次数の和)
    scalar A;       // 吸収率 = 1 - R - T (損失ゼロなら ≈ 0)

    // 各回折次数の反射率・透過率 (長さ = (2*nHx+1)*(2*nHy+1))
    // 格子: REF_orders[ny*i+j] は (m=i-Nx, n=j-Ny) 次の反射効率
    std::vector<scalar> REF_orders;
    std::vector<scalar> TRN_orders;
};

// 断面出力の指定 (任意)。指定時は最初の波長についてのみ CSV を書き出す。
// path / devicePath はそれぞれ独立に指定でき、空なら出力しない。
struct RCWAFieldRequest
{
    FieldComponent component = Ez;         // 出力する場成分
    SliceType      slice     = sliceXZ;    // スライス面 (場・構造で共通)
    scalar         coord     = 0.0;        // スライス位置 [μm]
    SaveOption     opt       = modulation; // |f| / Re / Im
    std::string    path;                   // 場イメージの出力 CSV パス
    std::string    devicePath;             // 誘電率分布 (実部) の出力 CSV パス

    bool wantsAnything() const { return !path.empty() || !devicePath.empty(); }
};

// prob を解き、波長ごとの R/T を返す。失敗時は空配列を返し err を設定する。
// fieldReq を指定すると最初の波長で場イメージ CSV も出力する。
std::vector<RCWAResult> runRCWA(const RCWAProblem& prob, std::string& err,
                                const RCWAFieldRequest* fieldReq = nullptr);
