#pragma once

#include "rcwa/defns.h"
#include "LayerSampler.h"
#include <functional>
#include <unordered_map>

class VarLayerSampler:
    public LayerSampler
{
private:
    Eigen::VectorXs x_;
    Eigen::VectorXs y_;
    Eigen::MatrixXcs eps_;

    std::unordered_map<int, std::function<scalar(scalar)>> xrules_;
    std::unordered_map<int, std::function<scalar(scalar)>> yrules_;
public: 
    VarLayerSampler(
        scalar lambda, 
        int nx, 
        int ny,
        const Eigen::VectorXs& x,
        const Eigen::VectorXs& y,
        const Eigen::MatrixXcs& eps,
        bool enablePML = true);
    
    void addXRule(int layout, std::function<scalar(scalar)> rule);
    void addYRule(int layout, std::function<scalar(scalar)> rule);

    void sample(scalar z, 
        Eigen::MatrixXcs& P,
        Eigen::MatrixXcs& Q) override;
    // GPU path not implemented for VarLayerSampler; inherits throwing default.
private:
    void assignSimulationRegion() override;
};