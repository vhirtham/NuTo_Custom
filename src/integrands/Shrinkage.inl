#pragma once

#include "Shrinkage.h"

#include <cmath>

namespace NuTo
{
namespace Integrands
{

template <int TDim>
Shrinkage<TDim>::Shrinkage(DofType dofTypeDisp, DofType dofTypeWV, DofType dofTypeRH, double weightingFactor)
    : mDofDisp{dofTypeDisp}
    , mDofWV{dofTypeWV}
    , mDofRH{dofTypeRH}
    , mWeightingFactor{weightingFactor}
{
}


template <int TDim>
Eigen::VectorXd Shrinkage<TDim>::Stress(const CellIpData& cellIpData, double deltaT)
{
    const Eigen::VectorXd WV = cellIpData.NodeValueVector(mDofWV);
    const Eigen::MatrixXd Nw = cellIpData.N(mDofWV);
    const double wv = (Nw * WV)[0];

    const Eigen::VectorXd RH = cellIpData.NodeValueVector(mDofRH);
    const Eigen::MatrixXd Ng = cellIpData.N(mDofRH);
    const double rh = (Ng * RH)[0];

    return wv * CapillaryPressure(rh) * ComponentVector() * mWeightingFactor;
}

template <int TDim>
Eigen::VectorXd Shrinkage<TDim>::Strain(const CellIpData& cellIpData, double deltaT)
{
    const Eigen::VectorXd WV = cellIpData.NodeValueVector(mDofWV);
    const Eigen::MatrixXd Nw = cellIpData.N(mDofWV);
    const double wv = (Nw * WV)[0];

    const Eigen::VectorXd RH = cellIpData.NodeValueVector(mDofRH);
    const Eigen::MatrixXd Ng = cellIpData.N(mDofRH);
    const double rh = (Ng * RH)[0];

    double porosity = 0.2;
    double K = 15.e9 / (3 * porosity);
    double K_s = 30.e9 / (3 * porosity);

    double compliance = 1 / (3 * porosity) * (1 / K - 1 / K_s);

    return wv * compliance * CapillaryPressure(rh) * ComponentVector() * (mWeightingFactor - 1.);
}


template <int TDim>
double Shrinkage<TDim>::CapillaryPressure(double rh)
{
    double p_g = 0.;
    double p_0 = 0.;
    double rho_w = 1.;
    double R = 8.31448;
    double T = 293;
    double M_w = 28.0;

    return p_g - p_0 - rho_w * R * T * std::log(rh) / M_w * 1.e6;
}

template <>
inline Eigen::VectorXd Shrinkage<1>::ComponentVector()
{
    return Eigen::VectorXd::Ones(1.);
}

template <>
inline Eigen::VectorXd Shrinkage<2>::ComponentVector()
{
    return (Eigen::VectorXd(3) << 1., 1., 0.).finished();
}

template <>
inline Eigen::VectorXd Shrinkage<3>::ComponentVector()
{
    return (Eigen::VectorXd(6) << 1., 1., 1., 0., 0., 0.).finished();
}


} // namespace Integrands
} // namespace NuTo
