#pragma once

#include "MoistureTransportBoundary.h"

template <int TDim, typename TDCw, typename TDCg, typename TWVEq>
NuTo::DofVector<double>
NuTo::Integrands::MoistureTransportBoundary<TDim, TDCw, TDCg, TWVEq>::Gradient(const NuTo::CellIpData& cellIpData,
                                                                               double deltaT)
{
    const Eigen::VectorXd WV = cellIpData.NodeValueVector(mDofTypeWV);
    const Eigen::MatrixXd Nw = cellIpData.N(mDofTypeWV);
    const Eigen::VectorXd RH = cellIpData.NodeValueVector(mDofTypeRH);
    const Eigen::MatrixXd Ng = cellIpData.N(mDofTypeRH);
    const double wv = (Nw * WV)[0];
    const double rh = (Ng * RH)[0];

    const double wveq = TWVEq::value(0., 0.4, 0., 0.);
    DofVector<double> gradient;
    gradient[mDofTypeWV] = -mExchangeRateWV * (wv - wveq) * Nw.transpose();
    gradient[mDofTypeRH] = -mExchangeRateRH * (rh - 0.4) * Ng.transpose();

    return gradient;
}


template <int TDim, typename TDCw, typename TDCg, typename TWVEq>
NuTo::DofMatrix<double>
NuTo::Integrands::MoistureTransportBoundary<TDim, TDCw, TDCg, TWVEq>::Stiffness(const NuTo::CellIpData& cellIpData,
                                                                                double deltaT)
{
    const Eigen::VectorXd WV = cellIpData.NodeValueVector(mDofTypeWV);
    const Eigen::MatrixXd Nw = cellIpData.N(mDofTypeWV);
    const Eigen::VectorXd RH = cellIpData.NodeValueVector(mDofTypeRH);
    const Eigen::MatrixXd Ng = cellIpData.N(mDofTypeRH);


    DofMatrix<double> stiffness;
    stiffness(mDofTypeWV, mDofTypeWV) = -mExchangeRateWV * Nw.transpose() * Nw;
    stiffness(mDofTypeWV, mDofTypeRH) = 0. * Nw.transpose() * Ng;
    stiffness(mDofTypeRH, mDofTypeWV) = 0. * Ng.transpose() * Nw;
    stiffness(mDofTypeRH, mDofTypeRH) = -mExchangeRateRH * Ng.transpose() * Ng;

    return stiffness;
}
