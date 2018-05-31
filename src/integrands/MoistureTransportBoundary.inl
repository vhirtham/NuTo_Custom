#pragma once

#include "MoistureTransportBoundary.h"
#include "MoistureTransport.h"

namespace NuTo
{

Integrands::MoistureTransportBoundary::MoistureTransportBoundary(
        const Integrands::MoistureTransport& moistureTransportIntegrand, DofType dofTypeWV, DofType dofTypeRH,
        double massExchangeRateWV, double massExchangeRateRH, double initialEnvironmentalRH)
    : mMoistureTransportIntegrand{moistureTransportIntegrand}
    , mDofTypeWV{dofTypeWV}
    , mDofTypeRH{dofTypeRH}
    , mExchangeRateWV{massExchangeRateWV}
    , mExchangeRateRH{massExchangeRateRH}
    , mEnvRH{initialEnvironmentalRH}
{
}

DofVector<double> Integrands::MoistureTransportBoundary::Gradient(const CellIpData& cellIpData, double deltaT)
{
    const Eigen::VectorXd WV = cellIpData.NodeValueVector(mDofTypeWV);
    const Eigen::MatrixXd Nw = cellIpData.N(mDofTypeWV);
    const Eigen::VectorXd RH = cellIpData.NodeValueVector(mDofTypeRH);
    const Eigen::MatrixXd Ng = cellIpData.N(mDofTypeRH);
    const double wv = (Nw * WV)[0];
    const double rh = (Ng * RH)[0];

    const double wveq = mMoistureTransportIntegrand.mWVEq.get()->value(0., mEnvRH, 0., 0.);

    DofVector<double> gradient;
    gradient[mDofTypeWV] = -mExchangeRateWV * (wv - wveq) * Nw.transpose();
    gradient[mDofTypeRH] = -mExchangeRateRH * (rh - mEnvRH) * Ng.transpose();

    return gradient;
}


DofMatrix<double> Integrands::MoistureTransportBoundary::Stiffness(const CellIpData& cellIpData, double deltaT)
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

void Integrands::MoistureTransportBoundary::SetEnvironmentalRelativeHumidity(double envRH)
{
    mEnvRH = envRH;
}
}
