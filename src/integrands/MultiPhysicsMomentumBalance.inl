#pragma once

#include "integrands/MultiPhysicsMomentumBalance.h"

namespace NuTo
{
namespace Integrands
{


template <int TDim>
MultiPhysicsMomentumBalance<TDim>::MultiPhysicsMomentumBalance(DofType dofDisp, DofType dofWV, DofType dofRH,
                                                               const Laws::MechanicsInterface<TDim>& mechanicsLaw, Shrinkage<TDim> &shrinkageLaw)
    : mDofDisp{dofDisp}
    , mDofWV{dofWV}
    , mDofRH{dofRH}
    , mLaw{mechanicsLaw}
    , mShrinkageIntegrand{shrinkageLaw}
{
}

template <int TDim>
Eigen::VectorXd MultiPhysicsMomentumBalance<TDim>::Strain(const CellIpData& cellIpData, double deltaT)
{
    return cellIpData.Apply(mDofDisp, Nabla::Strain());
}


template <int TDim>
Eigen::VectorXd MultiPhysicsMomentumBalance<TDim>::StrainMechanics(const CellIpData& cellIpData, double deltaT)
{
    return Strain(cellIpData,deltaT) - StrainShrinkage(cellIpData,deltaT);
}


template <int TDim>
Eigen::VectorXd MultiPhysicsMomentumBalance<TDim>::StrainShrinkage(const CellIpData& cellIpData, double deltaT)
{
    return mShrinkageIntegrand.Strain(cellIpData, deltaT);
}

template <int TDim>
Eigen::VectorXd MultiPhysicsMomentumBalance<TDim>::Stress(const CellIpData& cellIpData, double deltaT)
{
    return StressMechanics(cellIpData,deltaT) + StressShrinkage(cellIpData,deltaT);
}

template <int TDim>
Eigen::VectorXd MultiPhysicsMomentumBalance<TDim>::StressMechanics(const CellIpData& cellIpData, double deltaT)
{
    return mLaw.Stress(StrainMechanics(cellIpData,deltaT), deltaT, cellIpData.Ids());
}

template <int TDim>
Eigen::VectorXd MultiPhysicsMomentumBalance<TDim>::StressShrinkage(const CellIpData& cellIpData, double deltaT)
{
    return mShrinkageIntegrand.Stress(cellIpData, deltaT);
}

template <int TDim>
DofVector<double> MultiPhysicsMomentumBalance<TDim>::Gradient(const NuTo::CellIpData& cellIpData, double deltaT)
{
    DofVector<double> gradient;

    Eigen::MatrixXd B = cellIpData.B(mDofDisp, Nabla::Strain());

    gradient[mDofDisp] = B.transpose() * Stress(cellIpData, deltaT);

    return gradient;
}

template <int TDim>
DofMatrix<double> MultiPhysicsMomentumBalance<TDim>::Hessian0(const CellIpData& cellIpData, double deltaT)
{
    DofMatrix<double> hessian0;

    BMatrixStrain B = cellIpData.B(mDofDisp, Nabla::Strain());
    hessian0(mDofDisp, mDofDisp) =
            B.transpose() * mLaw.Tangent(cellIpData.Apply(mDofDisp, Nabla::Strain()), deltaT, cellIpData.Ids()) * B;

    return hessian0;
}


} // namespace Integrands
} // namespace NuTo
