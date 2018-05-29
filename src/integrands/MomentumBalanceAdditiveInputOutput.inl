#pragma once

#include "integrands/MomentumBalanceAdditiveInputOutput.h"

namespace NuTo
{
namespace Integrands
{

template <int TDim, class TAdditiveIOIntegrand>
template <class TMechanicsLaw>
MomentumBalanceAdditiveInputOutput<TDim, TAdditiveIOIntegrand>::MomentumBalanceAdditiveInputOutput(
        DofType dofDisplacements, const TMechanicsLaw& mechanicsLaw, const TAdditiveIOIntegrand& aIOIntegrand)
    : mDofDisplacement{dofDisplacements}
    , mMechanicsLaw{std::make_unique<TMechanicsLaw>(mechanicsLaw)}
    , mAIOIntegrand{aIOIntegrand}
{
}

template <int TDim, class TAdditiveIOIntegrand>
DofVector<double> MomentumBalanceAdditiveInputOutput<TDim, TAdditiveIOIntegrand>::Gradient(const CellIpData& cellIpData,
                                                                                           double deltaT)
{
    DofVector<double> gradient;

    Eigen::MatrixXd B = cellIpData.B(mDofDisplacement, Nabla::Strain());

    gradient[mDofDisplacement] = B.transpose() * Stress(cellIpData, deltaT);

    return gradient;
}

template <int TDim, class TAdditiveIOIntegrand>
DofMatrix<double> MomentumBalanceAdditiveInputOutput<TDim, TAdditiveIOIntegrand>::Hessian0(const CellIpData& cellIpData,
                                                                                           double deltaT)
{
    // Current implementation is only for staggered approaches! The combinations dofDisp - dofAdditive must be
    // implemented to use it in a non staggered approach!

    DofMatrix<double> hessian0;

    Eigen::MatrixXd B = cellIpData.B(mDofDisplacement, Nabla::Strain());
    hessian0(mDofDisplacement, mDofDisplacement) =
            B.transpose() *
            mMechanicsLaw->Tangent(cellIpData.Apply(mDofDisplacement, Nabla::Strain()), deltaT, cellIpData.Ids()) * B;

    return hessian0;
}


template <int TDim, class TAdditiveIOIntegrand>
DofMatrix<double> MomentumBalanceAdditiveInputOutput<TDim, TAdditiveIOIntegrand>::Hessian1(const CellIpData& cellIpData,
                                                                                           double deltaT)
{
    throw Exception(__PRETTY_FUNCTION__,"Not implemented");
}

template <int TDim, class TAdditiveIOIntegrand>
DofMatrix<double> MomentumBalanceAdditiveInputOutput<TDim, TAdditiveIOIntegrand>::Hessian2(const CellIpData& cellIpData,
                                                                                           double deltaT)
{
    throw Exception(__PRETTY_FUNCTION__,"Not implemented");
}

template <int TDim, class TAdditiveIOIntegrand>
Eigen::VectorXd MomentumBalanceAdditiveInputOutput<TDim, TAdditiveIOIntegrand>::Strain(const CellIpData& cellIpData,
                                                                                       double deltaT)
{
    return cellIpData.Apply(mDofDisplacement, Nabla::Strain());
}

template <int TDim, class TAdditiveIOIntegrand>
Eigen::VectorXd
MomentumBalanceAdditiveInputOutput<TDim, TAdditiveIOIntegrand>::StrainMechanics(const CellIpData& cellIpData,
                                                                                double deltaT)
{
    return Strain(cellIpData, deltaT) - StrainAdditive(cellIpData, deltaT);
}

template <int TDim, class TAdditiveIOIntegrand>
Eigen::VectorXd
MomentumBalanceAdditiveInputOutput<TDim, TAdditiveIOIntegrand>::StrainAdditive(const CellIpData& cellIpData,
                                                                               double deltaT)
{
    return mAIOIntegrand.Strain(cellIpData, deltaT);
}

template <int TDim, class TAdditiveIOIntegrand>
Eigen::VectorXd MomentumBalanceAdditiveInputOutput<TDim, TAdditiveIOIntegrand>::Stress(const CellIpData& cellIpData,
                                                                                       double deltaT)
{
    return StressMechanics(cellIpData, deltaT) + StressAdditive(cellIpData, deltaT);
}

template <int TDim, class TAdditiveIOIntegrand>
Eigen::VectorXd
MomentumBalanceAdditiveInputOutput<TDim, TAdditiveIOIntegrand>::StressMechanics(const CellIpData& cellIpData,
                                                                                double deltaT)
{
    return mMechanicsLaw->Stress(StrainMechanics(cellIpData, deltaT), deltaT, cellIpData.Ids());
}

template <int TDim, class TAdditiveIOIntegrand>
Eigen::VectorXd
MomentumBalanceAdditiveInputOutput<TDim, TAdditiveIOIntegrand>::StressAdditive(const CellIpData& cellIpData,
                                                                               double deltaT)
{
    return mAIOIntegrand.Stress(cellIpData, deltaT);
}

} // namespace Integrands
} // namespace NuTo
