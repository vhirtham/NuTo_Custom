#pragma once

#include "nuto/mechanics/constitutive/MechanicsInterface.h"
#include "nuto/mechanics/dofs/DofType.h"
#include "nuto/mechanics/dofs/DofVector.h"
#include "nuto/mechanics/dofs/DofMatrix.h"
#include "nuto/mechanics/interpolation/TypeDefs.h"
#include "nuto/mechanics/cell/CellIpData.h"
#include "nuto/mechanics/cell/CellData.h"

#include "integrands/Shrinkage.h"

namespace NuTo
{
namespace Integrands
{

template <int TDim>
class MultiPhysicsMomentumBalance
{
    DofType mDofDisp;
    DofType mDofWV;
    DofType mDofRH;
    const Laws::MechanicsInterface<TDim>& mLaw;
    Shrinkage<TDim> mShrinkageIntegrand;

public:
    MultiPhysicsMomentumBalance(DofType dofDisp, DofType dofWV, DofType dofRH,
                                const Laws::MechanicsInterface<TDim>& law);

    Eigen::VectorXd Strain(const NuTo::CellIpData& cellIpData, double deltaT);
    Eigen::VectorXd StrainMechanics(const NuTo::CellIpData& cellIpData, double deltaT);
    Eigen::VectorXd StrainShrinkage(const NuTo::CellIpData& cellIpData, double deltaT);

    Eigen::VectorXd Stress(const NuTo::CellIpData& cellIpData, double deltaT);
    Eigen::VectorXd StressMechanics(const NuTo::CellIpData& cellIpData, double deltaT);
    Eigen::VectorXd StressShrinkage(const NuTo::CellIpData& cellIpData, double deltaT);

    DofVector<double> Gradient(const CellIpData& cellIpData, double deltaT);

    DofMatrix<double> Hessian0(const CellIpData& cellIpData, double deltaT);

protected:
};


} // namespace Integrands
} // namespace NuTo

#include "integrands/MultiPhysicsMomentumBalance.inl"
