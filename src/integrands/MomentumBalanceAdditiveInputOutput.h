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

template <int TDim, class TAdditiveIOIntegrand>
class MomentumBalanceAdditiveInputOutput
{
    DofType mDofDisplacement;
    std::unique_ptr<Laws::MechanicsInterface<TDim>> mMechanicsLaw;
    TAdditiveIOIntegrand mAIOIntegrand;

public:
    template <class TMechanicsLaw>
    MomentumBalanceAdditiveInputOutput(DofType dofDisplacements, const TMechanicsLaw& mechanicsLaw,
                                       const TAdditiveIOIntegrand& aIOIntegrand);

    DofVector<double> Gradient(const CellIpData& cellIpData, double deltaT);

    //! @remark Current implementation is only for staggered approaches! The combinations dofDisp - dofAdditive must be
    //! implemented to use it in a non staggered approach!
    inline DofMatrix<double> Hessian0(const CellIpData& cellIpData, double deltaT);
    inline DofMatrix<double> Hessian1(const CellIpData& cellIpData, double deltaT);
    inline DofMatrix<double> Hessian2(const CellIpData& cellIpData, double deltaT);


    Eigen::VectorXd Strain(const NuTo::CellIpData& cellIpData, double deltaT) const;
    Eigen::VectorXd StrainMechanics(const NuTo::CellIpData& cellIpData, double deltaT) const;
    Eigen::VectorXd StrainAdditive(const NuTo::CellIpData& cellIpData, double deltaT) const;

    Eigen::VectorXd Stress(const NuTo::CellIpData& cellIpData, double deltaT) const;
    Eigen::VectorXd StressMechanics(const NuTo::CellIpData& cellIpData, double deltaT) const;
    Eigen::VectorXd StressAdditive(const NuTo::CellIpData& cellIpData, double deltaT) const;
};
} // namespace Integrands
} // namespace NuTo

#include "integrands/MomentumBalanceAdditiveInputOutput.inl"
