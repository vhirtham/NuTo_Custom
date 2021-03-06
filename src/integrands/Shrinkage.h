#pragma once

#include "nuto/mechanics/cell/CellIpData.h"
#include "nuto/mechanics/dofs/DofType.h"
#include "nuto/mechanics/dofs/DofMatrix.h"
#include "nuto/mechanics/dofs/DofVector.h"

namespace NuTo
{
namespace Integrands
{
template <int TDim>
class Shrinkage
{
    DofType mDofDisp;
    DofType mDofWV;
    DofType mDofRH;

    double mWeightingFactor = 0.0;
    double mStressScalingFactor = 1.0;

public:
    //! @brief ctor
    inline Shrinkage(DofType dofTypeDisp, DofType dofTypeWV, DofType dofTypeRH, double weightingFactor,
                     double stressScalingFactor = 1.0);


    inline Eigen::VectorXd Stress(const CellIpData& cellIpData, double deltaT) const;
    inline Eigen::VectorXd Strain(const CellIpData& cellIpData, double deltaT) const;


private:
    inline double CapillaryPressure(double rh) const;

    inline Eigen::VectorXd ComponentVector() const;
};
}
}

#include "Shrinkage.inl"
