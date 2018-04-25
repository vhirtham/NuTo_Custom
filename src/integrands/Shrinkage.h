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

public:
    //! @brief ctor
    inline Shrinkage(DofType dofTypeDisp, DofType dofTypeWV, DofType dofTypeRH, double weightingFactor);


    inline Eigen::VectorXd Stress(const CellIpData& cellIpData, double deltaT);
    inline Eigen::VectorXd Strain(const CellIpData& cellIpData, double deltaT);


private:
    inline double CapillaryPressure(double rh);

    inline Eigen::VectorXd ComponentVector();
};
}
}

#include "Shrinkage.inl"
