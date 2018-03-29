#pragma once

#include "nuto/mechanics/cell/CellIpData.h"
#include "nuto/mechanics/dofs/DofType.h"
#include "nuto/mechanics/dofs/DofMatrix.h"
#include "nuto/mechanics/dofs/DofVector.h"

namespace NuTo
{
namespace Integrands
{
template <int TDim, typename TDCw, typename TDCg, typename TWVEq>
class MoistureTransportBoundary
{
    double mExchangeRateWV = 2.;
    double mExchangeRateRH = 2.;
    DofType mDofTypeWV;
    DofType mDofTypeRH;

public:
    inline MoistureTransportBoundary(DofType dofTypeWV, DofType dofTypeRH)
        : mDofTypeWV(dofTypeWV)
        , mDofTypeRH(dofTypeRH){};

    inline DofVector<double> Gradient(const CellIpData& cellIpData, double deltaT);

    inline DofMatrix<double> Stiffness(const CellIpData& cellIpData, double deltaT);

    inline DofMatrix<double> Damping(const CellIpData& cellIpData, double deltaT);
};
}
}

#include "MoistureTransportBoundary.inl"
