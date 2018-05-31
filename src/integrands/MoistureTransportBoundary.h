#pragma once

#include "nuto/mechanics/cell/CellIpData.h"
#include "nuto/mechanics/dofs/DofType.h"
#include "nuto/mechanics/dofs/DofMatrix.h"
#include "nuto/mechanics/dofs/DofVector.h"

#include "integrands/MoistureTransportCoefficients.h"

namespace NuTo
{
namespace Integrands
{
class MoistureTransport;

class MoistureTransportBoundary
{
    const MoistureTransport& mMoistureTransportIntegrand;
    DofType mDofTypeWV;
    DofType mDofTypeRH;
    double mExchangeRateWV = 0.;
    double mExchangeRateRH = 0.;
    double mEnvRH = 1.0;

public:
    MoistureTransportBoundary(MoistureTransportBoundary&&) = default;

    inline MoistureTransportBoundary(const MoistureTransport& moistureTransportIntegrand, DofType dofTypeWV,
                                     DofType dofTypeRH, double massExchangeRateWV, double massExchangeRateRH,
                                     double initialEnvironmentalRH);

    inline DofVector<double> Gradient(const CellIpData& cellIpData, double deltaT);

    inline DofMatrix<double> Stiffness(const CellIpData& cellIpData, double deltaT);

    void SetEnvironmentalRelativeHumidity(double envRH);
};
}
}

#include "MoistureTransportBoundary.inl"
