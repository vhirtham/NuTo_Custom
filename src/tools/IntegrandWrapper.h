#pragma once

#include "nuto/mechanics/cell/CellIpData.h"

#include "nuto/mechanics/dofs/DofMatrix.h"
#include "nuto/mechanics/dofs/DofVector.h"

namespace NuTo
{

class IntegrandWrapperInterface
{
public:
    virtual DofVector<double> Residual(const CellIpData& cellIpData, double deltaT) = 0;
};


template <class TIntegrand>
class IntegrandWrapper : public IntegrandWrapperInterface
{
    typedef DofVector<double> (TIntegrand::*DofVectorFunction)(const CellIpData&, double);
    typedef DofMatrix<double> (TIntegrand::*DofMatrixFunction)(const CellIpData&, double);

    TIntegrand mIntegrand;
    DofVectorFunction mGradientFunction;

public:
    template <typename... Args>
    IntegrandWrapper(DofVector<double> (TIntegrand::*gradientFuncion)(const CellIpData&, double),
                     Args... args)
        : mIntegrand{std::forward<Args>(args)...}
        , mGradientFunction{gradientFuncion}
    {
    }

    virtual DofVector<double> Residual(const CellIpData& cellIpData, double deltaT) override
    {
        return (mIntegrand.*mGradientFunction)(cellIpData, deltaT);
    }
};
}
