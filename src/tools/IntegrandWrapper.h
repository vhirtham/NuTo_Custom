#pragma once

#include "nuto/mechanics/cell/CellIpData.h"

#include "nuto/mechanics/dofs/DofMatrix.h"
#include "nuto/mechanics/dofs/DofVector.h"

namespace NuTo
{

template <typename... TAdditionalParams>
class IntegrandWrapperInterface
{
public:
    virtual DofVector<double> Residual(const CellIpData&, TAdditionalParams...) = 0;
    virtual DofMatrix<double> Stiffness(const CellIpData&, TAdditionalParams...) = 0;
    virtual DofMatrix<double> Damping(const CellIpData&, TAdditionalParams...) = 0;
    virtual DofMatrix<double> Mass(const CellIpData&, TAdditionalParams...) = 0;
};


template <class TIntegrand, typename... TAdditionalParams>
class IntegrandWrapper : public IntegrandWrapperInterface<TAdditionalParams...>
{
    typedef DofVector<double> (TIntegrand::*DofVectorFunction)(const CellIpData&, TAdditionalParams...);
    typedef DofMatrix<double> (TIntegrand::*DofMatrixFunction)(const CellIpData&, TAdditionalParams...);

    TIntegrand mIntegrand;
    DofVectorFunction mResidualFunction = nullptr;
    DofMatrixFunction mStiffnessFunction = nullptr;
    DofMatrixFunction mDampingFunction = nullptr;
    DofMatrixFunction mMassFunction = nullptr;

public:
    IntegrandWrapper(TIntegrand&& integrand, DofVectorFunction gradientFuncion, DofMatrixFunction stiffnessFunction,
                     DofMatrixFunction dampingFunction = nullptr, DofMatrixFunction massFunction = nullptr)
        : mIntegrand{std::move(integrand)}
        , mResidualFunction{gradientFuncion}
        , mStiffnessFunction{stiffnessFunction}
        , mDampingFunction{dampingFunction}
        , mMassFunction{massFunction}
    {
    }

    template <typename... Args>
    IntegrandWrapper(DofVectorFunction gradientFuncion, DofMatrixFunction stiffnessFunction,
                     DofMatrixFunction dampingFunction = nullptr, DofMatrixFunction massFunction = nullptr,
                     Args... args)
        : mIntegrand{std::forward<Args>(args)...}
        , mResidualFunction{gradientFuncion}
        , mStiffnessFunction{stiffnessFunction}
        , mDampingFunction{dampingFunction}
        , mMassFunction{massFunction}
    {
    }

    virtual DofVector<double> Residual(const CellIpData& cellIpData, TAdditionalParams... additionalParams) override
    {
        return WrappedFunctionCall<DofVector<double>>(mResidualFunction, "residual", cellIpData,
                                                      std::forward<TAdditionalParams>(additionalParams)...);
    }

    virtual DofMatrix<double> Stiffness(const CellIpData& cellIpData, TAdditionalParams... additionalParams) override
    {
        return WrappedFunctionCall<DofMatrix<double>>(mStiffnessFunction, "stiffness", cellIpData,
                                                      std::forward<TAdditionalParams>(additionalParams)...);
    }

    virtual DofMatrix<double> Damping(const CellIpData& cellIpData, TAdditionalParams... additionalParams) override
    {
        return WrappedFunctionCall<DofMatrix<double>>(mDampingFunction, "damping", cellIpData,
                                                      std::forward<TAdditionalParams>(additionalParams)...);
    }

    virtual DofMatrix<double> Mass(const CellIpData& cellIpData, TAdditionalParams... additionalParams) override
    {
        return WrappedFunctionCall<DofMatrix<double>>(mMassFunction, "mass", cellIpData,
                                                      std::forward<TAdditionalParams>(additionalParams)...);
    }

    TIntegrand& GetIntegrand()
    {
        return mIntegrand;
    }

private:
    template <typename TFuncPtr>
    void CheckNullPtr(TFuncPtr funcPtr, std::string functionName)
    {
        if (funcPtr == nullptr)
            throw Exception(__PRETTY_FUNCTION__,
                            "No " + functionName + " function wrapped. Function pointer is nullptr.");
    }

    template <typename TReturnType, typename TFuncPtr>
    TReturnType WrappedFunctionCall(TFuncPtr funcPtr, std::string functionName, const CellIpData& cellIpData,
                                    TAdditionalParams... additionalParams)
    {
        CheckNullPtr(funcPtr, functionName);
        return (mIntegrand.*funcPtr)(cellIpData, std::forward<TAdditionalParams>(additionalParams)...);
    }
};
}
