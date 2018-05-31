#pragma once
#include "nuto/mechanics/dofs/DofVector.h"
#include "nuto/mechanics/dofs/DofMatrix.h"
#include "nuto/mechanics/cell/CellIpData.h"
#include "nuto/mechanics/cell/CellData.h"


namespace NuTo
{
class MTCoefficientInterface;

namespace Integrands
{


class MoistureTransport
{
    friend class MoistureTransportBoundary;


    DofType mDofTypeWV;
    DofType mDofTypeRH;
    double mRho_w;
    double mRho_g_sat;
    double mPV;

    std::unique_ptr<MTCoefficientInterface> mDCw = nullptr;
    std::unique_ptr<MTCoefficientInterface> mDCg = nullptr;
    std::unique_ptr<MTCoefficientInterface> mMeC = nullptr;
    std::unique_ptr<MTCoefficientInterface> mWVEq = nullptr;

public:
    MoistureTransport(MoistureTransport&&) = default;

    //! @brief ctor
    MoistureTransport(DofType dofTypeWV, DofType dofTypeRH, const MTCoefficientInterface& diffCoeffWV,
                      const MTCoefficientInterface& diffCoeffRH, const MTCoefficientInterface& massExchCoeff,
                      const MTCoefficientInterface& wvEquilibriumCoeff, double rho_w, double rho_g_sat, double PV);


    DofVector<double> Gradient(const CellIpData& cellIpData, double deltaT);

    DofMatrix<double> Stiffness(const CellIpData& cellIpData, double deltaT);

    DofMatrix<double> Damping(const CellIpData& cellIpData, double deltaT);

    DofMatrix<double> Mass(const CellIpData& cellIpData, double deltaT);

    void CheckValuesValid(double wv);
};
}
}
