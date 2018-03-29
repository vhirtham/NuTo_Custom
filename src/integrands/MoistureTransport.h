#pragma once
#include "nuto/mechanics/dofs/DofVector.h"
#include "nuto/mechanics/dofs/DofMatrix.h"
#include "nuto/mechanics/cell/CellIpData.h"
#include "nuto/mechanics/cell/CellData.h"

namespace NuTo
{
namespace Integrands
{

template <int TDim, typename TDCw, typename TDCg, typename TMeC, typename TWVEq>
class MoistureTransport
{
    DofType mDofTypeWV;
    DofType mDofTypeRH;
    double mRho_w;
    double mRho_g_sat;
    double mPV;

public:
    //! @brief ctor
    inline MoistureTransport(DofType dofTypeWV, DofType dofTypeRH, double rho_w, double rho_g_sat, double PV);


    inline DofVector<double> Gradient(const CellIpData& cellIpData, double deltaT);

    inline DofMatrix<double> Stiffness(const CellIpData& cellIpData, double deltaT);

    inline DofMatrix<double> Damping(const CellIpData& cellIpData, double deltaT);

    void CheckValuesValid(double wv);
};
}
}

#include "MoistureTransport.inl"
