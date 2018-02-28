#pragma once
#include "mechanics/dofs/DofVector.h"
#include "mechanics/dofs/DofMatrix.h"
#include "mechanics/cell/CellIpData.h"
#include "mechanics/cell/CellData.h"

namespace NuTo
{
namespace Integrands
{

template <int TDim, typename TDCw, typename TDCg, typename TMeC, typename TWVEq>
class MoistureTransport
{
    DofType mDofTypeRH;
    DofType mDofTypeWV;
    double mRho_w;
    double mRho_g_sat;
    double mPV;

public:
    //! @brief ctor
    MoistureTransport(DofType dofTypeRH, DofType dofTypeWV, double rho_w, double rho_g_sat, double PV);


    DofVector<double> Gradient(const CellIpData& cellIpData, double deltaT);

    DofMatrix<double> Stiffness(const CellIpData& cellIpData, double deltaT);

    DofMatrix<double> Damping(const CellIpData& cellIpData, double deltaT);
};
}
}

#include "MoistureTransport.inl"
