#pragma once

#include "nuto/mechanics/cell/CellIpData.h"
#include "nuto/mechanics/dofs/DofType.h"
#include "nuto/mechanics/dofs/DofMatrix.h"
#include "nuto/mechanics/dofs/DofVector.h"

namespace NuTo
{
namespace Integrands
{
class Shrinkage
{
public:

//! @brief ctor
inline Shrinkage(DofType dofTypeDisp, DofType dofTypeWV, DofType dofTypeRH);


inline DofVector<double> Gradient(const CellIpData& cellIpData, double deltaT);

inline DofMatrix<double> Stiffness(const CellIpData& cellIpData, double deltaT);

inline DofMatrix<double> Damping(const CellIpData& cellIpData, double deltaT);

};
}
}
