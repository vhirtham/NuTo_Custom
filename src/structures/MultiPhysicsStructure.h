#pragma once


#include "nuto/mechanics/cell/CellInterface.h"
#include "nuto/mechanics/constraints/Constraints.h"
#include "nuto/mechanics/dofs/DofInfo.h"
#include "nuto/mechanics/dofs/DofType.h"
//#include "nuto/mechanics/dofs/GlobalDofVector.h"
//#include "nuto/mechanics/dofs/GlobalDofMatrixSparse.h"
#include "nuto/mechanics/integrationtypes/IntegrationTypeTriangle.h"
#include "nuto/mechanics/integrationtypes/IntegrationTypeTensorProduct.h"
#include "nuto/mechanics/mesh/MeshFem.h"

#include "integrands/MoistureTransport.h"
#include "integrands/MoistureTransportBoundary.h"
#include "integrands/MoistureTransportCoefficients.h"

#include "boost/ptr_container/ptr_vector.hpp"


namespace NuTo {

class GlobalDofVector;
class GlobalDofMatrixSparse;

using MoistureTransportIntegrand = Integrands::MoistureTransport<2, MTCConst<1>, MTCConst<1, 100>, MTCConst<0>, MTCConst<1, 10>>;
using MoistureTransportBoundaryIntegrand = Integrands::MoistureTransportBoundary<2, MTCConst<1>, MTCConst<1, 100>, MTCConst<1, 10>>;

class MultiPhysicsStructure
{
    DofType mDofDisp = {"Displacements", 3};
    DofType mDofWV = {"waterVolumeFraction", 1};
    DofType mDofRH = {"relativeHumidity", 1};
    DofInfo mDofInfo;
    Constraint::Constraints mConstraints;
    IntegrationTypeTriangle integrationType = {2};
    IntegrationTypeTensorProduct<1> mIntegrationTypeBoundary = {3, eIntegrationMethod::GAUSS};
    MeshFem mMesh;
    Group<CellInterface> mGrpCellsMT;
    Group<CellInterface> mGrpCellsMTBoundary;
    boost::ptr_vector<CellInterface> mCellContainer;
    MoistureTransportIntegrand mMoistureTransport;
    MoistureTransportBoundaryIntegrand mMoistureTransportBoundary;

public:
    MultiPhysicsStructure(double rho_w = 1., double rho_g_sat = 0.1, double PV = 0.2)
        : mMoistureTransport(mDofWV, mDofRH, rho_w, rho_g_sat, PV)
        , mMoistureTransportBoundary(mDofWV, mDofRH)
    {
    }

    std::vector<DofType> GetDofs()
    {
        return {mDofDisp, mDofWV, mDofRH};
    }

    const Group<CellInterface> GetMoistureTransportCells() const
    {
        return mGrpCellsMT;
    }

    MeshFem& GetMesh()
    {
        return mMesh;
    }


    GlobalDofVector GradientMoistureTransport();
    GlobalDofVector GradientMoistureTransportBoundary();

    GlobalDofMatrixSparse StiffnessMoistureTransport();
    GlobalDofMatrixSparse StiffnessMoistureTransportBoundary();
    GlobalDofMatrixSparse DampingMoistureTransport();

    void CheckDofNumbering();
    GlobalDofVector CreateGlobalDofVector(std::vector<DofType> dofs);
    void GetDofVector(GlobalDofVector& dofVector, std::vector<DofType> dofs, int instance = 0);

    void MergeDofs(Eigen::VectorXd values, std::vector<DofType> dofs, int instance = 0);
    void MergeDofVector(GlobalDofVector& dofs, int instance = 0);

    void CreateUnitMesh(int numX, int numY, const InterpolationSimple& interpolationDispArg,
                        const InterpolationSimple& interpolationWVArg, const InterpolationSimple& interpolationRHArg, const InterpolationSimple& interpolationMTBoundaryArg);

    void AddBoundaryElementsAtAxis(int& cellId, const InterpolationSimple& interpolationMTBoundary, eDirection direction, double axisOffset = 0.0);

    Eigen::VectorXd ExtractDofs(std::vector<DofType> dofs, int instance = 0);
};

}
