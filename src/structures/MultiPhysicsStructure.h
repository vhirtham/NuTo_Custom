#pragma once


#include "nuto/mechanics/cell/CellInterface.h"
#include "nuto/mechanics/constraints/Constraints.h"
#include "nuto/mechanics/dofs/DofInfo.h"
#include "nuto/mechanics/dofs/DofType.h"
#include "nuto/mechanics/dofs/DofVector.h"
#include "nuto/mechanics/dofs/DofMatrixSparse.h"
#include "nuto/mechanics/integrationtypes/IntegrationTypeTriangle.h"
#include "nuto/mechanics/integrationtypes/IntegrationTypeTensorProduct.h"
#include "nuto/mechanics/mesh/MeshFem.h"

// Integrands
#include "integrands/MoistureTransport.h"
#include "integrands/Shrinkage.h"
#include "integrands/MoistureTransportBoundary.h"
#include "integrands/MoistureTransportCoefficients.h"
#include "nuto/mechanics/constitutive/LinearElastic.h"
#include "integrands/MultiPhysicsMomentumBalance.h"

#include "boost/ptr_container/ptr_vector.hpp"


namespace NuTo
{

class GlobalDofVector;
class GlobalDofMatrixSparse;


using MoistureTransportBoundaryIntegrand =
        Integrands::MoistureTransportBoundary<2, MTCConst<1>, MTCConst<1, 100>, MTCConst<1, 10>>;

class MultiPhysicsStructure
{
    DofType mDofDisp = {"Displacements", 2};
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
    Laws::LinearElastic<2> mLawLinearElastic;
    Integrands::Shrinkage<2> mShrinkage;
    Integrands::MultiPhysicsMomentumBalance<2> mMPMomentumBalance;
    Integrands::MoistureTransport mMoistureTransport;
    MoistureTransportBoundaryIntegrand mMoistureTransportBoundary;


public:
    MultiPhysicsStructure(double rho_w = 1., double rho_g_sat = 0.1, double PV = 0.2)
        : mLawLinearElastic(30.e9, 0.0)
        , mShrinkage(mDofDisp, mDofWV, mDofRH, 0.5)
        , mMPMomentumBalance(mDofDisp, mDofWV, mDofRH, mLawLinearElastic, mShrinkage)
        , mMoistureTransport(mDofWV, mDofRH, MTCoefficientConstant{1.}, MTCoefficientConstant{0.01},
                             MTCoefficientConstant{0.}, MTCoefficientConstant{0.1}, rho_w, rho_g_sat, PV)
        , mMoistureTransportBoundary(mDofWV, mDofRH)
    {
    }

    Integrands::MultiPhysicsMomentumBalance<2>& GetMultiPhysicsLaw()
    {
        return mMPMomentumBalance;
    }

    std::vector<DofType> GetDofs()
    {
        return {mDofDisp, mDofWV, mDofRH};
    }

    const Group<CellInterface> GetMoistureTransportCells() const
    {
        return mGrpCellsMT;
    }

    Constraint::Constraints& GetConstraints()
    {
        return mConstraints;
    }

    MeshFem& GetMesh()
    {
        return mMesh;
    }

    DofVector<double> GradientMechanics();
    DofVector<double> GradientMoistureTransport();
    DofVector<double> GradientMoistureTransportBoundary();

    DofMatrixSparse<double> StiffnessMechanics();
    DofMatrixSparse<double> StiffnessMoistureTransport();
    DofMatrixSparse<double> StiffnessMoistureTransportBoundary();
    DofMatrixSparse<double> DampingMoistureTransport();

    void CheckDofNumbering();
    DofVector<double> CreateGlobalDofVector(std::vector<DofType> dofs);
    void GetDofVector(DofVector<double>& dofVector, std::vector<DofType> dofs, int instance = 0);

    void MergeDofs(Eigen::VectorXd values, std::vector<DofType> dofs, int instance = 0);
    void MergeDofVector(DofVector<double>& dofs, int instance = 0);

    void CreateUnitMesh(int numX, int numY, const InterpolationSimple& interpolationDispArg,
                        const InterpolationSimple& interpolationWVArg, const InterpolationSimple& interpolationRHArg,
                        const InterpolationSimple& interpolationMTBoundaryArg);

    void AddBoundaryElementsAtAxis(int& cellId, const InterpolationSimple& interpolationMTBoundary,
                                   eDirection direction, double axisOffset = 0.0);

    Eigen::VectorXd ExtractDofs(std::vector<DofType> dofs, int instance = 0);
};
}
